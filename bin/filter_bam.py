#!/usr/bin/env python
"""Neutralize pairing flags for orphaned reads in a coordinate-sorted BAM.

Background
----------
blat_filter.py calls ``pysam.AlignmentFile.mate(read)`` for every read whose
flag claims it is in a proper pair (``read.alignment.is_proper_pair``). pysam's
``mate()`` performs a real coordinate-seek for the partner record and raises
``ValueError: mate not found`` if that record is not physically present in the
BAM. This crashes the whole job (see jobs/*/stderr_2):

    File "/app/blat_filter.py", line 709, in collect_mate_blat_queries
        mate = self._input_bam.mate(read.alignment)
    ValueError: mate not found

This happens when a read keeps the paired / proper-pair flags written at
alignment time, but its mate was physically removed afterwards (e.g. reads in
high-coverage regions dropped on purpose). The BAM is then internally
inconsistent.

What this script does
---------------------
Rather than dropping the orphaned reads (which would lose coverage used for
scoring), it *clears their pairing-related flags* so that downstream
``is_proper_pair`` / ``is_paired`` checks are False and ``mate()`` is never
called for them. The read sequence/position is preserved untouched.

A read's mate "exists" iff a record with the same query_name and the opposite
read1/read2 designation is present in the BAM. We determine this with two
streaming passes over the coordinate-sorted BAM:

  Pass 1  build the set of (query_name, is_read1) PRIMARY keys actually present
  Pass 2  for every paired read (primary, secondary or supplementary) whose
          mate key is absent, clear pairing flags

Why pass 1 only catalogues primaries, but pass 2 fixes every record
-------------------------------------------------------------------
``pysam.mate()`` only ever returns the *primary* mate alignment, so a read's
mate "exists" iff the opposite-member *primary* record is present. That is what
pass 1 records. But ``blat_filter.py`` inspects ``is_proper_pair`` on whatever
records the pileup yields, which can include secondary / supplementary
alignments. If one of those is orphaned and we leave its proper-pair flag set,
``mate()`` is still called and still raises. So pass 2 must neutralize the
pairing flags on *all* record types (primary, secondary, supplementary), each
judged against the primary present-set.

Flags cleared for an orphan (bit values):
  0x1   paired
  0x2   proper_pair
  0x8   mate_unmapped       (set, so anything that does inspect it is consistent)
  0x20  mate_reverse_strand (cleared)
The read is left as an effectively single-end record. We also reset the
mate-pointer fields (next_reference_id / next_reference_start / template_length)
so no stale coordinates remain.

Usage
-----
    python fix_orphan_mate_flags.py --in_bam in.bam --out_bam out.bam [--index]

The input must be coordinate-sorted (the standard recal BAM is). Output is
written in the same order; pass --index to also write a .bai.
"""

import argparse
import sys

import pysam


def mate_key(read):
    """Key identifying the *mate* this read points to: (name, mate_is_read1)."""
    # A proper read1's mate is read2 and vice-versa. is_read1/is_read2 are
    # derived from flags 0x40 / 0x80. If neither is set (rare, malformed),
    # fall back to treating the mate as "the other one" via is_read1 negation.
    return (read.query_name, not read.is_read1)


def present_key(read):
    """Key identifying this read's own identity: (name, is_read1)."""
    return (read.query_name, read.is_read1)


def build_present_set(bam_path):
    present = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            # Only primary, mapped, paired reads can be a "mate" that
            # mate() would look up. Secondary/supplementary records are not
            # what mate() returns, so don't let them mask a missing primary.
            if read.is_secondary or read.is_supplementary:
                continue
            if not read.is_paired:
                continue
            present.add(present_key(read))
    return present


def clear_pairing(read):
    """Turn a paired read into an effectively single-end record."""
    read.is_paired = False
    read.is_proper_pair = False
    read.mate_is_unmapped = True
    read.mate_is_reverse = False
    # Drop stale mate pointers.
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.template_length = 0
    return read


def main():
    ap = argparse.ArgumentParser(
        description="Clear pairing flags for orphaned reads so blat_filter.py's "
                    "mate() lookup never fails."
    )
    ap.add_argument("--in_bam", required=True, help="coordinate-sorted input BAM")
    ap.add_argument("--out_bam", required=True, help="output BAM path")
    ap.add_argument("--index", action="store_true",
                    help="also write a .bai index for the output BAM")
    args = ap.parse_args()

    sys.stderr.write("[fix_orphan_mate_flags] pass 1: cataloguing present reads\n")
    present = build_present_set(args.in_bam)
    sys.stderr.write("[fix_orphan_mate_flags] {0} paired primary records present\n".format(len(present)))

    n_total = 0
    n_fixed = 0
    sys.stderr.write("[fix_orphan_mate_flags] pass 2: writing filtered BAM\n")
    with pysam.AlignmentFile(args.in_bam, "rb") as in_bam:
        with pysam.AlignmentFile(args.out_bam, "wb", template=in_bam) as out_bam:
            for read in in_bam.fetch(until_eof=True):
                n_total += 1
                # Neutralize any paired record (primary, secondary OR
                # supplementary) whose primary mate is not present. mate() only
                # returns the primary mate, so presence is judged against the
                # primary present-set for every record type.
                if read.is_paired and mate_key(read) not in present:
                    clear_pairing(read)
                    n_fixed += 1
                out_bam.write(read)

    sys.stderr.write(
        "[fix_orphan_mate_flags] done: {0} records written, {1} orphans neutralized\n".format(
            n_total, n_fixed)
    )

    if args.index:
        sys.stderr.write("[fix_orphan_mate_flags] indexing {0}\n".format(args.out_bam))
        pysam.index(args.out_bam)


if __name__ == "__main__":
    main()
