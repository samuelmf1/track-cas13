#!/usr/bin/env python3
# offtarget_filter_cigar.py
import sys
import csv
import re
import argparse
from typing import List, Tuple, Dict

CIGAR_RE = re.compile(r"(\d*)([=XID])")


def parse_cigar(cigar: str) -> Tuple[Dict[str, int], List[Tuple[str, int, int]]]:
    """
    Parse a CIGAR string into:
      - counts of '=', 'X', 'I', 'D'
      - list of tuples (op, start_pos, end_pos), reference positions 1-indexed
    """
    counts = {"=": 0, "X": 0, "I": 0, "D": 0}
    positions: List[Tuple[str, int, int]] = []
    ref_pos = 1

    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str) if length_str else 1
        start = ref_pos
        end = ref_pos + (length - 1 if op in ("=", "X", "D") else 0)
        positions.append((op, start, end))

        if op in ("=", "X", "D"):
            ref_pos += length

        counts[op] += length

    return counts, positions


def valid_indel_positions(positions: List[Tuple[str, int, int]], debug: bool = False) -> Tuple[bool, str]:
    """
    Returns True if all insertions/deletions occur only at positions <=3 or >=21,
    allowing up to 2 consecutive indels of the same type.
    """
    i = 0
    n = len(positions)

    while i < n:
        op, start, _ = positions[i]
        if op in ("I", "D"):
            run_positions = [start]
            j = i + 1
            while j < n and positions[j][0] == op:
                run_positions.append(positions[j][1])
                j += 1

            # More than 2 consecutive I/D is invalid
            if len(run_positions) > 2:
                reason = f"{op} run too long at positions {run_positions}"
                return False, reason if debug else ""

            # All positions in run must be <=3 or >=21
            if any(pos > 3 and pos < 21 for pos in run_positions):
                reason = f"{op} at invalid position(s) {run_positions}"
                return False, reason if debug else ""

            i = j
        else:
            i += 1

    return True, ""


def log(debug: bool, *args):
    if debug:
        print(*args, file=sys.stderr)


def process_rows(reader, out, debug: bool = False):
    passed = failed = skipped = 0

    for row in reader:
        if not row or row[0].startswith("pat_id"):
            continue

        cigar = row[7]
        counts, positions = parse_cigar(cigar)
        total_ref_length = sum(counts.values()) - counts["I"]

        # Skip if reference length != 23
        if total_ref_length != 23:
            skipped += 1
            log(debug, f"[SKIP] {cigar}: total reference length != 23")
            continue

        x, i_count, d_count = counts["X"], counts["I"], counts["D"]
        if x + i_count + d_count > 2:
            failed += 1
            log(debug, f"[FAIL] {cigar}: too many edits (X={x}, I={i_count}, D={d_count})")
            continue

        # Rule 1: up to 2 mismatches, no indels
        if x <= 2 and i_count == 0 and d_count == 0:
            reason = f"{x} mismatch(es), no indels"
            passed += 1

        # Rule 2: any I/D, no mismatches, must be valid by positions
        elif x == 0 and (i_count > 0 or d_count > 0):
            ok, reason = valid_indel_positions(positions, debug)
            if ok:
                passed += 1
            else:
                failed += 1
                log(debug, f"[FAIL] {cigar}: {reason}")
                continue
        else:
            failed += 1
            log(debug, f"[FAIL] {cigar}: X={x}, I={i_count}, D={d_count} (invalid combination)")
            continue

        if passed:
            log(debug, f"[PASS] {cigar}: {reason}")
            out.write("\t".join(row) + "\n")

    if debug:
        log(debug, f"\n--- SUMMARY ---\nPassed: {passed}\nFailed: {failed}\nSkipped: {skipped}")


def main():
    parser = argparse.ArgumentParser(description="Filter SASSy results by CIGAR pattern and indel position.")
    parser.add_argument("input", nargs="?", default="-", help="Input TSV file (default: stdin)")
    parser.add_argument("--debug", action="store_true", help="Print debug information")
    args = parser.parse_args()

    reader = csv.reader(sys.stdin if args.input == "-" else open(args.input), delimiter="\t")
    process_rows(reader, sys.stdout, debug=args.debug)


if __name__ == "__main__":
    main()
