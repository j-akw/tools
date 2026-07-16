#!/usr/bin/env python3
"""
update_ind_pop.py — update the population column of an EIGENSTRAT .ind file.

Given a two-column mapping file (sample_id  new_pop), rewrite the .ind file
so that matching samples get their pop column replaced. Non-matching samples
are passed through unchanged.

Usage:
    python update_ind_pop.py -i input.ind -m pop_map.txt -o output.ind

    # In-place (overwrites input, backup written to <input>.bak):
    python update_ind_pop.py -i input.ind -m pop_map.txt --in-place
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import Dict


def load_pop_map(path: Path) -> Dict[str, str]:
    """Load a two-column (id, pop) mapping file. Whitespace-separated."""
    m: Dict[str, str] = {}
    with path.open() as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.split()
            if not parts:
                continue
            if len(parts) < 2:
                print(
                    f"warning: map line {lineno} has fewer than 2 columns, skipping: {line.rstrip()}",
                    file=sys.stderr,
                )
                continue
            m[parts[0]] = parts[1]
    return m


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--ind", type=Path, required=True, help="input .ind file")
    p.add_argument("-m", "--map", type=Path, required=True,
                   help="two-column file: sample_id new_pop")
    p.add_argument("-o", "--out", type=Path, help="output .ind (default: stdout unless --in-place)")
    p.add_argument("--in-place", action="store_true",
                   help="overwrite input; writes a .bak backup first")
    args = p.parse_args()

    if args.in_place and args.out:
        sys.exit("error: --in-place and -o are mutually exclusive")

    pop_map = load_pop_map(args.map)
    print(f"loaded {len(pop_map)} id->pop mappings", file=sys.stderr)

    # Read all input lines first so we can safely overwrite in-place.
    with args.ind.open() as fh:
        lines = fh.readlines()

    updated_lines = []
    updated_ids = set()
    unchanged_count = 0

    for line in lines:
        parts = line.split()
        if len(parts) < 3:
            # Malformed or blank — pass through as-is.
            updated_lines.append(line)
            continue
        sid, sex, pop = parts[0], parts[1], parts[2]
        new_pop = pop_map.get(sid)
        if new_pop is not None and new_pop != pop:
            pop = new_pop
            updated_ids.add(sid)
        else:
            unchanged_count += 1
        updated_lines.append(f"{sid:>20} {sex:1} {pop:>10}\n")

    # Write output
    if args.in_place:
        backup = args.ind.with_suffix(args.ind.suffix + ".bak")
        shutil.copy2(args.ind, backup)
        print(f"backup written to {backup}", file=sys.stderr)
        with args.ind.open("w") as fh:
            fh.writelines(updated_lines)
        dest = args.ind
    elif args.out:
        with args.out.open("w") as fh:
            fh.writelines(updated_lines)
        dest = args.out
    else:
        sys.stdout.writelines(updated_lines)
        dest = Path("<stdout>")

    # Report
    map_ids_not_in_ind = set(pop_map.keys()) - updated_ids - {
        p.split()[0] for p in updated_lines if len(p.split()) >= 3
        and pop_map.get(p.split()[0]) == p.split()[2]  # already correct
    }
    # Simpler: any map entry whose id never appeared in the .ind
    ind_ids = {p.split()[0] for p in updated_lines if len(p.split()) >= 3}
    orphan_map_ids = sorted(set(pop_map.keys()) - ind_ids)

    print(f"updated {len(updated_ids)} samples", file=sys.stderr)
    print(f"left unchanged: {unchanged_count}", file=sys.stderr)
    if orphan_map_ids:
        print(
            f"warning: {len(orphan_map_ids)} IDs in map file not found in .ind "
            f"(first 10: {orphan_map_ids[:10]})",
            file=sys.stderr,
        )
    print(f"wrote {dest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
