#!/usr/bin/env python3

"""
A tool to check two different EigenStrat databases for shared individuals,
and extract or remove individuals from an EigenStrat database.
"""

import sys
import argparse

VERSION = "2.0.0"


def file_len(fname):
    """Return the number of lines in a file using buffered counting."""
    count = 0
    with open(fname, "rb") as f:
        for _ in f:
            count += 1
    return count


def file_width(fname):
    """Return the number of genotypes per line in a .geno file."""
    with open(fname) as f:
        return len(f.readline().strip())


def validate_eigenstrat(genof, snpf, indf):
    """Check the consistency of an EigenStrat database."""
    geno_lines = file_len(genof)
    geno_width = file_width(genof)
    lines_snp = file_len(snpf)
    lines_ind = file_len(indf)

    if geno_lines != lines_snp:
        raise IOError("Input .snp and .geno files do not match.")
    if geno_width != lines_ind:
        raise IOError("Input .ind and .geno files do not match.")


def parse_ind_file(ind_fn):
    """Read the .ind file once and return a list of (name, sex, pop) tuples."""
    entries = []
    with open(ind_fn) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 3:
                entries.append((fields[0], fields[1], fields[2]))
    return entries


def check_duplicates(ind_fn, check_file):
    """Check two .ind files for duplicate individuals."""
    inds1 = set()
    with open(ind_fn) as f:
        for line in f:
            fields = line.strip().split()
            if fields:
                inds1.add(fields[0])

    inds2 = set()
    for line in check_file:
        fields = line.strip().split()
        if fields:
            inds2.add(fields[0])

    duplicates = inds1 & inds2
    for ind in sorted(duplicates):
        print("{:25s}".format(ind), "<---", "Duplicate individual", sep="\t", file=sys.stdout)
    print(
        "#Duplicate individual check finished.",
        len(duplicates),
        "duplicate individuals found.",
        sep=" ",
        file=sys.stdout,
    )


def build_sample_index(samples, ind_entries):
    """
    Build index, sex, and pop dicts for requested samples
    by looking them up in the pre-parsed ind entries.
    """
    # Build a lookup from the ind file (preserves first occurrence)
    ind_lookup = {}
    for idx, (name, sex, pop) in enumerate(ind_entries):
        if name not in ind_lookup:
            ind_lookup[name] = (idx, sex, pop)

    index = {}
    sex = {}
    pop = {}
    for sample in samples:
        if sample in ind_lookup:
            idx, s, p = ind_lookup[sample]
            index[sample] = idx
            sex[sample] = s
            pop[sample] = p

    return index, sex, pop


def extract_individuals(geno_fn, snp_fn, output_prefix, samples, index, sex, pop):
    """Extract selected individuals from the EigenStrat database."""
    # Pre-compute the ordered column indices for the samples that were found
    ordered = [(s, index[s]) for s in samples if s in index]

    with open(geno_fn) as geno_in, \
         open(output_prefix + ".geno", "w") as geno_out, \
         open(output_prefix + ".ind", "w") as ind_out, \
         open(output_prefix + ".snp", "w") as snp_out:

        # Write geno: one line at a time, picking only the columns we want
        for line in geno_in:
            fields = line.rstrip("\n")
            geno_out.write("".join(fields[idx] for _, idx in ordered) + "\n")

        # Write ind
        for s, _ in ordered:
            ind_out.write("{}\t{}\t{}\n".format(s, sex[s], pop[s]))

        # Copy snp unchanged
        with open(snp_fn) as snp_in:
            for line in snp_in:
                snp_out.write(line)

    # Report individuals requested but not found
    for s in samples:
        if s not in index:
            print(
                'Individual "', s, '" not found in ind file.',
                sep="",
                file=sys.stderr,
            )

    print("Extraction of ", len(ordered), " individuals complete.", sep="", file=sys.stderr)


def remove_individuals(geno_fn, snp_fn, ind_fn, output_prefix, samples, index):
    """Remove selected individuals from the EigenStrat database."""
    skip = set(index.values())
    sample_set = set(samples)

    with open(geno_fn) as geno_in, \
         open(output_prefix + ".geno", "w") as geno_out, \
         open(output_prefix + ".ind", "w") as ind_out, \
         open(output_prefix + ".snp", "w") as snp_out:

        # Write geno: skip columns at the indices to remove
        for line in geno_in:
            fields = line.rstrip("\n")
            geno_out.write("".join(c for i, c in enumerate(fields) if i not in skip) + "\n")

        # Write ind: skip matching sample names (re-read the file)
        with open(ind_fn) as ind_in:
            for line in ind_in:
                fields = line.strip().split()
                if fields and fields[0] not in sample_set:
                    ind_out.write(line)

        # Copy snp unchanged
        with open(snp_fn) as snp_in:
            for line in snp_in:
                snp_out.write(line)

    print("Removal of ", len(index), " individuals complete.", sep="", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s -g <GENO> -s <SNP> -i <IND> (-C <ind file> | -R | -E) "
              "[-L <SAMPLE LIST> | -S Ind [-S Ind2]] [-o <OUTPUT PREFIX>]",
        description="A tool to check two different EigenStrat databases for shared "
                    "individuals, and extract or remove individuals from an EigenStrat database.",
    )
    parser._optionals.title = "Available options"
    parser.add_argument("-g", "--genoFn", type=str, metavar="<GENO FILE>", required=True,
                        help="Path to the input .geno file.")
    parser.add_argument("-s", "--snpFn", type=str, metavar="<SNP FILE>", required=True,
                        help="Path to the input .snp file.")
    parser.add_argument("-i", "--indFn", type=str, metavar="<IND FILE>", required=True,
                        help="Path to the input .ind file.")
    parser.add_argument("-o", "--Output", type=str, metavar="<OUTPUT PREFIX>", required=False,
                        help="Output file prefix. Creates <PREFIX>.geno, <PREFIX>.snp, <PREFIX>.ind.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-C", "--Check", type=argparse.FileType("r"), metavar="<IND FILE>",
                       help="Check for duplicate individuals between -i and this .ind file.")
    group.add_argument("-E", "--Extract", action="store_true",
                       help="Extract selected individuals into a new EigenStrat database.")
    group.add_argument("-R", "--Remove", action="store_true",
                       help="Remove selected individuals, creating a new EigenStrat database.")

    group2 = parser.add_mutually_exclusive_group(required=False)
    group2.add_argument("-L", "--SampleList", type=argparse.FileType("r"), metavar="<LIST FILE>",
                        help="File with one sample name per line. Mutually exclusive with -S.")
    group2.add_argument("-S", "--Sample", action="append", metavar="<INDIVIDUAL>",
                        help="Individual to extract/remove. Can be called multiple times.")

    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(VERSION))

    args = parser.parse_args()

    # Validation
    if (args.Extract or args.Remove) and args.SampleList is None and args.Sample is None:
        parser.error("Sample (-S) or SampleList (-L) required for -R/-E functions.")
    if (args.Extract or args.Remove) and args.Output is None:
        parser.error("Output files (-o) must be specified when using -R/-E functions.")

    # --- Check mode ---
    if args.Check is not None:
        check_duplicates(args.indFn, args.Check)
        sys.exit(0)

    # --- Extract / Remove modes ---
    validate_eigenstrat(args.genoFn, args.snpFn, args.indFn)

    # Collect sample names
    samples = []
    if args.SampleList is not None:
        for line in args.SampleList:
            fields = line.strip().split()
            if fields and not fields[0].startswith("#"):
                samples.append(fields[0])

    if args.Sample is not None:
        for s in args.Sample:
            if s not in samples:
                samples.append(s)

    mode_label = "REMOVAL" if args.Remove else "EXTRACTION"
    print("Detected", len(samples), "individuals for {}.".format(mode_label), sep=" ", file=sys.stderr)

    # Parse the ind file once and build the index
    ind_entries = parse_ind_file(args.indFn)
    index, sex, pop = build_sample_index(samples, ind_entries)
    print("Indexed", len(index), "individuals.", sep=" ", file=sys.stderr)

    if args.Extract:
        extract_individuals(args.genoFn, args.snpFn, args.Output, samples, index, sex, pop)
    elif args.Remove:
        remove_individuals(args.genoFn, args.snpFn, args.indFn, args.Output, samples, index)


if __name__ == "__main__":
    main()
