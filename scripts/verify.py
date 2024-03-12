#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse


def verify_instance(fasta_path: str, k: int, complements: bool) -> bool:
    """
    Check if running ProphAsm2 on given fasta file produces the same set of k-mers as the original one.
    """
    args = ["./prophasm2", "-i", fasta_path, "-k", f"{k}", "-o", "./bin/simplitigs.fa", "-S"]
    if not complements:
        args.append("-u")
    subprocess.run(args)
    # in result; in original sequence; in result without complements; in original without complements; in merged file
    stats = [{}, {}, {}]
    runs = [
        (0, "./bin/simplitigs.fa", "simplitigs", complements),
        (1, fasta_path, "original", complements),
    ]
    for i, path, result, pass_complements in runs:
        args = ["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-o", f"./bin/{result}.jf", path]
        if pass_complements:
            args.insert(2, "-C")
        subprocess.run(args)
        with open(f"./bin/{result}_stats.txt", "w") as f:
            subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
        with open(f"./bin/{result}_stats.txt", "r") as f:
            for _ in range(4):
                key, value = f.readline().split()
                stats[i][key] = value
    # Count k-mers on merged file.
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/simplitigs.jf", "./bin/original.jf"])
    with open(f"./bin/merged_stats.txt", "w") as f:
        subprocess.run(["jellyfish", "stats", f"./bin/merged.jf"], stdout=f)
    with open(f"./bin/merged_stats.txt", "r") as f:
        for _ in range(4):
            key, value = f.readline().split()
            stats[2][key] = value
    distinct_key = "Distinct:"
    total_key = "Total:"
    if stats[0][distinct_key] != stats[1][distinct_key] or stats[0][distinct_key] != stats[2][distinct_key]:
        print("F")
        print(f"Failed: k={k}: expected orginal_distinct_count={stats[1][distinct_key]}, result_distinct_count={stats[0][distinct_key]} and merged_distinct_count={stats[2][distinct_key]} to be equal.")
        return False
    elif complements and stats[0][distinct_key] != stats[0][total_key]:
        print("W")
        print(f"Warning: k={k}: number of masked k-mers={stats[0][total_key]} is not minimal possible (minimum is {stats[0][distinct_key]}).")
    else:
        print(".", end="")
        sys.stdout.flush()
    return True


def main():
    # Initialize.
    if not os.path.exists("bin"):
        os.makedirs("bin")

    parser = argparse.ArgumentParser("check if ProphAsm2 outputs simplitigs which contain the same set of k-mers"
                                     "as the original sequence")
    parser.add_argument("--quick", help="if set do not check for full range of k", action="store_true")
    parser.add_argument("path", help="path to the fasta file on which ProphAsm2 is verified")
    args = parser.parse_args()

    success = True
    print("Testing ProphAsm2 outputs valid simplitigs on file " + args.path)
    for complements in [True, False]:
        for k in range(2, 33, 3 if args.quick else 1):
            success &= verify_instance(args.path, k, complements)
        print("")

    # Print status.
    if not success:
        print("Tests failed")
        exit(1)
    print("OK")

main()
