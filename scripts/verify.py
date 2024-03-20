#!/usr/bin/env python3
import subprocess
import sys
import os
import argparse

def run_jellyfish_count(path: str, k: int, complements: bool, pass_m: int, result: str) -> dict:
    args = ["jellyfish", "count", "-m", f"{k}", "-s", "100M", "-L", f"{pass_m}", "-o", f"./bin/{result}.jf", path]
    if complements:
        args.insert(2, "-C")
    subprocess.run(args)

def run_prophasm2(fasta_path: str, k: int, complements: bool, m: int, result: str, interpath=None) -> str:
    args = ["./prophasm2", "-i", fasta_path, "-k", f"{k}",  "-S", "-m", f"{m}"]
    if interpath is None:
        args += ["-o", f"./bin/{result}.fa"]
    else:
        args += ["-i", interpath, "-x", f"./bin/{result}.fa"]
    if not complements:
        args.append("-u")
    subprocess.run(args)


def read_jellyfish_stats(result: str) -> dict:
    stats = {}
    with open(f"./bin/{result}_stats.txt", "w") as f:
        subprocess.run(["jellyfish", "stats", f"./bin/{result}.jf"], stdout=f)
    with open(f"./bin/{result}_stats.txt", "r") as f:
        for _ in range(4):
            key, value = f.readline().split()
            stats[key] = value
    return stats


def check_equal_kmers(stats: dict, k: int, m: int, complements: bool) -> bool:
    distinct_key = "Distinct:"
    total_key = "Total:"
    if stats[0][distinct_key] != stats[1][distinct_key] or stats[0][distinct_key] != stats[2][distinct_key]:
        print("F")
        print(f"Failed: k={k}, m={m}: expected orginal_distinct_count={stats[1][distinct_key]}, result_distinct_count={stats[0][distinct_key]} and merged_distinct_count={stats[2][distinct_key]} to be equal.")
        return False
    elif complements and stats[0][distinct_key] != stats[0][total_key]:
        print("W")
        print(f"Warning: k={k}, m={m}: number of masked k-mers={stats[0][total_key]} is not minimal possible (minimum is {stats[0][distinct_key]}).")
    else:
        print(".", end="")
        sys.stdout.flush()
    return True


def verify_instance(fasta_path: str, k: int, complements: bool, m: int) -> bool:
    """
    Check if running ProphAsm2 on given fasta file produces the same set of k-mers as the original one.
    """
    run_prophasm2(fasta_path, k, complements, m, "simplitigs")
    # in result; in original sequence; in result without complements; in original without complements; in merged file
    stats = [{}, {}, {}]
    runs = [
        (0, "./bin/simplitigs.fa", "simplitigs", 1),
        (1, fasta_path, "original", m),
    ]
    for i, path, result, pass_m in runs:
        run_jellyfish_count(path, k, complements, pass_m, result)
        stats[i] = read_jellyfish_stats(result)
    # Count k-mers on merged file.
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/simplitigs.jf", "./bin/original.jf"])
    stats[2] = read_jellyfish_stats("merged")
    return check_equal_kmers(stats, k, m, complements)

def verify_intersection(fasta_path: str, interpath: str, k: int, complements: bool, m: int) -> bool:
    run_prophasm2(fasta_path, k, complements, m, "simplitigs", interpath)
    # This assumes that no k-mer appears twice.
    run_prophasm2(fasta_path, k, complements, m, "inter1")
    run_prophasm2(interpath, k, complements, m, "inter2")
    run_jellyfish_count("./bin/inter1.fa", k, complements, 1, "inter1")
    run_jellyfish_count("./bin/inter2.fa", k, complements, 1, "inter2")
    subprocess.run(["jellyfish", "merge", "-L", "2", "-o", f"./bin/inter.jf", "./bin/inter1.jf", "./bin/inter2.jf"])
    stats = [{}, {}, {}]
    run_jellyfish_count("./bin/simplitigs.fa", k, complements, 1, "simplitigs")
    stats[0] = read_jellyfish_stats("simplitigs")
    stats[1] = read_jellyfish_stats("inter")
    subprocess.run(["jellyfish", "merge", "-o", f"./bin/merged.jf", "./bin/simplitigs.jf", "./bin/inter.jf"])
    stats[2] = read_jellyfish_stats("merged")
    return check_equal_kmers(stats, k, m, complements)



def main():
    # Initialize.
    if not os.path.exists("bin"):
        os.makedirs("bin")

    parser = argparse.ArgumentParser("check if ProphAsm2 outputs simplitigs which contain the same set of k-mers"
                                     "as the original sequence")
    parser.add_argument("--quick", help="if set do not check for full range of k", action="store_true")
    parser.add_argument("path", help="path to the fasta file on which ProphAsm2 is verified")
    parser.add_argument("--interpath", help="path to a second fasta on which intersection is verified")
    args = parser.parse_args()

    success = True
    if args.interpath:
        print("Testing ProphAsm2 outputs valid intersection on files " + args.path + " and " + args.interpath)
        for complements in [True]:
            for m in range(1, 3):
                for k in range(2, 33, 7 if args.quick else 1):
                    success &= verify_intersection(args.path, args.interpath, k, complements, m)
                print("")

    print("Testing ProphAsm2 outputs valid simplitigs on file " + args.path)
    for complements in [True, False]:
        for m in range(1, 4):
            for k in range(2, 33, 3 if args.quick else 1):
                success &= verify_instance(args.path, k, complements, m)
            print("")


    # Print status.
    if not success:
        print("Tests failed")
        exit(1)
    print("OK")

main()
