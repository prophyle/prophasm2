# ProphAsm2

[![ProphAsm test](https://github.com/prophyle/prophasm2/actions/workflows/ci.yml/badge.svg)](https://github.com/prophyle/prophasm2/actions/)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisities](#prerequisities)
* [Getting started](#getting-started)
* [How to use](#how-to-use)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Contact](#contact)

<!-- vim-markdown-toc -->

## Introduction

ProphAsm2 is a versatile tool for computing simplitigs/SPSS
from *k-mer sets* and for *k-mer set operations*.
The new features compared to the original [ProphAsm](https://github.com/prophyle/prophasm)
include a largely speed and memory optimization, parallelization,
support for k-mer sizes up to 128 and support for minimum abundances.

Various types of sequencing datasets can be used as the input for
ProphAsm, including genomes, pan-genomes, metagenomes or sequencing reads.
Besides computing simplitigs, ProphAsm can also compute intersection
and set differences of k-mer
sets (while set unions are easy to compute simply by merging the source files).

Upon execution, ProphAsm first loads all specified datasets (see the `-i`
param) and the corresponding k-mer sets (see the `-k` param). If the `-x` param
is provided, ProphAsm then computes their intersection, subtracts the
intersection from the individual k-mer sets and computes simplitigs for the
intersection. If output files are specified (see the `-o` param), it computes
also set differences.



## Prerequisites

* GCC 4.8+ or equivalent
* ZLib


## Getting started

Download and compile ProphAsm:

```
git clone https://github.com/prophyle/prophasm2
cd prophasm2 && make -j
```

Compute simplitigs:

```
./prophasm -k 31 -i tests/test1.fa -o simplitigs.fa
```


## How to use

Set operations:
```
./prophasm -k 31 -i tests/test1.fa -i tests/test2.fa -o _out1.fa -o _out2.fa -x _intersect.fa -s _stats.tsv
   ```


## Command-line arguments

<!---
USAGE-BEGIN
-->
```


```
<!---
USAGE-END
-->


## Algorithm

In its core, ProphAsm2 uses the original algorithm for rapid computation of simplitigs as described in the [simplitig paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02297-z).

```python
def extend_simplitig_forward (K, simplitig):
	extending = True
	while extending:
		extending = False
		q = simplitig[-k+1:]
		for x in ['A', 'C', 'G', 'T']:
			kmer = q + x
			if kmer in K:
				extending = True
				simplitig = simplitig + x
				K.remove (kmer)
				K.remove (reverse_complement (kmer))
				break
	return K, kmer

def get_maximal_simplitig (K, initial_kmer):
	simplitig = initial_kmer
	K.remove (initial_kmer)
	K.remove (reverse_complement (initial_kmer))
	K, simplitig = extend_simplitig_forward (K, simplitig)
	simplitig = reverse_complement (simplitig)
	K, simplitig = extend_simplitig_forward (K, simplitig)
	return K, simplitig

def compute_simplitigs (kmers):
	K = set()
	for kmer in kmers:
		K.add (kmer)
		K.add (reverse_complement(kmer))
	simplitigs = set()
	while |K|>0:
		initial_kmer = K.random()
		K, simplitig = get_maximal_simplitig (K, initial_kmer)
		simplitigs.add (simplitig)
	return simplitigs
```



## Issues

Please use [Github issues](https://github.com/prophyle/prophasm2/issues).


## Changelog

See [Releases](https://github.com/prophyle/prophasm22/releases).


## Licence

[MIT](https://github.com/prophyle/prophasm/blob/master/LICENSE)


## Contact

[Ondrej Sladky](https://iuuk.mff.cuni.cz/~sladky/) \<ondra.sladky@gmail.com\>\
[Karel Brinda](https://brinda.eu) \<karel.brinda@inria.fr\>
