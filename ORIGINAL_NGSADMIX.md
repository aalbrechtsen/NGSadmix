# The Original NGSadmix Program

This page covers the original `NGSadmix` program in this repository (`NGSadmix.cpp`). For new analyses, we recommend [admixer](https://github.com/aalbrechtsen/admixer), a faster implementation of the same model (see the [README](README.md)). The original program is still available and works as before.

## Installation

Clone and compile:

```bash
git clone https://github.com/aalbrechtsen/NGSadmix.git
cd NGSadmix
g++ NGSadmix.cpp -O3 -lpthread -lz -o NGSadmix
```

NGSadmix is also distributed as part of [ANGSD](https://github.com/ANGSD/angsd).

## Quick Start

Run NGSadmix on a Beagle genotype-likelihood file:

```bash
./NGSadmix -likes inputBeagleFile.gz -K 3 -o outFileName -P 10
```

Main arguments:

- `-likes`: Beagle file with genotype likelihoods (see [Input](README.md#input-beagle-genotype-likelihood-files))
- `-K`: number of ancestral populations
- `-o` or `-outfiles`: output prefix
- `-P`: number of threads

Example:

```bash
./NGSadmix -likes input.gz -K 3 -P 4 -o myoutfiles -minMaf 0.05
```

## Output Files

- `.log`: run settings, progress, likelihoods, and timing
- `.qopt`: estimated admixture proportions for each individual
- `.fopt.gz`: estimated allele frequencies for each site and ancestral population
- `.filter`: per-site filtering summary

## Options

Run `./NGSadmix` without arguments for the full list. Useful options:

| option | meaning | default |
|---|---|---|
| `-seed` | seed for the random initialization | |
| `-fname` | starting ancestral population frequencies | |
| `-qname` | starting admixture proportions | |
| `-printInfo 1` | print marker ID and mean MAF for the analysed SNPs | |
| `-method 0` | disable EM acceleration | |
| `-misTol` | threshold for treating a genotype likelihood as missing | `0.05` |
| `-tolLike50` | stop when the likelihood changes less than this over 50 iterations | `0.1` |
| `-tol` | convergence tolerance | `1e-5` |
| `-dymBound` | use dynamic boundaries | enabled |
| `-maxiter` | maximum EM iterations | `2000` |
| `-minMaf` | minimum minor allele frequency | `0.05` |
| `-minLrt` | minimum likelihood-ratio threshold for MAF > 0 | |
| `-minInd` | minimum number of informative individuals at a site | |

A single run can end in a local optimum, so run several seeds and compare the best likelihoods. [CONVERGENCE_TUTORIAL.md](CONVERGENCE_TUTORIAL.md) shows how.

## Tutorials

- [TUTORIAL.md](TUTORIAL.md): the demo data with the original NGSadmix program.
- [CONVERGENCE_TUTORIAL.md](CONVERGENCE_TUTORIAL.md): checking convergence across multiple runs.
- [ADMIXER_TUTORIAL.md](ADMIXER_TUTORIAL.md): the same data with admixer, and a comparison of the two programs.
