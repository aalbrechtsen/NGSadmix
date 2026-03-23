# NGSadmix Tutorial

## Repository Layout

- `./NGSadmix` is the compiled executable
- `./NGSadmix.cpp` is the source code
- `./Demo/Data/` contains the tutorial input files
- `./Demo/Results/` is the recommended output directory for tutorial runs

The tutorial data included in the repository:

- `Demo/Data/Demo1input.gz`
- `Demo/Data/Demo1pop.info`
- `Demo/Data/Demo2input.gz`
- `Demo/Data/Demo2pop.info`

Generated output files such as `.log`, `.filter`, `.qopt`, `.fopt.gz`, and plot images should normally stay out of git.

## Build

Compile the program with:

```bash
g++ NGSadmix.cpp -O3 -lpthread -lz -o NGSadmix
```

## Setup

Run all commands below from the repository root.

Create an output directory:

```bash
mkdir -p Demo/Results
```

Check that the expected input files are present:

```bash
ls Demo/Data
```

<details>
<summary>Example output</summary>

```text
Demo1input.gz
Demo1pop.info
Demo2input.gz
Demo2pop.info
```

</details>

## Input Format

NGSadmix expects genotype likelihoods in Beagle format.

- Column 1: marker name
- Column 2: allele 1
- Column 3: allele 2
- Then 3 columns per individual:
  first homozygote likelihood for allele 1, then heterozygote, then homozygote for allele 2

The likelihoods for each individual at a site should sum to a positive value. In many Beagle files they are normalized to sum to 1, but they are still likelihoods, not posterior genotype probabilities.

Take a quick look at the first dataset:

```bash
gunzip -c Demo/Data/Demo1input.gz | head -n 10 | cut -f 1-10 | column -t
gunzip -c Demo/Data/Demo1input.gz | wc -l
```

<details>
<summary>Example output</summary>

```text
marker      allele1  allele2  Ind0      Ind0      Ind0      Ind1      Ind1      Ind1      Ind2
1_15000765  2        0        0.000070  0.333333  0.666596  0.000000  1.000000  0.000000  0.969666
1_15001337  3        1        0.941155  0.058845  0.000000  0.941099  0.058901  0.000000  0.888859
1_15001470  0        2        0.969684  0.030316  0.000000  0.888825  0.111175  0.000000  0.999024
1_15001480  3        1        0.969682  0.030318  0.000000  0.001054  0.995607  0.003339  0.998049
1_15001731  3        2        0.000000  0.058903  0.941097  0.000266  0.999204  0.000530  0.000000
1_15001817  0        2        0.969632  0.030368  0.000000  0.665549  0.333333  0.001117  0.999023
1_15002813  1        3        0.799979  0.200021  0.000000  0.001058  0.998937  0.000004  0.984610
1_15002941  2        0        0.888810  0.111190  0.000000  0.888706  0.111294  0.000000  0.941138
1_15002958  2        0        0.799990  0.200010  0.000000  0.000000  0.999999  0.000001  0.941165
5617
```

</details>

Summarize the population labels with:

```bash
cut -f 1 -d " " Demo/Data/Demo1pop.info | sort | uniq -c
```

<details>
<summary>Example output</summary>

```text
     10 CEU
     10 JPT
     10 YRI
```

</details>

## Example 1: Small Three-Population Dataset

This dataset contains 30 individuals:

- 10 CEU
- 10 JPT
- 10 YRI

Run NGSadmix with `K=3`:

```bash
./NGSadmix \
  -likes Demo/Data/Demo1input.gz \
  -K 3 \
  -minMaf 0.05 \
  -seed 1 \
  -o Demo/Results/Demo1NGSadmix
```

<details>
<summary>Example output</summary>

```text
Input: lname=Demo/Data/Demo1input.gz nPop=3, fname=(null) qname=(null) outfiles=Demo/Results/Demo1NGSadmix
Setup: seed=1 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=5616 nind=30
Input file has dim (AFTER filtering): nsites=5616 nind=30
	[ALL done] cpu-time used =  1.41 sec
	[ALL done] walltime used =  2.00 sec
best like=-113212.723533 after 324 iterations
```

</details>

This produces:

- `Demo1NGSadmix.log`: run settings and convergence summary
- `Demo1NGSadmix.filter`: site-level filtering summary
- `Demo1NGSadmix.qopt`: inferred admixture proportions per individual
- `Demo1NGSadmix.fopt.gz`: inferred allele frequencies per site and ancestral population

Inspect the outputs:

```bash
cat Demo/Results/Demo1NGSadmix.log
zcat Demo/Results/Demo1NGSadmix.fopt.gz | head -n 5
head -n 5 Demo/Results/Demo1NGSadmix.qopt
```

<details>
<summary>Example output</summary>

```text
Input: lname=Demo/Data/Demo1input.gz nPop=3, fname=(null) qname=(null) outfiles=Demo/Results/Demo1NGSadmix
Setup: seed=1 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=5616 nind=30
Input file has dim (AFTER filtering): nsites=5616 nind=30
	[ALL done] cpu-time used =  1.41 sec
	[ALL done] walltime used =  2.00 sec
best like=-113212.723533 after 324 iterations
0.32468940525360973082 0.00000001553018945903 0.19270115449468547264
0.43873052722421068683 0.58762409251919067721 0.33873459854532766977
0.33956235147554725273 0.46603492795355749845 0.29773789076870887937
0.40108656038516710129 0.54319294164788478607 0.35895680105348198863
0.47057389707330554707 0.05879453926525243790 0.17359600404874098167
0.00000000099999999998 0.00000000099999999998 0.99999999800000005656
0.00000000099999999994 0.00000000676328993035 0.99999999223671010018
0.77973425437547561057 0.22026574462452436221 0.00000000099999999998
0.99999999799999994554 0.00000000099999999990 0.00000000099999999990
0.00000000099999999998 0.61097831929236390280 0.38902167970763606997
```

</details>

### Plot Example 1 in R

Interactive R:

```r
pop <- read.table("Demo/Data/Demo1pop.info", as.is = TRUE)[, 1]
q <- read.table("Demo/Results/Demo1NGSadmix.qopt")
ord <- order(pop)
par(mar = c(7, 4, 1, 1))
barplot(
  t(q)[, ord],
  col = c(2, 1, 3),
  names.arg = pop[ord],
  las = 2,
  ylab = "Demo1 admixture proportions",
  cex.names = 0.75
)
```

Non-interactive PNG output:

```bash
Rscript -e 'png("Demo/Results/Demo1NGSadmix.png", width=1200, height=700); pop<-read.table("Demo/Data/Demo1pop.info", as.is=TRUE)[,1]; q<-read.table("Demo/Results/Demo1NGSadmix.qopt"); ord<-order(pop); par(mar=c(7,4,1,1)); barplot(t(q)[,ord], col=c(2,1,3), names.arg=pop[ord], las=2, ylab="Demo1 admixture proportions", cex.names=0.75); dev.off()'
```

<details>
<summary>Example output</summary>

```text
null device
          1
```

</details>

![Demo 1 admixture plot](Demo/Results/Demo1NGSadmix.png)

## Example 2: Larger Dataset

The second example uses 50,000 sites from 100 individuals from five populations:

- ASW
- CEU
- CHB
- MXL
- YRI

The Example 2 input files are:

- `Demo/Data/Demo2input.gz`
- `Demo/Data/Demo2pop.info`

Summarize the population labels with:

```bash
cut -f 1 -d " " Demo/Data/Demo2pop.info | sort | uniq -c
```

<details>
<summary>Example output</summary>

```text
     20 ASW
     20 CEU
     20 CHB
     20 MXL
     20 YRI
```

</details>

Run NGSadmix with `K=3`:

```bash
./NGSadmix \
  -likes Demo/Data/Demo2input.gz \
  -K 3 \
  -P 1 \
  -minMaf 0.05 \
  -seed 4 \
  -o Demo/Results/Demo2NGSadmixK3
```

<details>
<summary>Example output</summary>

```text
Input: lname=Demo/Data/Demo2input.gz nPop=3, fname=(null) qname=(null) outfiles=Demo/Results/Demo2NGSadmixK3
Setup: seed=4 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=50000 nind=100
Input file has dim (AFTER filtering): nsites=49475 nind=100
	[ALL done] cpu-time used =  34.20 sec
	[ALL done] walltime used =  34.00 sec
best like=-3865964.415563 after 225 iterations
```

</details>

Then compare with `K=4`:

```bash
./NGSadmix \
  -likes Demo/Data/Demo2input.gz \
  -K 4 \
  -P 1 \
  -minMaf 0.05 \
  -seed 4 \
  -o Demo/Results/Demo2NGSadmixK4
```

<details>
<summary>Example output</summary>

```text
Input: lname=Demo/Data/Demo2input.gz nPop=4, fname=(null) qname=(null) outfiles=Demo/Results/Demo2NGSadmixK4
Setup: seed=4 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=50000 nind=100
Input file has dim (AFTER filtering): nsites=49475 nind=100
	[ALL done] cpu-time used =  52.22 sec
	[ALL done] walltime used =  52.00 sec
best like=-3817738.015534 after 343 iterations
```

</details>

Inspect the run summaries:

```bash
tail -n 20 Demo/Results/Demo2NGSadmixK3.log
tail -n 20 Demo/Results/Demo2NGSadmixK4.log
```

<details>
<summary>Example output</summary>

```text
Input: lname=Demo/Data/Demo2input.gz nPop=3, fname=(null) qname=(null) outfiles=Demo/Results/Demo2NGSadmixK3
Setup: seed=4 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=50000 nind=100
Input file has dim (AFTER filtering): nsites=49475 nind=100
	[ALL done] cpu-time used =  34.20 sec
	[ALL done] walltime used =  34.00 sec
best like=-3865964.415563 after 225 iterations
Input: lname=Demo/Data/Demo2input.gz nPop=4, fname=(null) qname=(null) outfiles=Demo/Results/Demo2NGSadmixK4
Setup: seed=4 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=50000 nind=100
Input file has dim (AFTER filtering): nsites=49475 nind=100
	[ALL done] cpu-time used =  52.22 sec
	[ALL done] walltime used =  52.00 sec
best like=-3817738.015534 after 343 iterations
```

</details>

### Plot Example 2 in R

For `K=3`:

```r
pop <- read.table("Demo/Data/Demo2pop.info", as.is = TRUE)
q <- read.table("Demo/Results/Demo2NGSadmixK3.qopt")
ord <- order(pop[, 1])
barplot(
  t(q)[, ord],
  col = 2:10,
  space = 0,
  border = NA,
  xlab = "Individuals",
  ylab = "Demo2 admixture proportions for K=3"
)
text(tapply(1:nrow(pop), pop[ord, 1], mean), -0.05, unique(pop[ord, 1]), xpd = TRUE)
abline(v = cumsum(sapply(unique(pop[ord, 1]), function(x) sum(pop[ord, 1] == x))), col = 1, lwd = 1.2)
```

For `K=4`, replace `Demo2NGSadmixK3.qopt` with `Demo2NGSadmixK4.qopt`.

Non-interactive PNG output:

```bash
Rscript -e 'png("Demo/Results/Demo2NGSadmixK3.png", width=1200, height=700); pop<-read.table("Demo/Data/Demo2pop.info", as.is=TRUE); q<-read.table("Demo/Results/Demo2NGSadmixK3.qopt"); ord<-order(pop[,1]); par(mar=c(5,4,2,1)); barplot(t(q)[,ord], col=2:10, space=0, border=NA, xlab="Individuals", ylab="Demo2 admixture proportions for K=3"); text(tapply(1:nrow(pop), pop[ord,1], mean), -0.05, unique(pop[ord,1]), xpd=TRUE); abline(v=cumsum(sapply(unique(pop[ord,1]), function(x) sum(pop[ord,1]==x))), col=1, lwd=1.2); dev.off()'

Rscript -e 'png("Demo/Results/Demo2NGSadmixK4.png", width=1200, height=700); pop<-read.table("Demo/Data/Demo2pop.info", as.is=TRUE); q<-read.table("Demo/Results/Demo2NGSadmixK4.qopt"); ord<-order(pop[,1]); par(mar=c(5,4,2,1)); barplot(t(q)[,ord], col=2:10, space=0, border=NA, xlab="Individuals", ylab="Demo2 admixture proportions for K=4"); text(tapply(1:nrow(pop), pop[ord,1], mean), -0.05, unique(pop[ord,1]), xpd=TRUE); abline(v=cumsum(sapply(unique(pop[ord,1]), function(x) sum(pop[ord,1]==x))), col=1, lwd=1.2); dev.off()'
```

<details>
<summary>Example output</summary>

```text
null device
          1
null device
          1
```

</details>

![Demo 2 admixture plot for K=3](Demo/Results/Demo2NGSadmixK3.png)

![Demo 2 admixture plot for K=4](Demo/Results/Demo2NGSadmixK4.png)

## Practical Notes

- The output directory should usually not be committed to git.
- `-seed` controls the random initialization and makes runs reproducible.
- `-K` is the assumed number of ancestral populations.
- `-P` sets the number of CPU threads.
- `-minMaf 0.05` filters out low-frequency sites before inference.

To compare models more carefully, run several seeds for the same `K` and compare the best likelihoods in the `.log` files.
