# Unique molecular identifiers (UMI) - de-duplication
Processes a raw FASTQ file by truncating UMIs and collapsing unique ones.


## Overview
Unique molecular identifiers (UMIs)are short, random sequences attached to RNA molecules at the beginning of library preparation. They allow to identify duplicates of distinct RNA molecules and therefore remove PCR amplification bias. Reads with the same UMI and aligned to the same genomic regions can be collapsed to asingle representative read. However, errors in base-calling of the UMI-sequence causes increased number of unique sequences found in the sample and reduces precision of original transcripts abundance quantification.

## Installation
To run UMI de-duplication, just clone the git:

```
git clone https://github.com/leachim/project-krakow.git
```

change to the directory:

```
cd project-krakow
```

and install the python package:

```
pip install -e .
```

Alternatively, you can also directly install the package from github:

```
pip install git+ssh://git@github.com/leachim/project-krakow.git#egg=umidedup
```

## Usage
### Input

* Called raw reads in FASTQ format. The file extension needs to be '.fq' or '.fastq'.
UMIs can be either truncated and added to the read description as 'UMI_<SEQUENCE>' or be still part of the read. In the latter case, the UMI base pair size needs to be given as an argument.

### Output

* A modified FASTQ file, containing only the unique read without UMIs.
* A statistics file, summarizing how many UMIs were collapsed.

### Arguments

* input [str] - Absolute or relative path to the FASTQ file
* length \[int\] (optional) - UMI size in base pairs
* algorithm [str] - De-duplication algorithm
* treshold \[int\] - Cutoff threshold for the chosen algorithm
* output [str] - Absolute or relative path to an output directory. If the directory does not exist, it will be created.

### Example script

Run the script `example/example-comparison.py` to see how to use the methods in practice.

### Run the example wrapper

Change to the right directoy (if not done yet):

```
cd project-krakow
```

and run the UMI de-duplication example:

```
python3 example/example-wrapper.py -i example/example_data/reads.fq -l 7 -a alevin -o example/example_output
```

## De-duplication - Algorithms
**Hamming-distance**

Compute pairwise hamming distance between UMI sequences and merge those UMIs that have a hamming distance smaller or equal to the threshold value, the merging is done in a greedy manner, starting and growing those UMIs that have the most reads, no reads are discarded.

> ### Note
> Treshold argument required

**Frequency**

Remove UMIs with less than threshold reads from the fastq file.

> ### Note
> Treshold argument required

**Alevin**

 This de-duplication algorithm is based on the work from [Srivastava et al.](https://doi.org/10.1186/s13059-019-1670-y). A detailed description of their full implementation can be found [here](https://salmon.readthedocs.io/en/latest/alevin.html). In short, it constructs a UMI graph, where UMIs are nodes and two UMIs are connected via an edge if their Hamming distance is 1. Subsequently, motivated by the principle of parsimony, edges are collapsed until a minimum number of UMIs can explain the data.


## Simulation of test data
To simulate test data, you can run the following command:

```
python3 tests/simulate_data.py -n <NO_READS> -d <AVG_READ_DEPTH> -seq_e <SEQ_ERROR> -o <OUT_DIR>
```

It will simulate reads with PCR amplification and sequencing error, generate the corresponding FASTQ file, and save the true reads with corresponding simulated depth.


> ### Note
> To get more help on the simulate_data.py arguments, run:
> ```
> python3 tests/simulate_data.py --help
> ```
