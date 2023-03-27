# Motif Finder

Welcome to Motif Finder!
This is a command line utility that allows you to take a FASTA file, specify a few parameters, and (hopefully) get some motifs prevalent in the sequences.

## Installation

If you have the Rust toolchain installed, you can install `motif_finder` with:

`cargo install motif_finder`

If you do not have the Rust toolchain installed, you can install it (here)[rustup.rs]
If you don't want to install it, you can also use the precompiled binaries in the releases tab on the right for your platform

If your platform isn't included, you can build it for your platform by cloning this repository:

`git clone https://github.com/nithishbn/MotifFinder.git`

and running

`cargo build --release`

in the source directory.
This will leave an executable in the `target/release/` directory which you can then run in the command line:
`motif_finder`

## Data format

This tool technically accepts all FASTA files, but the way it's meant to be used is to use an interesting approach in motif finding.

### RNASeq

By using RNASeq data and aligning it back to a reference genome, we can identify the alignment sites of transcripts. Using these alignment sites, we can generate the set of sequences _x_ bp upstream of the site in which to look for motifs, specifically for transcription factor binding sites.

This method involves finding an organism with RNASeq data, a reference genome, and a few bioinformatics tools including [samtools](https://www.htslib.org/), [bamtools](https://github.com/pezmaster31/bamtools/wiki), and [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html).

## Examples

You can try to find the motifs present in `promoters.fasta`, a set of 4 promoters known in _P. tricornutum_, a relatively unknown diatom species.

### Gibbs Sampler

Gibbs Sampler is an algorithm that iteratively searches for the best set of motifs in a set of sequences and throws out motifs at random until all iterations are finished.

`motif_finder -i promoters.fasta -e 4 -k 10 -o promotifs.txt gibbs -t 100 -r 100`

### Randomized Motif Search

Randomized Motif Search is an algorithm that iteratively searches for the best set of motifs in a set of sequences and throws out motifs at random until the score cannot be improved anymore.

`motif_finder -i promoters.fasta -e 4 -k 10 -o promotifs.txt randomized -r 100`

### Median String

Median String is an algorithm that checks the hamming distance from each kmer from each sequence and returns the minimized kmer from all strings. This algorithm is incredibly slow but can result in very accurate but short kmers.
Be warned when using large k values.

`motif_finder -i promoters.fasta -e 4 -k 8 -o promotifs.txt median`

### Alignment

If you wish to align the motifs you've generated back to the sequences from which they were generated to identify the highest locally scored motif over all sequences, you can run the same commands as above but with the `-a` flag

`motif_finder -i promoters.fasta -e 4 -k 8 -o promotifs.txt randomized -r 100`

This will generate alignments for the motifs after identifying the motifs.
