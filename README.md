# HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts using zero-inflated negative binomial regression

## Source of biases in metagenomic Hi-C experiments
In addition to chromosomal contacts of interests, several other factors unrelated to chromosomal contacts can influence the number of Hi-C interactions between contigs. We refer to such factors as biases. Two kinds of biases with substantial influences on metagenomic Hi-C contact maps were reported: explicit biases and implicit biases. Explicit biases include three potential factors: the number of enzymatic restriction sites on contigs, the length and the coverage of contigs. Implicit biases contain unobserved interactions and spurious inter-species contacts. Unobserved interactions mean that some chimerical DNA fragments are missed due to the factors such as mappability of contigs and in vivo constraints on accessibility. Spurious inter-species contacts arise from the ligation of DNA fragments between closely related species. It is necessary to correct these biases before downstream analysis.

## HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts
HiCzin is a normalization method designed for metagenomic Hi-C data. HiCzin is based on the zero-inflated negative binomial regression frameworks. We model the population of the intra-species contacts using the negative binomial distribution. HiCzin combines the counting distribution of the intra-species contacts with a mass distribution of unobserved contacts. The residues of the counting part serve as normalized contacts. HiCzin can also be used to detect and remove spurious contacts.


## Usage 
HiCzin is available through the Rscript HiCzin.R. 
