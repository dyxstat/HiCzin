# HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts using zero-inflated negative binomial regression


Please refer to the folder 'RECOMB paper related materials' for supplementary materials, related datasets and codes in our HiCzin RECOMB paper (Yuxuan Du, Sarah M. Laperriere, Jed Fuhrman, and Fengzhu Sun.Journal of Computational Biology.ahead of printhttp://doi.org/10.1089/cmb.2021.0439).



## HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts
HiCzin is a normalization method designed for metagenomic Hi-C data. HiCzin is based on the zero-inflated negative binomial regression frameworks. We model the population of the intra-species contacts using the negative binomial distribution. HiCzin combines the counting distribution of the intra-species contacts with a mass distribution of unobserved contacts. The residues of the counting part serve as normalized contacts. HiCzin can also be used to detect and remove spurious contacts.


## Usage 
HiCzin is available through the Rscript HiCzin.R. 
