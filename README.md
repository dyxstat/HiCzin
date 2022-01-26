# HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts using zero-inflated negative binomial regression


Please refer to the folder 'RECOMB paper related materials' for supplementary materials, datasets and codes in our HiCzin RECOMB paper (Yuxuan Du, Sarah M. Laperriere, Jed Fuhrman, and Fengzhu Sun.Journal of Computational Biology.ahead of printhttp://doi.org/10.1089/cmb.2021.0439).



## HiCzin: Normalizing metagenomic Hi-C data and detecting spurious contacts
HiCzin is a normalization method designed for metagenomic Hi-C data. HiCzin is based on the zero-inflated negative binomial regression frameworks. We model the population of the intra-species contacts using the negative binomial distribution. HiCzin combines the counting distribution of the intra-species contacts with a mass distribution of unobserved contacts. The residues of the counting part serve as normalized contacts. HiCzin can also be used to detect and remove spurious contacts.

## Install and Setup
### conda
We recommend using conda to run HiCzin.

After installing Anaconda (or miniconda), Users can clone the repository with git:
```
git clone https://github.com/dyxstat/HiCzin.git
```

Once complete, you can enter the repository folder and then create a HiCzin environment using conda:
```
# enter the HiCzin folder
cd HiCzin
# Construct environment
conda env create -f HiCzin_linux_env.yaml 
or
conda env create -f HiCzin_osx_env.yaml
# Enter the environment
conda activate HiCzin_env
```

Normalization method in HiCzin depends on R package 'glmmTMB'. Though the R package can be installed by 'conda install -c conda-forge r-glmmtmb', you will meet one potential warning derived from the dependency version (https://github.com/glmmTMB/glmmTMB/issues/615): 'Error in .Call("FreeADFunObject", ptr, PACKAGE = DLL) : "FreeADFunObject" not available for .Call() for package "glmmTMB"' and we are not sure whether this warning would influence the noramlization results. To get rid of this warning, we strongly recommend you to install the source version of package 'glmmTMB' directly in R:

```
# Enter the R
R
# Download the R package and you may need to select a CRAN mirror for the installation
install.packages("glmmTMB", type="source")
```

Finally, you can test the pipeline, and testing result are in test/out/hiczin.log:
```
python ./hiczin.py test test/out
```


## HiCzin analysis
### Implement the HiCzin normalization
```
python ./hiczin.py norm [Parameters] FASTA_file BAM_file TAX_file COV_file OUTPUT_directory
```
#### Parameters
```
-e: Case-sensitive enzyme name. Use multiple times for multiple enzymes 
    (Optional; If no enzyme is input, HiCzin_LC mode will be automatically employed to do normalization)
--min-len: Minimum acceptable contig length (default 1000)
--min-mapq: Minimum acceptable alignment quality (default 30)
--min-match: Accepted alignments must be at least N matches (default 30)
--min-signal: Minimum acceptable signal (default 2)
--thres: Maximum acceptable fraction of incorrectly identified valid contacts in spurious contact detection (default 0.05)
-v: Verbose output
```
#### Input File

* ***FASTA_file***: a fasta file of the assembled contig (e.g. final.contigs.fa)
* ***BAM_file***: a bam file of the Hi-C alignment (e.g. MAP_SORTED.bam)
* ***TAX_file***: a csv file of contigs' taxonomy assignment by TAXAassign (e.g. contig_tax.csv)
* ***COV_file***: a txt file of contigs' coverage information computed by script: ‘jgi summarize bam contig depths’ from MetaBAT2 (e.g. depth.txt)


#### Example
```
python ./hiczin.py pipeline -e HindIII -e NcoI -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```
If the restriction enzymes employed in the experiment are unspecified, use
```
python ./hiczin.py pipeline -v final.contigs.fa MAP_SORTED.bam contig_tax.csv depth.txt out
```


### Output file
All outputs of HiCzin normalization pipeline are located in the OUT_directory you specified when running the software. 

* ***hiczin.log***: the specific implementation information of HiCzin pipeline
* ***contig_info.csv***: information of contigs, containing four columns (contigs' name, the number of restriction sites on contigs, contigs' length 
and coverage) or three columns (contigs' name, length and coverage)
* ***valid_contact.csv***: information of valid intra-species contacts, containing three columns (index of the first contig, index of the second contig, and value of raw valid  intra-species contacts)
* ***Normalized_contact_matrix.npz***: a sparse matrix of normlized Hi-C contact maps in csr format and can be reloaded using ***scipy.sparse.load_npz('Normalized_contact.npz')***
* ***HiCzin_normalized_contact.gz***: Compressed format of the normalized contacts and contig information by pickle. This file can further serve as the input of HiCBin and HiFine.


## Contacts and bug reports
If you have any questions or suggestions, welcome to contact Yuxuan Du (yuxuandu@usc.edu).


## Copyright and License Information
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.







