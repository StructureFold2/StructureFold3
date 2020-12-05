# StructureFold3 <img src='assets/sf3_logo.png' align='right' width='400px' />

This is a work in progress, porting Structurefold2 into Python3.

The manual will be reworked/updated into both markdown and pdf formats. Modules 
will be added one at a time as they are ported 
over/updated, so it could be a while before SF3 is ready for use; SF2 is still 
the 'official' version of StructureFold for the time being. 

All file formats will be fully interchangable between SF2 and SF3; data
processed by either should be cross-compatible.

**Software Dependencies**
+ [Python 3](https://www.python.org/)
+ [BioPython](https://biopython.org/)
+ [Numpy](https://numpy.org/)
+ [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
+ [SAMtools](http://samtools.sourceforge.net/)
+ [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

**Recommended Software**
+ [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
+ [RNAStructure](https://rna.urmc.rochester.edu/RNAstructure.html)
+ [R](https://www.r-project.org/)

## Updates and Errata

## Planned Updates

## Changes from SF2 to SF3

+ react_correlation and rtsc_correlation are now one module, sf3_rx_correlation.
It is now required to provide the accompanying <.fasta> file, rather than an option.
