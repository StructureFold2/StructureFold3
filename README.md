# StructureFold3 <img src='assets/sf3_logo.png' align='right' width='400px' />

StructureFold2 is being ported into Python3. 

+ SF2 and SF3 will be cross-compatible.
+ Modules will be added as they are ported/updated.
+ The manual will be reworked/updated into both markdown and pdf formats.
+ SF2 is still the 'official' version of StructureFold for the time being.
+ Modules are provisional and have not gone through extended testing yet. Use at your own risk!

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

## Changes from SF2 to SF3

+ react_correlation and rtsc_correlation are now one module, sf3_rx_correlation.
It is now required to provide the accompanying <.fasta> file used to generate the 
<.react> or <.rtsc> files, rather than this being an option.

+ rtsc_abundances is now sf3_rtsc_abundances. Much like other modules, it now
compresses the analysis of multiple files into a single unfied csv. 

+ Help menus have been further subdivided for increased clarity, options have
adopted a more unified framework. Most modules will no longer utilize glob
to pick input files automatically.

+ fasta_composition has been reworked into fasta_statistics. It is now
agnostic to the type of sequence (nucleotide/amino acid), and takes a
threshold value (default 1%) of the entire file to determine 
the columns to report.

+ rtsc_specificity has been reworked into sf3_rtsc_specificity. Input
options have been harmonized with the newer modules, and the output 
file is more consistently organized. 

## Upcoming Changes

+ Manual will be worked on, easier to import modules will come first.

+ Comprehensive testing against the SF2 versions of the modules.

+ Merge sam_filter and sam_to_rtsc, making this step one module. Although
samtools is ubiquitous, the future module will probably be pure Python. Therefore
samtools will probably be removed from the dependencies. The future module
will filter mappings while generating <.rtsc>. 

+ Rework batch_fold_rna entirely. The input files and options need updating,
and a consistency check needs to be added between the fasta/react. 
Although the module supports multi-threading, it currently distributes
one thread per RNA regardless of the RNA's length. The future version
will call Fold-smp and use all specified cores to fold one RNA at a time.