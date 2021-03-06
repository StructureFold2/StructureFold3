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
threshold value (default 1%) of the entire file composition to determine 
the columns to report.

+ rtsc_specificity has been reworked into sf3_rtsc_specificity. Input
options have been harmonized with the newer modules, and the output 
file is more consistently organized. 

+ sam_filter and sam_to_rtsc have been merged into sf3_sam_to_rtsc, 
encapsulating this step into a single module. Samtools has been removed
from the dependencies, as the new module performs all filtering entirely 
in Python without calling samtools or generating intermediate files.
Logfiles now have a column for number of reads excluded due to bitflag
based parameters, where in SF2, this would be the difference between
sam_lines and filtered_lines. See the new manual entry for more details.
The current StructureSeq2 protocol does not typically generate paired-end reads,
but building them into this rework would be much eaiser than using the
old methodology.

+ batch_fold_rna has been reworked into sf3_batch_fold.
This module has been entirely re-written from scratch. Input(s) are
pre-screened for all potential errors, details of which are logged to
a single file in the out directory. Run parameters are saved to a separate
parameter log in the out directory, rather than as the name of the out directory.
Instead of using 'Fold' and 'partition', the '-smp' versions of these
programs are called. Allocating one core per RNA turned out to be
inefficient, as one core was devoted to both long and short RNAs, hence
there was a lot of time waiting on the longer folds. This module is now
much more maintainable and upgradeable. Vienna package support is not
yet implemented pending any requests for it. See the manual entry for more details.

## Upcoming Changes

+ Manual will be worked on, easier to import modules will come first.

+ Comprehensive testing against the SF2 versions of the modules.
