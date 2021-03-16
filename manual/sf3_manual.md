# StructureFold3
Tools to facilitate rapid analysis of whole transcriptome chemical probing data

## Introduction
StructureFold3 is the Python3 successor to [Structurefold2](https://github.com/StructureFold2/StructureFold2)
(@tack2018structurefold2), aiming to continue building functionality and ease of use.

**Software Dependencies**

+ [Python 3](https://www.python.org/)
+ [BioPython](https://biopython.org/)
+ [Numpy](https://numpy.org/)
+ [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (@martin2011cutadapt)
+ [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (@langmead2012fast)
+ [STAR](https://github.com/alexdobin/STAR) (@dobin2013star)

**Recommended Software**

+ [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
+ [RNAStructure](https://rna.urmc.rochester.edu/RNAstructure.html) (@reuter2010rnastructure)
+ [R](https://www.r-project.org/) (@team2013r)

StructureFold3 is a set of Python3 scripts designed to process Structure-seq
(@ding2015genome, @ding2014vivo, @ritchey2017structure) or other high-throughput libraries generated to yield 
reverse transcription stops via chemically probing RNA structure into reactivity values at the single nucleotide level (structurome).
These reactivity values may be used to guide folding with RNAstructure (@reuter2010rnastructure) or similar software, 
greatly enhancing structure prediction with *in-vivo* probing data. This package allows a user with a minimal computational
or statistical background to execute a number of mostly automated processes to assemble *in-vivo* reverse transcriptase stops collected by
chemical probing and derive reactivity values indicative of the probability of base pair strandedness, thus generating
a reasonable facsimile of the actual *in-vivo* RNA structurome.

Ultimately, the final design and analysis of the data is up to the experimenter; data may be piped into RNAstructure
to obtain predicted folds of transcripts, or the experimenter may choose to directly analyze the derived chemical
reactivity or raw RT stop values in any way that suits their experimental design with typical statistical packages, such as R
(@team2013r). StructureFold3 is amenable to a variety of chemical probes (e.g. DMS, glyoxal, SHAPE), provided a negative 
reagent control library is included. Structurefold3 is designed to be modularl; steps and programs can be modified to accomodate
unforseen experimental, chemical, or biological experimental constraints.

## Walkthrough

### Before you begin...

An accurate StructureFold3 analysis hinges on a prudent selection of the transcriptome of interest. We highly recommend
that the user take the time to carefully consider the version, quality of annotation, and size of the transcriptome to be 
used; the <.fasta> you initially choose that contains this information may not be changed during the course of an analysis.
Some Eukaryotes contain an excessive number of transcript isoforms per gene, which may complicate downstream analyses, and
it may be wise to prune to no more than two isoforms per gene to map against. Selecting the most abundant isoform based
on other techniques or experiments is an alternative approach to reduce the number of included transcripts per gene and
enhance the precision of reactivity values and the simplicity of all downstream analyses. Ensembl (yates2020ensembl) is a great place to
look for transcriptome (cDNA) files if one has not already been selected to map against, although absolutely any <.fasta>
file containing transcripts will work provided the transcripts are of sufficient length.

**Ensembl Transcriptome Sources**

+ [Ensembl](https://www.ensembl.org/info/data/ftp/index.html)
+ [Ensembl Plants](http://plants.ensembl.org/info/website/ftp/index.html)
+ [Ensembl Protists](http://protists.ensembl.org/info/website/ftp/index.html)
+ [Ensembl Bacteria](http://bacteria.ensembl.org/info/website/ftp/index.html)

We recommend mapping to the whole transcripts before partitioning the resulting derived information
into transcript features (e.g. 5'UTR, CDS, 3'UTR) in later steps so that reads overlapping two of these
features will not encounter trouble mapping. To ensure that transcript annotation is sufficient to later
subdivide features for individual analysis, check to ensure the start and stop coordinates of each of
these features are included with the annotation or an associated <.gtf> or <.gff> file. At this time
we are unable to provide an automated script to handle annotation and automatically perform this
subdividison due to the lack of a consistent annotation file standard, but it could be included in a future version.

StructureFold3 is not a computationally demanding package by itself. Bowtie2, the recommended short read aligner,
is likewise lightweight and not exceptionally demanding on most hardware, even on larger transcriptomes and data sets.
During the course of data processing however, many intermediate files will be generated; thus, having at least ~4 TB of
space available is advised to store these files until an analysis is completed. Using RNAStructure 
to interpret the restraints/constraints generated by StructureFold3 with the included scripts to batch-fold thousands
of transcripts is computationally demanding; this can tie up multiple cores for several days. Thus, access to a server
or a dedicated personal machine is advised if carrying out this particular task. We recommend running StructureFold3
on either any mainstream Linux system or MACOS. It is possible to run StructureFold3 on a Windows install within an
emulated environment, but all features may not work properly and no additional support can be given.

### Installation and Support
Everything is hosted on the github [repository](https://github.com/StructureFold2/StructureFold3),
including this manual. A typical install would be accomplished via git.
```
git clone https://github.com/StructureFold2/StructureFold3
```
Next, this directory would be added to the path, typically by adding information to your shell profile
```
PATH=$PATH:/home/path/to/new/folder/StructureFold3
```
Lastly, whatever permissions are most applicable must now be set.
```
chmod -R 755 StructureFold3/
```
Make sure to check the dependencies (requirements) of StructureFold3 and install these according to any instructions
provided on respective websites and documentation. Many of these packages are commonly available from package managers
native to either Linux or MacOS systems. For a listing of these dependencies, check the top of page 1. Depending how
these are installed, permissions and the path may need to be set up for each dependency. StructureFold3 will automatically
call some of these programs. To check that  you have successfully installed StructureFold3 and all of the core dependencies,
try executing the test script, sf3_test_sf3.py, in a new shell. The script will notify the user of any missing dependencies
that are not available in the path or for which the user does not have adequate permissions to run. If you plan to customize
the StructureFold3 pipeline and forgo the use of suggested trimming software, these dependencies do not need to be installed and set up.

This manual makes the best attempt possible to detail the present StructureFold3 tool set and demonstrate its use. 
However, it is nearly impossible to plan for all contingencies. While the software is provided ‘as is’, we are eager 
to improve and further streamline our analysis pipeline. If you encounter a bug or require support, or have an idea for 
a feature that would improve your overall StructureFold3 experience, please do not hesitate to contact via GitHub. At some
point in the future, SF3 may be available via PIP, etc.


### Segment 1, Generating Reactivity Profiles
The intial parts of StructureFold3 are linear. Reads are trimmed, mapped, then converted into a custom format
(<.rtsc>, reverse transcriptase stop counts), before these are used to calculate reactivity (<.react>). Other 
trimming and/or mapping software can be subsituted readily, as the first two modules are of convenience; they 
run programs under the reccomended settings and all output is properly named and organized. Any other trimming/mapping
protocols can be substituted, provided the end result is in <.sam> format.

|Module                                         | Function                         |In   |Out  |
|:----------------------------------------------|:---------------------------------|-----|-----|
|[sf3_fastq_trim](#sf3_fastq_trim.py)           | Batch runs read trimming         |fastq|fastq|
|[sf3_fastq_map](#sf3_fastq_map.py)             | Batch runs read mapping          |fastq|sam|
|[sf3_sam_to_rtsc](#sf3_sam_to_rtsc.py)         | Extracts read stops              |sam  |rtsc|
|[sf3_rtsc_specificity](sf3_rtsc_specificity.py)| Calculates Reagent Specificity   |rtsc |csv|
|[sf3_rtsc_coverage](sf3_rtsc_coverage.py)      | Calculates transcript coverage   |rtsc |csv,txt|
|[sf3_rtsc_to_react](#sf3_rtsc_to_react.py)     | Calculates reactivity            |rtsc |react|


**Trimming**

Trimming is accomplished via the [sf3_fastq_trim](#sf3_fastq_trim.py) module, which drives
cutadapt (@martin2011cutadapt). A typical analysis would simply run the module within the directory
containing the raw <.fastq> files, generating a trimmed version of each and writing cutadapt's detailed
output as a text file.
```
sf3_fastq_trim.py -directory 
```
The actual run of reverse transcriptase is infered by extacting the insert between the
StructureSeq adapters. Inspect the log as a diagnostic. When the 3' sequencing adapter
is trimmed, everything downstream of it is removed. Excessive loss of read length after 3'
adapter trimming may suggest the RNA was low quality or partially degrated before the adapters were ligated.
While trimming in general for most RNA based techniques is becoming sort of obselete given the option
to softclip when mapping, it retains important for structure probing.

**Mapping**
Mapping is accomplished via the [sf3_fastq_map](#sf3_fastq_map.py) module, which drives
Bowtie2 (@langmead2012fast). A typical analysis would simply run the module within the directory 
containing the raw <.fastq> files, generating a mapped <.sam> of each and writing a detailed log 
output as a text file.
```
sf3_fastq_map.py -directory 
```

**Extract Read Stops**
Reverse transcriptase stops are extracted from mapped reads
```
sf3_sam_to_rtsc.py 
```

**Check Reagent Specificity**
The specificity of libraries can be ascertained such that it matches the predicted pattern
```
sf3_rtsc_specificity.py
```

**Calculate Read Coverage**
The per-transcript coverage for one or multiple files is calculated
```
sf3_rtsc_coverage.py
```

**Calculate Reactivity**


### Segment 2: Reactivity Analysis

**Segment 2: Reactivity Analysis**

|Module                                          | Function                              |
|:-----------------------------------------------|:--------------------------------------|
|[sf3_react_statistics](#sf3_react_statistics.py)| Summarizes reactivity metrics         |


### Segment 3: Folding and Subsequent Analysis

**Segment 3: Folding and Subsequent Analysis**

|Module                                         | Function                              |
|:----------------------------------------------|:--------------------------------------|
|[sf3_batch_fold](#sf3_batch_fold.py)           | Batch runs folding programs           |




## Modules

### sf3_fastq_trim.py
This module drives cutadapt (@martin2011cutadapt) under the reccomended settings. By
default, the Structure-Seq2 5' and 3' adapterstrimmed, but alternative adapters can be removed (-tp, -fp).
Simularly, the minimum quality and minimum length parameters are set by default, but can be changed (-minlen, -minqual).
The maximum length of a read to be accepted can also be set (-maxlen); filtering reads by the maximum length after adapter
trimming is a circuitous way to only accept reads where adapters were both found and removed. For example, by setting 
the maximum allowed read length under the actual lenght of the raw reads, only reads where adapters have been both detected
and removed are accepted, ensuring that data going forward.....

```
Batch runs folding programs

optional arguments:
  -h, --help        show this help message and exit

Input:
  -fastq f [f ...]  Input Specific <.fastq> files
  -directory        Select all <.fastq> in a target directory

Settings:
  -fp FP            [default = TGAACAGCGACTAGGCTCTTCA] 5' adapter
  -tp TP            [default = GATCGGAAGAGCACACGTCTG] 3' adapter
  -minlen MINLEN    [default = 20] Minimum accepted sequence length
  -minqual MINQUAL  [default = 30] Minumum accepted base quality
  -maxlen MAXLEN    [default = None] Maximum seq length
  -cores CORES      [default = 4] Number of cores to use
  -nextseq          Use NextSeq/NovaSeq quality scores

Output:
  -suffix SUFFIX    [default = trimmed] Trimmed <.fastq> file suffix
  -logname LOGNAME  [default = trim_log] Name of the trim log file

```

### sf3_fastq_map.py
This module does not yet have an SF3 implementation.

### sf3_sam_to_rtsc.py
Generates <.rtsc> from <.sam> files.

**Usage**
```
Converts <.sam> into reverse transcriptase stop files <.rtsc>

optional arguments:
  -h, --help            show this help message and exit

Input:
  fasta                 Index Fasta File
  sam                   Input SAM file(s)

Settings:
  -mismatches <number>  [default = 3] Maximum allowed mismatches/indels
  -firstmm              Accept alignments with first base mismatches
  -reverse              Accept alignments to the reverse strand
  -rm_secondary         Remove secondary alignments

Output:
  -logname LOGNAME      [default = filter_log.csv] Log file name
```

### sf3_rtsc_coverage.py
Calculates reagent coverage on one or more <.rtsc> files

**Usage**
```
Calculates RT stop coverage from <.rtsc> file(s)

optional arguments:
  -h, --help   show this help message and exit

Input:
  fasta        Reference Fasta
  f            Input <.rtsc> files

Settings:
  -bases ACGT  [default = AC] Coverage Specificity
  -ot OT       [default = 1.0] Overlap file threshold

Output:
  -name NAME   Output file name
  -ol          Create an overlap file
  -on          Overlap file name
```

### sf3_fasta_statistics.py
This module takes any fasta and creates a per sequence composition
report. It is intended to be used on the transcriptome file used for
mapping, thus producing a <.csv> detailing transcript composition. This
can be matched to other processed data from the same transcripts to investigate
any potential relationships between composition and reactivity.

The module is agnostic to the actual characters used in the sequences
so it could be used on a fasta containing amino acids to obtain similar
per sequence composition frequencies. By default, the column containing
the sequence names will be headed with 'transcript' but this can be
specified by the user with the -description flag. By default, it
will generate columns for any character in the sequences at or above
1% of the grand total; configurable via the -threshold flag.

**Usage**
```
Creates a composition report for a <.fasta> file

optional arguments:
  -h, --help           show this help message and exit

Input:
  fasta                Input Fasta File

Settings:
  -threshold <number>  [default = 1.0] Minimum percent to report
  -description <name>  [default = transcript] Sequence column heading

Output:
  -affix AFFIX         [default = composition] <.csv> file affix
```

### sf3_react_heat_correct.py
Probing reagents can be more reactive at higher temperatures, thus one may wish to 
normalize this effect out from any two pairwise <.react>s derived from experiments
at different temperatures. For any two such <.react> files generated against the same
transcriptome, the sum of all base reactivties will be summed for both files, and used
to scale the the lower/higher temperature file up/down respectively, such that the total
amount of reactivity on all bases in both files is the same. 

**Usage**
```
Corrects two <.react>s for differential temperature

optional arguments:
  -h, --help      show this help message and exit

Input:
  lower           lower temp <.react> file
  higher          higher temp <.react> file

Output:
  -suffix SUFFIX  [default = corrected] Suffix for out files
```

### sf3_react_to_csv.py
For close or custom analysis, researchers may find the <.react> format clunky
or just not particularly amenable to individual analysis. This module converts <.react>
files into <.csv> files. A single transcript may be extracted to a <.csv> file 
or the entire file may be batch converted into a directory of <.csv> files. The
corresponding <.fasta> used to generate the <.react> files must be provided so
nucelotides can be matched to the reactivity values in the output <.csv>.

**Usage**
```
Reformats reactivity values to <.csv> format

optional arguments:
  -h, --help        show this help message and exit

Input:
  react             File to pull values from
  fasta             File used to generate the <.react>

Settings:
  -mode {S,M}       Single or Multi-Transcript mode
  -transcript NAME  [S] Specific Transcript to reformat
  -restrict <.txt>  [M] Limit analysis to these specific transcripts
  -outdir OUTDIR    [M] Name of the out directory, overrides default
```

### sf3_rtsc_abundances.py
It is possible to infer releative transcript abundance from <.rtsc> files.
This module will take one or more <.rtsc> files and report a selected abundance
metric (RPKM/TPM) for each transcript based on the total number of mapped RT stops.
All <.rtsc> being calculated together should be generated from mapping against the
same reference files, and thus should contain the exact same set of transcript entries.
However, by default transcripts without an entry in a file will result in an 'NA', 
unless changed by the user to zeros (-zero).

**Usage**
```
Determines approximate transcript abundance using <.rtsc> files

optional arguments:
  -h, --help        show this help message and exit

Input:
  rtsc              Input <.rtsc> file(s)

Settings:
  -mode {RPKM,TPM}  Abundance Metric to Calculate
  -zero             Set missing values to zero

Output:
  -name NAME        Specify output file name
```

### sf3_rtsc_specificity.py
This module calculates the reverse transcriptase stop specificity using
<.rtsc> files and the fasta they were generated against. Hence, comparing 
treated and untreated samples will detail the 

**Usage**
```
Analyzes native/reagent nucleotide RT stop specificity

optional arguments:
  -h, --help      show this help message and exit

Input:
  index           Fasta used to generate the <.rtsc>
  rtsc            Input <.rtsc>

Settings:
  -report REPORT  Include these nucelotides in report
  -round DIGITS   [default = 5] Decimal places to report

Output:
  -name NAME      Specify output file name
```

### sf3_sam_to_rtsc.py
This module filters mapped reads in SAM format, extracting the implied
reverse transcriptase (RT) stops that pass the default and/or user
defined criteria. These stops are then written to a Reverse Transcriptase
Stop Count <.rtsc> file to represent these stops in subsequent analysis steps.

**Usage**
```
Converts <.sam> into reverse transcriptase stop files <.rtsc>

optional arguments:
  -h, --help            show this help message and exit

Input:
  fasta                 Index Fasta File
  sam                   Input SAM file(s)

Settings:
  -mismatches <number>  [default = 3] Maximum allowed mismatches/indels
  -firstmm              Accept alignments with first base mismatches
  -reverse              Accept alignments to the reverse strand
  -rm_secondary         Remove secondary alignments

Output:
  -logname LOGNAME      [default = filter_log.csv] Log file name
```

### sf3_fasta_stepwise.py
This module is designed to break up large sequences into smaller
'sliding windows', so that small, granular segments can be folded
and further studied. 

**Usage**
```
Splits sequences into sliding window segments

optional arguments:
  -h, --help  show this help message and exit

Input:
  fasta       Input Fasta File

Settings:
  -size SIZE  [default = 10] Size of the Windows
  -step STEP  [default = 5] Step of the Windows

Output:
  -name NAME  Specify output file name
  -sort_file  Sort the output by name

```

### sf3_batch_fold.py
This module batch runs RNAStructure's Fold or partition on transcript or oligo
sets, with the option to use <.react> files as restraints/constraints. Thus,
RNAStructure must be installed, configured, and accesssible in the path. Input
sequences are prescreened to prevent errors when folding; the min and max length
tolerated can be assigned (-minlen, -maxlen), and any discrepancies between a
sequence and any restraints/constraints will be recorded (-errorname to change log name)
as a desynch error.

**Usage**
```
Batch runs folding programs

optional arguments:
  -h, --help            show this help message and exit

Input:
  fasta                 Reference Fasta File

Settings:
  -mode {R,V}           RNAStructure or Vienna
  -react <.react>       React file to use as restraints/constraints
  -restrict <.txt>      Limit folding to these specific transcripts
  -temp TEMP            [default = 310.15] Temperature to use for folding
  -cores CORES          [default = 4] Number of cores to use
  -distance DISTANCE    [default = 99999] Maximum pairing distance
  -minlen MINLEN        [default = 10] Minimum length to fold
  -maxlen MAXLEN        [default = 5000] Maximum length to fold
  -truncate TRUNCATE    [default = 0] Ignore <.react> for last N nucleotides
  -threshold TH         Apply hard constraints using this threshold

RNAStructure Settings:
  -slope SLOPE          [default = 1.8] Parameter for RNAstructure
  -intercept INTERCEPT  [default = -0.6] Parameter for RNAstructure
  -partition            Use the partition function
  -multiple             Output all structures rather than just MFE

Output:
  -errorname ERRORNAME  Name the error log
  -paraname PARANAME    Name the parameter log
  -outdir OUTDIR        Name the out directory
  -ct CT                Name the ct folder
  -ps PS                Name the ps folder
  -bp BP                Name the bp folder (Partition Only)
  -pfs PFS              Name the pfs folder (Partition Only)

```

### sf3_rtsc_to_react.py
This module uses two <.rtsc> files (reagent -/+) and the sequences they were generated
with (<.fasta>) to calculate per base reactivity. By default, the sample will be normalized
to its own 2-8% scale; invoking the -scale flag followed by a <.scale> file
will apply the 2-8% scale generated by another sample, thus calibrating both to a common normalization
scale between samples. The specificity of the reactivity caluculation is determined by the -bases 
flag (default AC).

**Usage**
```
Generates a <.react> file from two <.rtsc> files

optional arguments:
  -h, --help            show this help message and exit

Input:
  control               Control <.rtsc> file
  treatment             Reagent <.rtsc> file
  fasta                 Transcript <.fasta> file

Settings:
  -restrict <.txt>      Limit analysis to these specific transcripts
  -scale <.scale>       Provide a normalizaiton file for calculation
  -threshold THRESHOLD  [default = 7.0] Maximum Reactivity Cap
  -ln_off               Do not take the natural log of the stop counts
  -nrm_off              Do not perform final 2-8% reactivity normalization
  -bases ACGT           [default = AC] Reaction Specificity

Output:
  -save_fails           Log transcripts with zero or missing scales
  -name NAME            Specify output file name
```

## References

