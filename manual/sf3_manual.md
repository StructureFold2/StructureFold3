# Work in progress!

* All modules to be detailed

* Walkthrough of a typical analysis

* All the information there is on StructureSeq/Fold

## Modules

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
  -nrm_off              Do not perform final 2-8％ reactivity normalization
  -bases ACGT           [default = AC] Reaction Specificity

Output:
  -save_fails           Log transcripts with zero or missing scales
  -name NAME            Specify output file name
```