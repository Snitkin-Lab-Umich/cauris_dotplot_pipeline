# cauris_dotplot
This script takes input assemblies (in the form of .fasta files) and generates synteny dot plots from them. 

## Installation

In order to run this script on the Great Lakes cluster, load these modules first:

```
module load Bioinformatics mummer/4.0.0rc1
module load R/4.4.0
```

You will also need to install ggplot2 in R before running the script. 

## Generating Plots

Run the script with the following command:

```
python dotplot.py --query [query_path] --subject [subject_path]
```

This will generate your plots in the results/ directory. An explanation of each input is below.

### --query

(Required) This should be a path to a directory of .fasta files for the assemblies you want to align to your reference. Make sure these end with .fasta, .fa, or .fna. All queries will be aligned to the same reference.

### --subject

(Required) This should be a path to a .fasta file for the reference sequence you want to align everything to. Only one file is supported at the moment, so you will need to run the script multiple times if you have multiple reference sequences. A sample reference file for testing is included in the reference/ directory (from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002759435.3/).

### --name

(Optional) Provide a name for this batch. A directory with this name will be generated in results/.

### --highlight

(Optional) Provide a path to a table of regions to highlight. The selected regions of each contig will be indicated on the plot with a red box. This requires a specific format, and an example file is provided with highlight_data.tsv. You will also need to know the name of the contig or scaffold you want to highlight a region of, which can be determined via BLAST or with this tool: https://github.com/Snitkin-Lab-Umich/caurisblast. 

### --alignments

(Optional) Provide a path to an existing results/ directory. This will skip the Nucmer alignments, speeding up the script. This intended for cases where you want to remake plots from existing data. This won't work if the file names or structure of the results/ directory is changed. This option is mutually exclusive with providing your own query and subject.

## Output

If the script ran successfully, a few different output directories will be generated:

### plots/
This contains all of your dotplots as .pdf files.
### results/
This contains the Nucmer results (.coord and .delta files).
### temp/
This contains tables with information about the contigs in your assemblies, such as names and lengths.
### logs/
This contains a log file for the most recent run.

