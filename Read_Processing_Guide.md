# ChIP-Seq Read Processing Guide

This guide was made to instruct beginners to process their own data for visual inspection and downstream analysis.

Here, we will describe how to process your raw fastq files and get peaks and track files for IGV.

## Notes on bash scripts

First, a note on bash scripts using the Sun Grid Engine (SGE) environment. 
At the top of every script are the qsub options:

```
#!/bin/bash
#$ -cwd
#$ -t 1-12
#$ -V
#$ -hold_jid fastqc
#$ -N trimgalore

```

```#!/bin/bash``` allows this to be submitted as a batch job using ```qsub```.
```-cwd``` enforces that the script is run in the current working directory.
```-t``` is the task IDs to run, we use a prefix file to cycle through these IDs creating a job array in this part of each script :

```
myFile=`cat prefix.txt`
myFile=`head -n $SGE_TASK_ID prefix.txt | tail -1`
```

```-V``` allows bash to find programs in your path
```-hold_jid``` allows you to hold the job in the queue until the specified job has completed
```-N``` allows you to name the job

If your jobs are queued but not running for a long time try adding ```-l h_rt=06:00:00``` after qsub for a max run time of 6 hours,
the default option reserves a run time of 72 hours.

You may also add additional options when you run ```qsub```, reference http://gridscheduler.sourceforge.net/htmlman/manuals.html for more options.

1) mkdir for your project

2) cd to project and mkdir fastq, trim, map and peakcall

3) copy over all scripts and change the option for -t in each script to match the
number of samples you have (1-12 for 12, 1-8 for 8 etc.)

4) cd to fastq and make symlinks of your fastq files.

Make a file with all symlink commands called ln.sh, then run it:

```
ls ~/groupso/SO_DATA/2018-03/C*JC_R*.fastq.gz | column -t | awk '{print "ln -s "$1}' > ln.sh

bash ln.sh
```
If you mess up the symlinks use ```find -type l -delete``` to delete all symlinks

5) make prefix.txt as follows:

```
ls *R1.fastq.gz | cut -c1-8 > prefix.txt
```

6) run fastQC.sh 


7) cd to trim, qsub trimgalore.sh
8) cd to map, open bowtie2.sh check the genome is the correct file for your project
and that it exists
9) qsub bowtie2.sh
10) qsub rmdup.sh to remove duplicates
- This is optional but advised in most cases,
if you decide to keep duplicates, you may skip but you will have to edit
the following scripts so that they can find your bam files.

11) cd to peakcall, make two files:

"prefix_test.txt" - contains all prefixes for chip samples (no input)
        "prefix_ctrl.txt" - contains all prefixes for input controls (or igg),
        must be in the same order corresponding to the chip samples
        eg: C00005CL is the input for the first 2 samples and C00006CL is the
        input for the second 2

        prefix_test.txt:        prefix_ctrl.txt:
        C00001CL                C00005CL
        C00002CL                C00005CL
        C00003CL                C00006CL
        C00004CL                C00006CL

        If the number of lines in prefix_ctrl.txt is shorter than prefix_test.txt
        the array job will still run. For example If you have only 1 input for your whole experiment

12) qsub macs2.sh or macs2_broad.sh for broad peak calling

13) qsub makebigwigs.sh - for this you require a file "chrom.sizes.txt" containing the chromosome
sizes of your organism. You can find this online or create it from your bam
file header! Extract with samtools view and then use sed to remove the
unwanted characters e.g.:

samtools view -H R00202CL.bam | head -n -1 | sed 's/@SQ\t//g; s/SN://g; s/\t/ /g; s/LN://g' | tail -n +2 > chrom.sizes.txt

14) *Optional 1*  qsub filterpeaks.sh - for this you require a bed file
containing blacklisted regions (includes repetitive regions which cause
mapping errors and artefacts). You can find these in here:
https://sites.google.com/site/anshulkundaje/projects/blacklists
They may need to have the chromsome names edited to fit ensembl annotation (remove "chr").

15) *Optional 2*: Use bedtools to look for differences between peak files. You need to
decide which are the best options to use as this is project specific.

http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
eg.
bedtools intersect -wa -a C00012RK_peaks.broadPeak -b C00016RK_peaks.broadPeak > wa.12v16.broadPeak

http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html
eg.
bedtools subtract -A -a C00012RK_peaks.broadPeak -b C00016RK_peaks.broadPeak > 12v16.subtract.broadPeak

16) *Optional 3*: Make heatmaps using Deeptools!

17) Any other steps will now be completed in R using bioconductor packages
such as diffbind or csaw. This is highly project specific!
       

