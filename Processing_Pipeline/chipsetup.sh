#!/bin/bash
#$ -cwd

# get the number of samples to process
number=`wc -l prefix.txt | cut -d ' ' -f1`

# get the geneome and organism options
while getopts g:o: option
do
case "${option}"
in
g) genome=${OPTARG};;
o) org=${OPTARG};;

esac
done

mkdir fastqc
mkdir trim
mkdir map
mkdir peakcall

echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -N fastqc

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

fastqc $myFile\_R1.fastq.gz
fastqc $myFile\_R2.fastq.gz

exit 0;' > fastqc/fastqc.sh

echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid fastqc
#$ -N trimgalore

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

trim_galore --fastqc --paired --retain_unpaired --quality 20 \
-stringency 1  --gzip -e 0.05  \
../fastq/$myFile\_R1.fastq.gz ../fastq/$myFile\_R2.fastq.gz

exit 0;' > trim/trimgalore.sh


echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -pe smp 8
#$ -hold_jid trimgalore
#$ -N bowtie2

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`

#So, the PBS ID corresponds to that line in the files

#check is this correct?
genome='$genome'


bowtie2 -X 400 -N 1 -p 8 --un $myFile\_unaligned.fq -x $genome  \
		-1 ../trim/$myFile\_R1_val_1.fq.gz \
		-2 ../trim/$myFile\_R2_val_2.fq.gz 2>> $myFile\_res.summary \
			| samtools view -F 0xC -bS - \
			| samtools sort -T /tmp/$myFile.sorted -o $myFile.sorted.bam - \
				samtools index $myFile.sorted.bam

gzip $myFile\_unaligned.fq


exit 0;' > map/bowtie2.sh

echo '#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid bowtie2
#$ -N rmdup


myFile=`cat ../prefix.txt`

myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the PBS ID corresponds to that line in the files

java -Xmx512m -XX:ParallelGCThreads=1 -jar ~/groupso/bin/picard.jar MarkDuplicates \
											M=$myFile.2.met \
											I=$myFile.sorted.bam \
											O=$myFile.rmdup.bam \
											REMOVE_DUPLICATES=true  \
											VALIDATION_STRINGENCY=SILENT  \
											MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000


exit 0;' > map/rmdup.sh


echo '#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid rmdup
#$ -N macs2

myFile=`cat prefix_test.txt`
myCtrl=`cat prefix_ctrl.txt`
myFile=`head -n $SGE_TASK_ID prefix.txt | tail -1`
myCtrl=`head -n $SGE_TASK_ID prefix_ctrl.txt | tail -1`
#So, the PBS ID corresponds to that line in the files

# change -g mm to -g hs if human

macs2 callpeak -B -t ../map/$myFile.rmdup.bam \
				-c ../map/$myCtrl.rmdup.bam \ 
				-f BAMPE -g '$org' -n $myFile --keep-dup=all

exit 0;' > peakcall/macs2.sh


echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid rmdup
#$ -N macs2

myFile=`cat prefix_test.txt`
myCtrl=`cat prefix_ctrl.txt`
myFile=`head -n $SGE_TASK_ID prefix.txt | tail -1`
myCtrl=`head -n $SGE_TASK_ID prefix_ctrl.txt | tail -1`
#So, the PBS ID corresponds to that line in the files

macs2 callpeak -B -t ../map/$myFile.rmdup.bam \
				-c ../map/$myCtrl.rmdup.bam  \
				-f BAMPE -g '$org' --broad -q 0.05 \
				--broad-cutoff 0.1 -n $myFile --keep-dup=all

exit 0;' > peakcall/macs2_broad.sh


echo \
'#!/bin/bash
#$ -cwd
#$ -V
#$ -t 1-'$number'
#$ -l h_vmem=2G
#$ -N bigwigs


myFile=`cat prefix_test.txt`
myFile=`head -n $SGE_TASK_ID prefix_test.txt | tail -1`
#So, the PBS ID corresponds to that line in the files

bedGraphToBigWig  $myFile\_treat_pileup.bdg chrom.sizes.txt $myFile\_treat_pileup.bw
bedGraphToBigWig  $myFile\_control_lambda.bdg chrom.sizes.txt $myFile\_control_lambda.bw

exit 0;' > peakcall/makebigwigs.sh

