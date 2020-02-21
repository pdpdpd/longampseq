#!/bin/bash


#modified from 0915 verision
#end with fastq
#start with merging
#use full hg19
#use q30 only for the initial alignment and filtering (same for 0915)
#removed the deduct part

for file in *R1_001.fastq.gz
do
  flash -M600 $file ${file/R1/R2}
  mv ./out.extendedFrags.fastq ${file/.fastq.gz/_merged.fastq}
done

#chr index: F primer:	AGTGGGGCTGGAATAAAAGTAGAAT, Chr11:5250918-5250942 - strand
#chr index: F primer:	TTTTTCCTTTTGTTGCCTTTGCTTC, Chr11:5245453-5245477 + strand

#HBB5kb
#  samtools view -b -q 30 ${file/.fastq/_sorted.bam} chr11:5245453-5248229 > ${file/.fastq/_sortedleft.bam}
#  samtools index ${file/.fastq/_sortedleft.bam}
#  samtools view -b -q 30 ${file/.fastq/_sorted.bam} chr11:5248230-5250918 > ${file/.fastq/_sortedright.bam}
#done

echo 02/20/2020 > log.txt

for file in *merged.fastq
do
  bwa mem -A2 -E1 /home/yp11/Desktop/genomes/hg19/hg19.fa ${file} > ${file/.fastq/.sam}
  samtools view -S -b -q 30 ${file/.fastq/.sam} | samtools sort -o ${file/.fastq/_sorted.bam}
  samtools index ${file/.fastq/_sorted.bam}
  samtools view -b ${file/.fastq/_sorted.bam} chr2:60720359-60722401 > ${file/.fastq/_sortedleft.bam}
  samtools index ${file/.fastq/_sortedleft.bam}
  samtools view -b ${file/.fastq/_sorted.bam} chr2:60722402-60724626 > ${file/.fastq/_sortedright.bam}
  samtools index ${file/.fastq/_sortedright.bam}
  samtools view -F 4 ${file/.fastq/_sortedleft.bam} | cut -f1 | sort -u > ${file/.fastq/_leftID.txt}
  samtools view -F 4 ${file/.fastq/_sortedright.bam} | cut -f1 | sort -u > ${file/.fastq/_rightID.txt}
  comm -12 ${file/.fastq/_leftID.txt} ${file/.fastq/_rightID.txt} > ${file/.fastq/_bothID.txt}
  seqtk subseq ${file} ${file/.fastq/_bothID.txt} > ${file/.fastq/30_filtered.fastq}
#now we have filtered fastq

  bwa mem -A2 -E1 /home/yp11/Desktop/genomes/hg19/hg19.fa ${file/.fastq/30_filtered.fastq} >${file/.fastq/30_filtered.sam}
  samtools view -S -b ${file/.fastq/30_filtered.sam} -o ${file/.fastq/30_filtered.bam}
  samtools sort ${file/.fastq/30_filtered.bam} -o ${file/.fastq/30_filteredsorted.bam}
  bedtools bamtobed -i ${file/.fastq/30_filteredsorted.bam} > ${file/.fastq/30_filtered.bed}
#get the filtered, non-deducted reads

#extract the reads that have supplementary alignments
  samtools view -F 4 ${file/.fastq/30_filtered.bam} | cut -f1 | uniq -d > ${file/.fastq/filtered_2+alignID.txt}
  samtools view -F 4 ${file/.fastq/30_filtered.bam} | cut -f1 | uniq -u > ${file/.fastq/filtered_1alignID.txt}
  seqtk subseq ${file/.fastq/30_filtered.fastq} ${file/.fastq/filtered_2+alignID.txt} > ${file/.fastq/30_filtered_2+.fastq}

  bwa mem -A2 -E1 /home/yp11/Desktop/genomes/hg19/hg19.fa ${file/.fastq/30_filtered_2+.fastq} >${file/.fastq/filtered_2+.sam}
  samtools view -S -b ${file/.fastq/filtered_2+.sam} -o ${file/.fastq/filtered_2+.bam}
  bedtools bamtobed -i ${file/.fastq/filtered_2+.bam} > ${file/.fastq/filtered_2+.bed}

  #convert to csv
  cat ${file/.fastq/filtered_2+.bed} | tr "\\t" "," > ${file/.fastq/filtered_2+.csv}
  python ~/Documents/scripts/0220_bedfile.py ${file/.fastq/filtered_2+.csv} ${file/.fastq/largedel_output.csv} ${file/.fastq/largedel_group.csv}

  # just find the starting position and length of large deletion
  # ignore strand and reads with multiple fragments
  # the format of bed files are more like this:
  # chr11	59136655	59136956	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
  # chr11	59138619	59138657	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
  # call a python script for bed file analysis


  # append read numbers into the log file.
  # total aligned events:
  echo ${file/_merged.fastq/} >> log.txt
  echo $(cat ${file/.fastq/30_filtered.fastq}|wc -l)/4|bc >> log.txt
done