#############################################
#### Trim adaptors and low quality reads ####
#############################################

##Set variables
input_dir=mypath/Apis_RNA_seq/input/
output_dir=mypath/Apis_RNA_seq/output/trimmed_fastqc/


for i in `seq 571 578`;
do
	java -jar mypath/biosoft/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33  \
	${input_dir}DP8400013540BL_L01_${i}_1.fq.gz ${input_dir}DP8400013540BL_L01_${i}_2.fq.gz \
	${input_dir}trimmed_${i}_1P.fastq.gz ${input_dir}trimmed_${i}_1U.fastq.gz \
	${input_dir}trimmed_${i}_2P.fastq.gz ${input_dir}trimmed_${i}_2U.fastq.gz \
	ILLUMINACLIP:mypath/biosoft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
	SLIDINGWINDOW:5:20 LEADING:10 TRAILING:10 MINLEN:50 ;
done

for i in `seq 581 584`;
do
	java -jar mypath/biosoft/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33  \
	${input_dir}DP8400013540BL_L01_${i}_1.fq.gz ${input_dir}DP8400013540BL_L01_${i}_2.fq.gz \
	${input_dir}trimmed_${i}_1P.fastq.gz ${input_dir}trimmed_${i}_1U.fastq.gz \
	${input_dir}trimmed_${i}_2P.fastq.gz ${input_dir}trimmed_${i}_2U.fastq.gz \
	ILLUMINACLIP:mypath/biosoft/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
	SLIDINGWINDOW:5:20 LEADING:10 TRAILING:10 MINLEN:50 ;
done
#######################################
#### Quality control after trimming####
#######################################

for i in `seq 571 578`;
do
	fastqc -o ${output_dir} -t 2 ${input_dir}trimmed_${i}_1P.fastq.gz ${input_dir}trimmed_${i}_2P.fastq.gz
done
## Set variables
trimmed_fastq_dir=mypath/project/Apis_OSN_bulk_RNA-seq/input/
output_dir=mypath/Apis_RNA_seq/output/STAR/

for i in `seq 581 584`;
do
	fastqc -o ${output_dir} -t 2 ${input_dir}trimmed_${i}_1P.fastq.gz ${input_dir}trimmed_${i}_2P.fastq.gz
done
## Set variables
trimmed_fastq_dir=mypath/project/Apis_OSN_bulk_RNA-seq/input/
output_dir=mypath/Apis_RNA_seq/output/STAR/
#####################################
#### Mapping reads to the genome ####
#####################################

for i in `seq 571 578`;
do
  /public/home/jinxu/software/STAR-master/bin/Linux_x86_64/STAR --runThreadN 8 \
  --genomeDir /public/home/shipy3/DB/Amel/Amel_HAv3.1/STAR_index \
  --readFilesIn ${trimmed_fastq_dir}trimmed_${i}_1P.fastq.gz ${trimmed_fastq_dir}trimmed_${i}_2P.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix ${output_dir}${i} \
  --outSAMtype BAM Unsorted
done

for i in `seq 581 584`;
do
  /public/home/jinxu/software/STAR-master/bin/Linux_x86_64/STAR --runThreadN 8 \
  --genomeDir /public/home/shipy3/DB/Amel/Amel_HAv3.1/STAR_index \
  --readFilesIn ${trimmed_fastq_dir}trimmed_${i}_1P.fastq.gz ${trimmed_fastq_dir}trimmed_${i}_2P.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix ${output_dir}${i} \
  --outSAMtype BAM Unsorted
done

#####################
####samtools sort####
#####################

## Set variables
input_dir=mypath/Apis_RNA_seq/output/STAR/
output_dir=mypath/Apis_RNA_seq/output/STAR/Sorted/

## Sorted by samtools
for i in `seq 571 578`;
do
   samtools sort -@ 20 -o ${output_dir}${i}.sorted.bam \
   ${input_dir}${i}Aligned.out.bam
done

for i in `seq 581 584`;
do
   samtools sort -@ 20 -o ${output_dir}${i}.sorted.bam \
   ${input_dir}${i}Aligned.out.bam
done