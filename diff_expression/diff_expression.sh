
###################
####Count reads####
###################

## Set variables
bam_dir=mypath/Apis_RNA_seq/output/STAR/Sorted/
htseq_count_output_dir=mypath/Apis_RNA_seq/output/htseq_count/
gtf_file=/public/home/shipy3/DB/Amel/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf

## Count reads by htseq-count
for i in `seq 571 584`;
do
        mypath/miniconda3/bin/htseq-count -f bam -r pos -s no \
        -a 10 -t exon -i gene_id -m union \
        ${bam_dir}${i}.sorted.bam \
        ${gtf_file} \
        > ${htseq_count_output_dir}${i}.count ;
done


bam_dir=mypath/Apis_RNA_seq/output/STAR/Sorted/
htseq_count_output_dir=mypath/Apis_RNA_seq/output/htseq_count/
gtf_file=/public/home/shipy3/DB/Amel/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf


