###############
#### rMATS ####
###############

###Step 1: make rMATS comparison text files with paths to sample bam Files
touch F-An.txt
touch D-An.txt
touch NE-An.txt
touch N-An.txt
touch OLD-An.txt
echo mypath/Apis_RNA_seq/output/STAR/sorted/571.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/572.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/573.sorted.bam > mypath/bee/diffgroup/F-An.txt
echo mypath/Apis_RNA_seq/output/STAR/sorted/576.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/577.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/578.sorted.bam > mypath/bee/diffgroup/D-An.txt
echo mypath/Apis_RNA_seq/output/STAR/sorted/574.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/575.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/584.sorted.bam > mypath/bee/diffgroup/NE-An.txt
echo mypath/Apis_RNA_seq/output/STAR/sorted/581.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/582.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/583.sorted.bam > mypath/bee/diffgroup/N-An.txt
echo mypath/Apis_RNA_seq/output/STAR/sorted/571.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/576.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/581.sorted.bam > mypath/bee/diffgroup/OLD-An.txt

### step 2: Run rMATS
#f-d
python mypath/anaconda3/envsmypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/F-An.txt \
--b2 mypath/bee/diffgroup/D-An.txt \
--gtf mypath/project/Apis_OSN_bulk_RNA-seq/00ref/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/F-D \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/F-D
#n-f
nohup
python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/N-An.txt \
--b2 mypath/bee/diffgroup/F-An.txt  \
--gtf mypath/ref/bee/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/N-F \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/N-F &
#n-d
python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/N-An.txt \
--b2 mypath/bee/diffgroup/D-An.txt \
--gtf mypath/project/Apis_OSN_bulk_RNA-seq/00ref/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/N-D  \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/N-D

#ne-n
python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/NE-An.txt \
--b2 mypath/bee/diffgroup/N-An.txt \
--gtf mypath/project/Apis_OSN_bulk_RNA-seq/00ref/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/NE-N \
-t paired \
--variable-read-length \
--readLength 101 \
--nthread 10 \
--tmp mypath/bee/result/NE-N
#ne-f
nohup python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/NE-An.txt \
--b2 mypath/bee/diffgroup/F-An.txt \
--gtf mypath/project/Apis_OSN_bulk_RNA-seq/00ref/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/NE-F \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/NE-F &
#ne-d
nohup python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/NE-An.txt \
--b2 mypath/bee/diffgroup/D-An.txt \
--gtf mypath/project/Apis_OSN_bulk_RNA-seq/00ref/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/NE-D \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/NE-D &
#ne-old
nohup python mypath/anaconda3/envs/your_env_name/bin/rmats.py --b1 mypath/bee/diffgroup/NE-An.txt \
--b2 mypath/bee/diffgroup/OLD-An.txt \
--gtf mypath/ref/bee/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf \
--od mypath/bee/result/NEOLD \
-t paired \
--readLength 101 --variable-read-length \
--nthread 10 \
--tmp mypath/bee/result/NEOLD &

###########################
#### rmats2sashimiplot ####
###########################

### step1: merge bam fiels of replicates into one bam file
nohup /public/home/jinxu/software/samtools-1.10/samtools merge mypath/Apis_RNA_seq/output/STAR/sorted/NE.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/574.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/575.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/584.sorted.bam &
nohup /public/home/jinxu/software/samtools-1.10/samtools merge mypath/Apis_RNA_seq/output/STAR/sorted/F.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/571.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/572.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/573.sorted.bam &
nohup /public/home/jinxu/software/samtools-1.10/samtools merge mypath/Apis_RNA_seq/output/STAR/sorted/D.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/576.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/577.sorted.bam \
mypath/Apis_RNA_seq/output/STAR/sorted/578.sorted.bam &
nohup /public/home/jinxu/software/samtools-1.10/samtools merge mypath/Apis_RNA_seq/output/STAR/sorted/N.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/581.sorted.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/582.sorted.bam  \
mypath/Apis_RNA_seq/output/STAR/sorted/583.sorted.bam &

###step2: change 'group' to 'chr' in bam files and input file
nohup    /public/home/jinxu/software/samtools-1.10/samtools view -H mypath/Apis_RNA_seq/output/STAR/sorted/NE.bam | sed -e 's/SN:Group\([0-9]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' -e 's/SN:GroupUN_\([0-9]\)/SN: chrUN_\1/'| samtools reheader - mypath/Apis_RNA_seq/output/STAR/sorted/NE.bam > mypath/Apis_RNA_seq/output/STAR/sorted/NE.chr.bam &
nohup    /public/home/jinxu/software/samtools-1.10/samtools view -H mypath/Apis_RNA_seq/output/STAR/sorted/N.bam | sed -e 's/SN:Group\([0-9]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' -e 's/SN:GroupUN_\([0-9]\)/SN: chrUN_\1/'| samtools reheader - mypath/Apis_RNA_seq/output/STAR/sorted/N.bam > mypath/Apis_RNA_seq/output/STAR/sorted/N.chr.bam &
nohup    /public/home/jinxu/software/samtools-1.10/samtools view -H mypath/Apis_RNA_seq/output/STAR/sorted/F.bam | sed -e 's/SN:Group\([0-9]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' -e 's/SN:GroupUN_\([0-9]\)/SN: chrUN_\1/'| samtools reheader - mypath/Apis_RNA_seq/output/STAR/sorted/F.bam > mypath/Apis_RNA_seq/output/STAR/sorted/F.chr.bam &
nohup   /public/home/jinxu/software/samtools-1.10/samtools view -H mypath/Apis_RNA_seq/output/STAR/sorted/D.bam | sed -e 's/SN:Group\([0-9]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' -e 's/SN:GroupUN_\([0-9]\)/SN: chrUN_\1/'| samtools reheader - mypath/Apis_RNA_seq/output/STAR/sorted/D.bam > mypath/Apis_RNA_seq/output/STAR/sorted/D.chr.bam &


sed -e 's/chrGroup/chr/g' mypath/bee/result/NEN_typical.txt > mypath/bee/result/NEN_typical_chr.txt
sed -e 's/chrGroup/chr/g' mypath/bee/result/FD_typical.txt > mypath/bee/result/FD_typical_chr.txt

###step3: run rmats2sashimiplot
# NEN
nohup mypath/anaconda3/bin/rmats2sashimiplot --b1 mypath/Apis_RNA_seq/output/STAR/sorted/NE.chr.bam \
--b2 mypath/Apis_RNA_seq/output/STAR/sorted/N.chr.bam \
 -t SE \
 -e mypath/bee/result/NEN_typical_chr.txt \
 --exon_s 1 --intron_s 5 --l1 NE --l2 N \
 -o mypath/bee/result/rmatsplot_all/ATP_NEN &


# FD
nohup mypath/anaconda3/bin/rmats2sashimiplot --b1 mypath/Apis_RNA_seq/output/STAR/sorted/F.chr.bam \
--b2 mypath/Apis_RNA_seq/output/STAR/sorted/D.chr.bam \
 -t SE \
 -e mypath/bee/result/FD_typical_chr.txt \
 --exon_s 1 --intron_s 5 --l1 F --l2 D \
 -o mypath/bee/result/rmatsplot_all/ATP_FD &

###########################
#### isoforms analysis ####
###########################

### set variables
NE=mypath/Apis_RNA_seq/output/STAR/sorted/574.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/575.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/584.sorted.bam
N=mypath/Apis_RNA_seq/output/STAR/sorted/581.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/582.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/583.sorted.bam
F=mypath/Apis_RNA_seq/output/STAR/sorted/571.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/572.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/573.sorted.bam
D=mypath/Apis_RNA_seq/output/STAR/sorted/576.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/577.sorted.bam,mypath/Apis_RNA_seq/output/STAR/sorted/578.sorted.bam

### cufflinks install
#conda install -c bioconda cufflinks

### run cufflinks
nohup cuffdiff -o mypath/bee/isoform/refall/NEN -p 7 --min-reps-for-js-test 3 mypath/ref/bee/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf $NE $N &
nohup cuffdiff -o mypath/bee/isoform/refall/FD -p 7 --min-reps-for-js-test 3 mypath/ref/bee/Amel_HAv3.1/annotation/GCF_003254395.2_Amel_HAv3.1_genomic.gtf $F $D &

