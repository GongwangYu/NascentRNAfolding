

#sample="test"


index="/mnt/data/home/yugongwang/disk3/Nascent_RNA/reference/yeast/yeast_bowtie2_full-gene_index/Saccharomyces_cerevisiae.R64-1-1.95.full-gene"
trimmomatic="java -jar /mnt/data/home/yugongwang/APP/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 40"

for sample in vivo_rep1 vivo_rep2 dmso_rep1 dmso_rep2
do
###trimed adapter
mkdir -p ${sample}
$trimmomatic ./test/${sample}/${sample}_1.fq  ./test/${sample}/${sample}_2.fq \
-baseout ./${sample}/trimedAdapter_${sample}.fq  \
ILLUMINACLIP:/mnt/data/home/yugongwang/APP/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:1:true \
MINLEN:21


### trimed barcode
cutadapt -u 3 -m 18 -j 40 -o ./${sample}/trimedBarcode_${sample}_1.fq ./${sample}/trimedAdapter_${sample}_1P.fq
cutadapt -u -3 -q 15,0 -m 18 -j 40  -o ./${sample}/trimedBarcode_${sample}_2.fq ./${sample}/trimedAdapter_${sample}_2P.fq


### mapped to full-gene
#read1

bowtie2 -x $index -U ./${sample}/trimedBarcode_${sample}_1.fq  -S ./${sample}/mapped_${sample}_1.sam \
--un ./${sample}/unMapped_${sample}_1.fq  -p 20 --norc  &>./${sample}/align_summary_${sample}_1.txt

#read2

bowtie2 -x $index -U ./${sample}/trimedBarcode_${sample}_2.fq  -S ./${sample}/mapped_${sample}_2.sam \
--un ./${sample}/unMapped_${sample}_2.fq -p 20 --nofw &>./${sample}/align_summary_${sample}_2.txt


###mapping quality filter
#read1
samtools view -@ 20 -F 260 -q 30 ./${sample}/mapped_${sample}_1.sam -o ./${sample}/mapped_${sample}_flt_1.sam

#read2
samtools view -@ 20 -F 260 -q 30 ./${sample}/mapped_${sample}_2.sam -o ./${sample}/mapped_${sample}_flt_2.sam



###remove PCR duplicate
python3 readCollapse.py ./test/${sample}/${sample}_1.fq ./${sample}/mapped_${sample}_flt_1.sam \
./${sample}/mapped_${sample}_flt_2.sam ./${sample}/rmdup_${sample}_1.sam ./${sample}/rmdup_${sample}_2.sam

###
less ./${sample}/rmdup_${sample}_1.sam|cut -f1-4,6 > ./${sample}/${sample}_1.txt
less ./${sample}/rmdup_${sample}_2.sam|cut -f1-4,6 > ./${sample}/${sample}_2.txt
done

