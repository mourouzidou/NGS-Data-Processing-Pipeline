list=$1
reference=$2
for entry in $(cat $list);
do
    # Setting variable "var" to get every entry in the given file with NA..
    # Using awk to keep only the entries of interest and print the column $1 and Base Pair Counts to a new file for each individual 
    
    awk -v var=$entry -F'\t' ' { if ($10 =='var' )  print$--NF"\t\t "$1 } ' 20130502.phase3.analysis.sequence.index > ${entry}_allSequences.txt

    # Find the greatest number of BP for each individual ..| using xargs match this number with its ftp link
    # |xargs create a file containing the ftp links of interest
    awk  'BEGIN{max=0}{if(($1)>max)  max=($1)}END {print max}' ${entry}_allSequences.txt  | xargs -I {} awk -v var={}  '{ if ($1=='var') print$2}' ${entry}_allSequences.txt | xargs -I {} echo {} >> ftplinks.txt
    
    # Download the links above
    

done

for link in $(cat ftplinks.txt);
do

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/$link

done



# index the reference genome
#bwa index $reference

# create a file with SRR/ERR IDs
sed -E 's/.*read\/([A-Z]{1}RR[0-9]+)_.*/\1/' ftplinks.txt > SRRlist.txt

# mapping 
for id in $(cat SRRlist.txt);
do

    bwa mem -t 10 $reference ${id}_1.filt.fastq.gz ${id}_2.filt.fastq.gz  -o $id.sam
    samtools fixmate -O bam ${id}.sam ${id}.bam
    samtools sort -o ${id}.sort.bam $id.bam
    bcftools mpileup -g 10 -Oz -o ${id}.gvcf.gz -f ${reference} ${id}.sort.bam
    bcftools index ${id}.gvcf.gz
done

 
bcftools merge -Oz --gvcf $reference --merge all -o merged.vcf.gz *gvcf.gz
bcftools call -mv merged.vcf.gz -o merged.vcf.gz.mv.vcf
