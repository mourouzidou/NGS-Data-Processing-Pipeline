# Ensure entries.txt exists
if [ ! -f entries.txt ]; then
    echo "entries.txt not found in the current directory!"
    exit 1
fi

# Trimming with Sickle
for i in $(cat entries.txt); do
    sickle pe -f ERR0091${i}_1.fastq -r ERR0091${i}_2.fastq -t sanger -o trimmedERR0091${i}_1.fastq -p trimmedERR0091${i}_2.fastq -s ERR0091${i}trimmed_s.fastq -q 12 -l 15
done

# Run FastQC on trimmed files
find . -iname "tr*.fastq" | xargs -I {} fastqc {}

# Run MultiQC to aggregate FastQC reports
multiqc . || { echo "multiqc command not found!"; exit 1; }

# Index genome using STAR
STAR --runMode genomeGenerate --genomeDir ../ --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v43.chr_patch_hapl_scaff.annotation.gtf --genomeSAindexNbases 10 --sjdbOverhang 90 --runThreadN 10 || { echo "STAR command not found or failed!"; exit 1; }

# Mapping with STAR
for e in $(cat entries.txt); do
    STAR --genomeDir .. --readFilesIn ./ERR0091${e}trimmed_s.fastq --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $e || { echo "STAR mapping failed for $e"; exit 1; }
done

# Count features using HTSeq
ls *Aligned* | xargs -n1 -P6 -I {} sh -c 'htseq-count -s no -r pos -t exon -i pacid -f bam "$1" ../Homo_sapiens.GRCh38.109.chromosome.19.gff3 > "$1.counts"' -- {}

# Index BAM files
for i in $(ls *.bam); do
    samtools index "$i" || { echo "Failed to index BAM file $i"; exit 1; }
done

# Trimming again (if needed)
for i in $(cat entries.txt); do
    sickle pe -f ERR0091${i}_1.fastq -r ERR0091${i}_2.fastq -t sanger -o trimmedERR0091${i}_1.fastq -p trimmedERR0091${i}_2.fastq -s ERR0091${i}trimmed_s.fastq -q 12 -l 15
done

# Run FastQC on trimmed files (again, if needed)
find . -iname "tr*.fastq" | xargs -I {} fastqc {}

# Multiple quality control with MultiQC
./multiqc . || { echo "multiqc command not found!"; exit 1; }

# Index genome using STAR (repeated, ensure you only run it if necessary)
STAR --runMode genomeGenerate --genomeDir ../ --genomeFastaFiles GRCh38.p13.genome.fa --sjdbGTFfile gencode.v43.chr_patch_hapl_scaff.annotation.gtf --genomeSAindexNbases 10 --sjdbOverhang 90 --runThreadN 10 || { echo "STAR command not found or failed!"; exit 1; }

# Mapping with STAR (repeated, only if necessary)
for e in $(cat entries.txt); do
    STAR --genomeDir .. --readFilesIn ./ERR0091${e}trimmed_s.fastq --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $e || { echo "STAR mapping failed for $e"; exit 1; }
done

# HTSeq counting (repeated, only if necessary)
ls *Aligned* | xargs -n1 -P6 -I {} sh -c 'htseq-count -s no -r pos -t exon -i pacid -f bam "$1" ../Homo_sapiens.GRCh38.109.chromosome.19.gff3 > "$1.counts"' -- {}

# Index BAM files (repeated, only if necessary)
for i in $(ls *.bam); do
    samtools index "$i" || { echo "Failed to index BAM file $i"; exit 1; }
done
