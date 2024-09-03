list=$1

mkdir RnaERRs
cd RnaERRs

for entry in $(cat ../$list);
do
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR009/ERR0091${entry}/ERR0091${entry}_1.fastq.gz 
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR009/ERR0091${entry}/ERR0091${entry}_2.fastq.gz 

done

gunzip *.gz

fastqc ERR*.fastq
