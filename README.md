# HIFI Q30 SSU Sequences

## Q30 reads

### Extracting reads from ```BAM``` files

```bash
ccs filename.bam --min-length 1000 --min-rq 0.999 -j 24 filename.ccs.1000.bam
samtools fasta -@ 20 filename.ccs.1000.bam > ft5_mf2.hifi_q30.fa
```

### Extracting reads from ```CCS.BAM``` files

```bash
samtools view -e '[rq]>=0.999' -@ 20 ft5_mf2.ccs.all.bam \
    | samtools fasta -@ 20 - \
    | seqkit --threads 20 seq -m 1000 > ft5_mf2.hifi_q30.fa
```

### Extracting reads from ```FASTQ``` files
```bash
ls *fastq | sed 's/.fastq//' 、
    | rush -j 12 'fastp -i {1} -o {1}.hifi_q30.fq.gz \
        --qualified_quality_phred 30 --unqualified_percent_limit 1 --length_required 1000 --thread 4 --html {1}.fastp_report.html --json {1}.fastp_report.json &> {1}.fastp.log' &> fastp.log
```

## Reference loci DB
homology-based strategy to find related reads. First, constructing a reference blastn DB
```bash

read -r -d '' URLS <<'EOF'
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.28SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.23SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.5SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.23SrRNA.fna.gz
https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.5SrRNA.fna.gz
https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/bac120_ssu_reps_r226.fna.gz
https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/ar53_ssu_reps_r226.fna.gz
EOF

echo "$URLS" | tr ' ' '\n' | while read -r u; do 
    wget -c "$u"
done

cat ar53_ssu_reps_r226.fna | awk '{print $1}' | awk '{ if ($1 ~ />/) {print $1"_ar53"} else {print }}' > ar53_simple_id.16s.fna

cat bac120_ssu_reps_r226.fna | awk '{print $1}' | awk '{ if ($1 ~ />/) {print $1"_bac120"} else {print }}' > bac120_simple_id.16s.fna

ls *fna | grep -v -e ar53 -e bac120 | sed 's/.fna//' | while read a;do cat ${a}.fna | awk '{print $1}' | awk -v loci="${a}" '{ if ($1 ~ />/) {print $1"_"loci} else {print }}';done > ncbi.archaea.fungi.bacteria.loci.fna

cat ncbi.archaea.fungi.bacteria.loci.fna ar53_simple_id.16s.fna bac120_simple_id.16s.fna >  fungi.bacteria.archaea.all.loci.fna

ls *fna | grep -v fungi.bacteria.archaea.all.loci.fna | xargs -I {} rm {}
cd-hit -i fungi.bacteria.archaea.all.loci.fna -o fungi.bacteria.archaea.all.loci.0.97.fna -c 0.97 -T 12

sed -i 's/bac120/bacteria120/' fungi.bacteria.archaea.all.loci.0.97.fna
sed -i 's/ar53/archaea53/' fungi.bacteria.archaea.all.loci.0.97.fna
sed -i 's/bacteria120/bacteria.16SrRNA/' fungi.bacteria.archaea.all.loci.0.97.fna
sed -i 's/archaea53/archaea.16SrRNA/' fungi.bacteria.archaea.all.loci.0.97.fna

makeblastdb -in fungi.bacteria.archaea.all.loci.0.97.fna -out fungi.bacteria.archaea.all.loci.0.97.DB -dbtype nucl
```



