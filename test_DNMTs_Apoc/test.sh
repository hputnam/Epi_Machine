#!/bin/bash
#SBATCH --job-name="blastp"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery

echo "This script is to test the full blast pipeline using the DNMTs fasta file as query, to the Apoculata CDS sequences as a database."
echo "To be executed on the URI Andromenda HPC. If using another server, please double-check if the module version needs to be updated. Use 'module show' to see latest version."
echo "Starting test..." $(date)

module load all/BLAST+/2.10.1-gompi-2020a BLAST+
module list

echo "Now creating desired output files for blast search 1... to be used for STable 1, blastp, and phylogenetic analysis" $(date)
blast_formatter -archive blast_res/DNMTs/TEST_Apoc_DNMTs.asn -outfmt 10 -out blast_res/DNMTs/TEST_Apoc_DNMTs.csv
blast_formatter -archive blast_res/DNMTs/TEST_Apoc_DNMTs.asn -sorthits 1 -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out blast_res/DNMTs/TEST_Apoc_DNMTs.NOTfasta
awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' blast_res/DNMTs/TEST_Apoc_DNMTs.NOTfasta > blast_res/DNMTs/TEST_Apoc_DNMTs.fasta
rm blast_res/DNMTs/TEST_Apoc_DNMTs.NOTfasta

echo "Now running blastp with the output of the tblastn search" $(date)
mv blast_res/DNMTs/TEST_Apoc_DNMTs.fasta blast_res/DNMTs/TEST_Apoc_DNMTs1.fasta
sed -e 's/-//g' blast_res/DNMTs/TEST_Apoc_DNMTs1.fasta > blast_res/DNMTs/TEST_Apoc_DNMTs.fasta
rm blast_res/DNMTs/TEST_Apoc_DNMTs1.fasta
blastp -query blast_res/DNMTs/TEST_Apoc_DNMTs.fasta -db nr -remote \
	-outfmt 11 -evalue 1e-05 -max_hsps 1 -culling_limit 1 -entrez_query "Scleractinia [Organism]" \
	-out blast_res/DNMTs/TEST2_Apoc_DNMTs.asn

echo -e "\n Now running blastp with the output of the tblastn search" $(date)
mv blast_res/DNMTs/TEST_Apoc_DNMTs.fasta blast_res/DNMTs/TEST_Apoc_DNMTs1.fasta
sed -e 's/-//g' blast_res/DNMTs/TEST_Apoc_DNMTs1.fasta > blast_res/DNMTs/TEST_Apoc_DNMTs.fasta
rm blast_res/DNMTs/TEST_Apoc_DNMTs.NOTfasta
blastp -query blast_res/DNMTs/TEST_Apoc_DNMTs.fasta -db nr -remote \
	-evalue 1e-05 -max_hsps 1 -culling_limit 1 -entrez_query "Scleractinia [Organism]" -out blast_res/DNMTs/TEST2_Apoc_DNMTs.out
RID=($(grep RID blast_res/DNMTs/TEST2_Apoc_DNMTs.out | head -1 | sed -e 's/RID: //g'))

echo -e "\n Now creating desired output files for blast search 2... to be used for STable 2 and phylogenetic analysis" $(date)
blast_formatter -rid "$RID" -outfmt "10 qseqid sacc stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out blast_res/DNMTs/TEST2_Apoc_DNMTs.csv

blast_formatter -rid "$RID" -sorthits 1 -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out blast_res/DNMTs/TEST2_Apoc_DNMTs.NOTfasta
awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' blast_res/DNMTs/TEST2_Apoc_DNMTs.NOTfasta > blast_res/DNMTs/TEST2_Apoc_DNMTs.fasta
rm blast_res/DNMTs/TEST2_Apoc_DNMTs.NOTfasta

echo "Test complete!" $(date)
