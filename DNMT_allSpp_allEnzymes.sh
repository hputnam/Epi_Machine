#!/bin/bash
#SBATCH --job-name="Fellous_test"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery

echo "This script is to test the full blast pipeline using the DNMTs fasta file as query, to the CDS sequences of all coral custom sequence databases."
echo "To be executed on the URI Andromenda HPC. If using another server, the module version may need to be updated. Use module show to see latest version."
echo "Starting test..." $(date)

module load all/BLAST+/2.10.1-gompi-2020a #BLAST+
module list

echo "STARTING TEST $(date)"
echo -e "Creating array of Species Spp and setting variable resDir to easily place results in the directory blast_res/test_res" 
mapfile -t Spp <species #Make a string of queries for the query variable
resDir="blast_res/test_res"

echo "Running tblastn of query DNMTs against all species databases $(date)"
for Sp in "${!Spp[@]}"; do
	echo -e "Now conducting local alignments of protein queries to ${Spp[Sp]} nucleotides with tblastn. $(date)"
	tblastn -query queries/DNMTs -db db/"${Spp[Sp]}"/"${Spp[Sp]}" \
		-outfmt 11 -evalue 1e-05 -out "$resDir"/"${Spp[Sp]}"_DNMTs.asn
	echo "Now creating desired output files from tblastn search" $(date)
	blast_formatter -archive "$resDir"/"${Spp[Sp]}"_DNMTs.asn -outfmt 10 -out "$resDir"/"${Spp[Sp]}"_DNMTs.csv
	blast_formatter -archive "$resDir"/"${Spp[Sp]}"_DNMTs.asn -sorthits 1 -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out "$resDir"/"${Spp[Sp]}"_DNMTs.NOTfasta 
	awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' "$resDir"/"${Spp[Sp]}"_DNMTs.NOTfasta > "$resDir"/"${Spp[Sp]}"_DNMTs_p.fasta
	awk -F',' '{print $2}' "$resDir"/"${Spp[Sp]}"_DNMTs.NOTfasta > "$resDir"/"${Spp[Sp]}"_new_query.fasta
	awk 'FNR==NR{a[$0];next} /^>/{val=$0;sub(/^>/,"",val);flag=val in a?1:0} flag' queries/"${Spp[Sp]}"_blastx_query.fasta  ref/"${Spp[Sp]}".fasta
	rm "$resDir"/"${Spp[Sp]}"_DNMTs.NOTfasta
	blast_formatter -archive "$resDir"/"${Spp[Sp]}"_DNMTs.asn -sorthits 1 -max_target_seqs 1 -outfmt "6 qseqid sstart send sseqid" -out "$resDir"/"${Spp[Sp]}"_DNMTs_blastn.tsv
	while read -r QUERY START STOP SUBJECT; do 
	echo "$QUERY" "$START" "$STOP" "$SUBJECT"; blastdbcmd -entry "$QUERY" -db db/"${Spp[Sp]}"/"${Spp[Sp]}" -target_only -range "$START"-"$STOP" >> queries/"${Spp[Sp]}"_blastn_query.NOTfasta; 
	sed 's/>/"$QUERY"|"${Spp[Sp]}""$SUBJECT"/g' queries/"${Spp[Sp]}"_blastn_query.NOTfasta > "$resDir"/"${Spp[Sp]}"_DNMTs_n.fasta
	rm queries/"${Spp[Sp]}"_blastn_query.NOTfasta
done < "$resDir"/"${Spp[Sp]}"_DNMTs_blastn.tsv
done

cat > "$resDir"/*_DNMTs_n.fasta > DNMTs_allspp.nucl.fasta
cat > "$resDir"/*_DNMTs_p.fasta > DNMTs_allspp.prot.fasta
cat > "$resDir"/*_DNMTs.csv > DNMTs_allspp_alignment_report.csv

echo -e "\n Now running blastp with the output of the tblastn search $(date)"
for Sp in "${!Spp[@]}"; do
	echo "Running blastp on ${Spp[Sp]} versus nr database $(date)"
	mv "$resDir"/"${Spp[Sp]}"_DNMTs.fasta "$resDir"/"${Spp[Sp]}"_DNMTs1.fasta
	sed -e 's/-//g' "$resDir"/"${Spp[Sp]}"_DNMTs1.fasta > "$resDir"/"${Spp[Sp]}"_DNMTs.fasta
	rm "$resDir"/"${Spp[Sp]}"_DNMTs1.fasta
	blastp -query "$resDir"/"${Spp[Sp]}"_DNMTs.fasta -db nr -remote \
		-evalue 1e-05 -max_hsps 1 -culling_limit 1 -out "$resDir"/"${Spp[Sp]}"_DNMTs_nr_hits.out
	RID=($(grep RID "$resDir"/"${Spp[Sp]}"_DNMTs_nr_hits.out | head -1 | sed -e 's/RID: //g'))
	echo -e "\n Now creating desired output files for blastp search $(date)"
	blast_formatter -rid "$RID" -outfmt "10 qseqid sacc sscinames stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out "$resDir"/"${Spp[Sp]}"_DNMTs_nr_hits.csv
#	blast_formatter -rid "$RID" -sorthits 1 -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out "$resDir"/"${Spp[Sp]}"_DNMTs_nrScler_hits.NOTfasta
#	awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' "$resDir"/"${Spp[Sp]}"_DNMTs_nrScler_hits.NOTfasta > "$resDir"/"${Spp[Sp]}"_DNMTs_nrScler_hits.fasta
#	rm "$resDir"/"${Spp[Sp]}"_DNMTs_nrScler_hits.NOTfasta
done

cat "$resDir"/*_DNMTs_nr_hits.csv > allspp_nr_hits.csv
echo "TEST COMPLETE!" $(date)


