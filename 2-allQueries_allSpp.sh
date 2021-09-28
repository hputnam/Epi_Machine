#!/bin/bash
#SBATCH --job-name="Fellous_test"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=erin_chille@uri.edu
#SBATCH --exclusive
#SBATCH -D /data/putnamlab/erin_chille/Fellous_Epi_Machinery

echo "JOB STARTED $(date)"
echo "This script is to test the full blast pipeline using all queries to the CDS sequences of all coral custom sequence databases."
echo "To be executed on the URI Andromenda or Bluewaves HPC. If using another server, the module version may need to be updated. Use module show to see latest version."

module load all/BLAST+/2.10.1-gompi-2020a #BLAST+
module list

echo -e "Creating array of Species Spp and setting variable resDir to easily place results in the directory blast_res/test_res" 
mapfile -t Spp <species #Make a string of queries for the query variable
resDir="blast_res" #CHANGE THIS TO WHEREEVER YOU WANT THE RESULTS TO LIVE

echo "Showing loop combinations $(date)"

for Sp in "${!Spp[@]}"; do
        mapfile -t Qrys <query
        for Qry in "${!Qrys[@]}"; do
                printf '%s  %s\n' "${Spp[Sp]}" "${Qrys[Qry]}"
        done
done

echo "Running tblastn of query "${Qrys[Qry]}" against all species databases. \
Output will inlcude nucleotide and protein fasta files and an alignment report for each combination $(date)"

for Sp in "${!Spp[@]}"; do
	mapfile -t Qrys <query
	for Qry in "${!Qrys[@]}"; do
		echo -e "Now conducting local alignments of protein queries to ${Spp[Sp]} nucleotides with tblastn. $(date)"
		tblastn -query queries/"${Qrys[Qry]}" -db db/"${Spp[Sp]}"/"${Spp[Sp]}" \
			-outfmt 11 -evalue 1e-05 -out "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".asn
		echo "Now creating desired output files from tblastn search" $(date)
		blast_formatter -archive "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".asn -outfmt 10 -out "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".csv
		blast_formatter -archive "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".asn -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".NOTfasta 
		awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".NOTfasta > "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}"_p.fasta
		rm "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".NOTfasta
		blast_formatter -archive "$resDir"/"${Spp[Sp]}"_"${Qrys[Qry]}".asn -max_target_seqs 1 -outfmt "6 qseqid sstart send sseqid" -out "$resDir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_blastn.tsv
		while read -r QUERY START STOP SUBJECT; do 
			echo "Extracting $QUERY for  $SUBJECT $(date)"
			blastdbcmd -entry "$SUBJECT" -db db/"${Spp[Sp]}"/"${Spp[Sp]}" -target_only -range "$START"-"$STOP" >> queries/"${Spp[Sp]}"_blastn_query.NOTfasta 
			sed "s/>/>$QUERY|${Spp[Sp]}/g" queries/"${Spp[Sp]}"_blastn_query.NOTfasta >> "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_n.fasta
			rm queries/$"${Spp[Sp]}"_blastn_query.NOTfasta
		done < "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_blastn.tsv
	done
done

echo "Combining results $(date)"
mapfile -t Qrys <query
for Qry in "${!Qrys[@]}"; do
	cat "$resDir"/"${Qrys[Qry]}"/*_"${Qrys[Qry]}"_n.fasta > "$resDir"/"${Qrys[Qry]}"_allspp.nucl.fasta
	cat "$resDir"/"${Qrys[Qry]}"/*_"${Qrys[Qry]}"_p.fasta > "$resDir"/"${Qrys[Qry]}"_allspp.prot.fasta
	cat "$resDir"/"${Qrys[Qry]}"/*_"${Qrys[Qry]}".csv > "$resDir"/"${Qrys[Qry]}"_allspp_alignment_report.csv
	cat "$resDir"/"${Qrys[Qry]}"_allspp.nucl.fasta >> "$resDir"/allqueries_allspp.nucl.fasta
	cat "$resDir"/"${Qrys[Qry]}"_allspp.prot.fasta >> "$resDir"/allqueries_allspp.prot.fasta
	cat "$resDir"/"${Qrys[Qry]}"_allspp_alignment_report.csv >> "$resDir"/allqueries_allspp_alignment_report.csv
done

echo -e "\n Now running a blastp search of all queries against the nr database with the output of the tblastn search  \
Output will inlcude a protein fasta file and an alignment report for each combination $(date)"

for Sp in "${!Spp[@]}"; do
	mapfile -t Qrys <query
	for Qry in "${!Qrys[@]}"; do
		echo "Running blastp on ${Spp[Sp]} versus nr database $(date)"
		mv "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_p.fasta "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"1.fasta
		sed -e 's/-//g' "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"1.fasta > "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_p.fasta
		rm "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"1.fasta
		blastp -query "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_p.fasta -db nr -remote \
			-evalue 1e-05 -max_hsps 1 -culling_limit 1 -entrez_query "Scleractinia [Organism] OR Porifera [Organism]" -out "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.out
		RID=($(grep RID "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.out | head -1 | sed -e 's/RID: //g'))
		echo -e "\n Now creating desired output files for blastp search $(date)"
		blast_formatter -rid "$RID" -outfmt "10 qseqid sacc sscinames stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.csv
		blast_formatter -rid "$RID" -max_target_seqs 1 -outfmt "6 qseqid sseqid sseq" -out "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.NOTfasta
		awk 'BEGIN { OFS = "\n" } { print ">"$1"|"$2, $3 }' "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.NOTfasta > "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.fasta
		rm "$resdir"/"${Qrys[Qry]}"/"${Spp[Sp]}"_"${Qrys[Qry]}"_nr_hits.NOTfasta
	done
done

echo "Combining results $(date)"
for Qry in "${!Qrys[@]}"; do
	cat "$resDir"/"${Qrys[Qry]}"/*_"${Qrys[Qry]}"_nr_hits.csv > "$resDir"/"${Qrys[Qry]}"_allspp_nr_hits.csv
	cat "$resDir"/"${Qrys[Qry]}"/*_"${Qrys[Qry]}"_nr_hits.fasta > "$resDir"/"${Qrys[Qry]}"_allspp_nr_hits.fasta
	cat "$resDir"/"${Qrys[Qry]}"_allspp_nr_hits.csv >> "$resDir"/allqueries_allspp_nr_hits.csv
	cat "$resDir"/"${Qrys[Qry]}"_allspp_nr_hits.fasta >> "$resDir"/allqueries_allspp_nr_hits.fasta
done

echo -e "Blastp complete /n JOB COMPLETE! $(date)"
