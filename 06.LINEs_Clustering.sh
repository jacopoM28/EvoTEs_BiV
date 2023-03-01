#Cluster all LINE-derived RVT in based on their nucleotide sequence. Then build up a protein consensus sequences for each cluster with more than 5 members.
#These consensus will be used for phylogenetic inference

cd ALL_LINES

#-----------------------------------------------------------Cluster elements-----------------------------------------------------------#

for i in *.fasta; do 
  cd-hit-est -g 1 -o "$i"_nr -i "$i" -n 4 -c 0.80 -t 1 -aS 0.8 -d 0 -G 0; 
done;

#-----------------------------------sort clusters by size and retrive corresponding sequences (min cluster size = 5)-------------------#

for i in *.fasta; do 
  varSpecie=$( echo "$i" | cut -d"_" -f1); 
  clstr_sort_by.pl < "$varSpecie"_LINEs_RVT.ORF.fasta_nr.clstr no > "$varSpecie"_LINEs_RVT.ORF.fasta_nr.ordered.clstr; 
  make_multi_seq.pl "$varSpecie"_LINEs_RVT.ORF.pep "$varSpecie"_LINEs_RVT.ORF.fasta_nr.ordered.clstr multi-seq 5; 
  mv multi-seq "$varSpecie"_Min5_Clusters;
  make_multi_seq.pl "$varSpecie"_LINEs_RVT.ORF.fasta "$varSpecie"_LINEs_RVT.ORF.fasta_nr.ordered.clstr multi-seq 1;
  mv multi-seq "$varSpecie"_nucl_ALL_Clusters;
  cd "$varSpecie"_Min5_Clusters; 
  for j in $( ls ); do 
    mv "$j" "$varSpecie"_"$j".fasta; 
  done; 
  cd ../; 
  cd "$varSpecie"_nucl_ALL_Clusters;
  for j in $( ls ); do 
    mv "$j" "$varSpecie"_"$j"_nucl.fasta;
  done;
  cd ../
done;

#-----------------------Retrive RVT segment domain for each sequence, align, trim them and finally build up a consensus--------------------#

for i in $( ls | grep "_Min5_Clusters"); do 
  cd "$i"; 
  for j in *fasta; do
    rpsblast -query "$j" -evalue 1e-5 -db /media/storage/jacopomartelossi/dbs/Pfam_rpsblast/Pfam -outfmt 6 -out "$j".rpsblast -num_threads 30
    awk -v OFS=" " -F"\t" '{print$1,$2";"$11,$12,$7"-"$8}' "$j".rpsblast  > "$j"_reformatted.rpsblast; 
    /media/storage/jacopomartelossi/Software/cath-tools/cath-resolve-hits.ubuntu-20.04 --input-format raw_with_scores "$j"_reformatted.rpsblast > "$j"_ResolvedOverlap.rpsblast; 
  done; 
  cd ../; 
done;

for i in $( ls | grep "_Min5_Clusters"); do
  cd "$i";
  for j in *_ResolvedOverlap.rpsblast; do 
    varFAM=$( echo "$j" | cut -d"." -f1,2); 
    grep -Ff /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/RVT_Simplified_CDD.txt "$j" | awk -v OFS="\t" '{print$1,$4}' | awk 'BEGIN{FS=OFS="\t"} {gsub(/-/, "\t", $2)} 1' > "$varFAM"_RVT.bed; 
  done;
  for k in *.bed; do 
    varName=$( echo "$k" | cut -d"_" -f1,2); 
    bedtools getfasta -bed "$k" -fi "$varName".fasta > "$varName"_RVT.fasta;
    mafft --localpair --maxiterate 1000 --thread 30 --anysymbol "$varName"_RVT.fasta > "$varName"_RVT.mafft;
    trimal -in "$varName"_RVT.mafft  -out "$varName"_RVT.mafft.trimal -gt 0.5; 
    cons -sequence "$varName"_RVT.mafft.trimal -plurality 3 -outseq "$varName"_cons.fasta -name "$varName"_LINE; 
  done;
cd ../
done;
