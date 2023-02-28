#!bin/bash
#Reformat merged species-specific repeat libraries, apply a last copy-number based filter and run RepeatMasker

#-------------------------------Rename MITEs elements and Helitrons-----------------------------------------------#

for i in $( ls | grep "_TransposableElements"); do 
  varSpecie=$(echo "$i" | cut -d"_" -f1); 
  cd "$i"/Libraries; sed 's/|/#MITE /' *plus_Mollusca_nr.fasta > "$varSpecie"_plus_Mollusca_nr_renamed.fasta; 
  sed -i 's/SUB.*\]/RC\/Helitron/g' "$varSpecie"_plus_Mollusca_nr_renamed.fasta;
  cd ../../; 
done;

#----------------------------------Blast filter-------------------------------------------------------------------#

for i in $(ls | grep "TransposableElements"); do 
  varSpecie=$( echo "$i" | cut -d"_" -f1); 
  cd "$i"/Libraries; 
 #Print list of consensus with at least 5 good back-blast hits against the source genome
  python ../../Blast_Filter.py --lib "$varSpecie"_plus_Mollusca_nr_renamed.fasta --out "$varSpecie" --num_threads 30 --blast_query_cov 70 --blast_identity 70 --genome /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/ALL_Genomes/"$varSpecie"_genome.fasta --min_Blast_Hits 5; 
  cd ../../;
done;

#Recover good consensus sequences
for i in $( ls | grep "TransposableElements"); do
  varSpecie=$(echo "$i" | cut -d"_" -f1);
  cd "$i"/Libraries;
  awk -F"\t" '{print$1}' "$varSpecie"*BlastFiltered.txt > Consensus_ToKeep.txt;
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "$varSpecie"_plus_Mollusca_nr_renamed.fasta > "$varSpecie"_plus_Mollusca_nr_renamed.oneliner
  grep -A1 -Ff Consensus_ToKeep.txt "$varSpecie"_plus_Mollusca_nr_renamed.oneliner > "$varSpecie"_FinalFiltered.fasta;
  sed -i 's/--//g' "$varSpecie"_FinalFiltered.fasta;
  sed -i '/^$/d' "$varSpecie"_FinalFiltered.fasta;
  cd ../../
done;

#------------------------------------RepeatMasker-----------------------------------------------------------------#

for i in $(ls | grep "_TransposableElements"); do 
  varSpecie=$( echo "$i" | cut -d"_" -f1); 
  mkdir "$i"/RepeatMasker; 
  cd "$i"/RepeatMasker; 
  RepeatMasker -lib ../Libraries/*FinalFiltered.fasta -pa 20 -gff -a -no_is -nolow -norna -s /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/ALL_Genomes/"$varSpecie"*.fasta;
cd ../../; 
done;
