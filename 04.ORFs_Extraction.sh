#!bin/bash
#Improve TE annotation with RepeatCraft; extend to 1000bp at both ends all corrected insertions, extract them and screen for ORF longer than 900 nucleotides.
#Finally run RPSBlast against Pfam database.

export LD_LIBRARY_PATH=/home/jacopomartelossi/.conda/envs/ORFinder/lib

#------------------------------------------------RepeatCraft for better TEs annotation-----------------------------------------#

for i in $( ls | grep "TransposableElements"); do 
  varSpecie=$( echo "$i" | cut -d"_" -f1); cd "$i"/RepeatMasker;
  /media/storage/jacopomartelossi/Software/repeatcraftp/repeatcraft.py -r *gff -u *out -o "$varSpecie"_RepeatCraft -c /media/storage/jacopomartelossi/Software/repeatcraftp/repeatcraft.cfg --mode loose; 
  cd ../../; 
done;

#---------------------------------------------------Extraction of each TE copy-------------------------------------------------#

for i in $( ls | grep "TransposableElements"); do 
  varSpecie=$( echo "$i" | cut -d"_" -f1); 
  echo "$i"; cd "$i"/RepeatMasker; 
  awk -v OFS="\t" -F"\t" '{print$1,$4,$5,$7,$6,$9}' *rmerge.gff | sed 's/Tstart.*ID=//g' | sed 's/;.*$//g' > RepMasker.bed; 
  cat RepMasker.bed | awk -v OFS="\t" -F"\t" '{print$1,$2,$3,$6,$5,$4}' > RepMasker.modified.bed; 
  bedtools slop -i RepMasker.modified.bed -g /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/ALL_Genomes/"$varSpecie"_genome.fasta.fai -b 1000 | bedtools merge -s -c 4,6 -o distinct > RepMasker.modified_1k.bed  
  bedtools getfasta -s -fi /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/ALL_Genomes/"$varSpecie"_genome.fasta -bed RepMasker.modified_1k.bed > RepMasker_ExtractedSequences.fasta; 
  cd ../../; 
done;

#--------------------------------Screen all insertions for ORFs longer than 900 nucleotides-----------------------------------#

for i in $( ls | grep "TransposableElements"); do 
  varSpecie=$( echo "$i" | cut -d"_" -f1);
  cd "$i"/RepeatMasker;
  sed -i 's/:/_/g' RepMasker_ExtractedSequences.fasta; 
  mkdir Autonomous_Search; 
  /media/storage/jacopomartelossi/Software/ORFfinder -ml 900 -strand plus -n true -outfmt 1 -in RepMasker_ExtractedSequences.fasta -out "$varSpecie"_ExtractedSequences_ORF.fasta
  /media/storage/jacopomartelossi/Software/ORFfinder -ml 900 -strand plus -n true -outfmt 0 -in RepMasker_ExtractedSequences.fasta -out "$varSpecie"_ExtractedSequences_ORF.pep
  cd Autonomous_Search
  
#------------------------------------------------------RPSBlast against Pfam---------------------------------------------------#
  rpsblast -query "$varSpecie"_ExtractedSequences_ORF.pep -evalue 1e-5 -db /media/storage/jacopomartelossi/dbs/Pfam_rpsblast/Pfam -outfmt 6 -out "$varSpecie"_ORFSearch.rpsblast -num_threads 30  
  cd ../../../; 
done;
