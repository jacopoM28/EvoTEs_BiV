#!bin/bash
##Automatic filtering of species-specific repeat libraries

#----------------------------------------------Merge all libraries--------------------------------------------------------#

for i in $(ls | grep "_TransposableElements"); do 
 varSpecie=$(echo "$i" | cut -d"_" -f1); 
 mkdir "$i"/Libraries; 
 cat "$i"/MITEs/families_nr.fasta "$i"/RepeatModeler/*families.fa "$i"/Helitrons/*hel.fa > "$i"/Libraries/"$varSpecie"_AllLibraries.fa; 
done;

#-----------------------------Blastx against a clean set of proteins (i.e. without proteins related to TEs)----------------------#

for i in $(ls | grep "_TransposableElements"); do 

varSpecie=$(echo "$i" | cut -d"_" -f1); 
mkdir "$i"/Libraries/RemoveGeneFragments; 
cp "$i"/Libraries/*_AllLibraries.fa "$i"/Libraries/RemoveGeneFragments; 
blastx -query "$i"/Libraries/RemoveGeneFragments/*AllLibraries.fa -db /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/Cleaned_proteoms.fa -evalue 1e-10 -num_descriptions 10 -out "$i"/Libraries/RemoveGeneFragments/"$varSpecie"_Blastx.out -num_threads 20; 
#NB: The database was previously cleaned from TE-related proteins using Blastp

done;

#--------------------------ProtExcluder to remove genes and gene fragments--------------------------------------------------#

for i in $(ls | grep "_TransposableElements"); do 
 cd "$i"/Libraries/RemoveGeneFragments/; 
 perl /media/storage/jacopomartelossi/Software/ProtExcluder1.1/ProtExcluder.pl *Blastx.out *AllLibraries.fa; 
 cd ../../../; 
done;

#------------------------------------Remove tandem repeats-----------------------------------------------------------------#

for i in $(ls | grep "_TransposableElements"); do 
 cd "$i"/Libraries/; 
 mkdir Clean_Tandem; 
 varSpecie=$( echo "$i" | cut -d"_" -f1); 
 perl /media/storage/jacopomartelossi/Software/cleanup_tandem.pl -f RemoveGeneFragments/*fanoProtFinal -minlen 50 -misschar n -nr 0.5 > Clean_Tandem/"$varSpecie"_NoProt_NoTandem.fasta; 
 cd ../../; 
done;

#-------------------------------------------Merge with Mollusca library and remove redundancy following 80-80 rule-------#

for i in $(ls | grep "_TransposableElements"); do 
 varSpecie=$( echo "$i" | cut -d"_" -f1); 
 cd "$i"/Libraries/; 
 cat Clean_Tandem/*NoTandem.fasta /media/storage/jacopomartelossi/dbs/Mollusca_RepBaseDfam_14-09-21.fasta > "$varSpecie"_plus_Mollusca.fasta; 
 cd-hit-est -T 0 -i *Mollusca.fasta -M 1000000 -o "$varSpecie"_plus_Mollusca_nr.fasta -c 0.8 -n 5 -aS 0.8 -g 1 -G 0
 cd ../../; 
done;
