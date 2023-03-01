#!bin/bash
#Script to identify DDE/D related ORFs based on a collection of HMM profiles builded up from  https://www.pnas.org/doi/10.1073/pnas.1104208108 dataset.

#---------------------------------------Identify DDE/D related signatures in previously identified ORFs------------------------------#
for i in *.pep; do
  hmmscan -E 1e-05 --cpu 8 --domtblout "${i/.pep/ClassII.domtblout}" --tblout "${i/.pep/ClassII.tblout}" ../DDE_D.hmm "$i"
  #Retrive best hits
  awk '!x[$4]++' "${i/.pep/ClassII.domtblout}" > "${i/.pep/ClassII.domtblout.BestHits}"
done;

#------------------------------------------------Retrive DDE/D related ORFs sequences------------------------------------------------#
for i in *domtblout.BestHits; do
  varSpecie=$( echo "$i" | cut -d"_" -f1);
  while read line; do 
    #Create a variable with name of the best-hitted superfamily
    varName=$( echo "$line" | awk '{print$1}'); 
    echo "$varName"; 
    #Create variable with name of the sequence
    varSeq=$(echo "$line" | awk '{print$4}'); 
    echo "$varSeq";
    #Retrive each hitted sequence and add specie name and Superfamily classification to the header
    grep  -w -A1 "$varSeq" ../"$varSpecie"_ExtractedSequences_ORF.pep | sed 's/>.*/>'"$varSpecie"'_'"$varName"'/'  >> "$varSpecie"_DDE_D.pep  
  done< "$i"
  #Add a progressive number at the end of each header
  awk '/^>/{$0=$0"_"(++i)}1' "$varSpecie"_DDE_D.pep  > "$varSpecie"_DDE_D.Renamed.pep
done

#--------------------------------------Confirm DDE/D hits against the RepeatPeps lib------------------------------------------------------#
for i in *.pep; do 
  blastp -evalue 1e-05 -query "$i" -db /media/storage/jacopomartelossi/dbs/RepeatPeps.lib -outfmt 6 -out "$i".Blastp -num_threads 10 -max_target_seqs 1 -max_hsps 1; 
 done
