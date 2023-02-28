#!bin/bash
#Identify LINEs derived RVT-containiing ORFs (i.e. excuding LTRs)

mkdir ALL_LINES

for i in $(ls | grep "_TransposableElements"); do
  varSpecie=$( echo "$i" | cut -d"_" -f1);
  cd "$i"/RepeatMasker/Autonomous_Search;
  mkdir LINEs

#-----------------------------Reformat rpsblast output and resolve overlapping hits-------------------------------------------------#
  
  awk -v OFS=" " -F"\t" '{print$1,$2";"$11,$12,$7"-"$8}' "$varSpecie"_ORFSearch.rpsblast  > "$varSpecie"_ORFSearch_reformatted.rpsblast
  /media/storage/jacopomartelossi/Software/cath-tools/cath-resolve-hits.ubuntu-20.04 --input-format raw_with_scores "$varSpecie"_ORFSearch_reformatted.rpsblast > "$varSpecie"_ORFSearch_ResolvedOverlap.rpsblast

#--------------------------------Extract LINEs that match the conditions (i.e. exhibit and RVT domain) ------------------------------#
 
 #RVT_Simplified_CDD.txt is a tab-delimited list of CDD entries related to RVT domains. We also remove 
 awk -F"\t" '{print$1}' /media/storage/jacopomartelossi/Evolution_Transposable.Elements/data/RVT_Simplified_CDD.txt | grep -Ff - "$varSpecie"_ORFSearch_ResolvedOverlap.rpsblast | cut -d" " -f1 > LINEs/"$varSpecie"_LINES.txt

#-------------------------------Extract all  LINEs autonomous copies and confirm LINE RVT domain through hmmscan---------------------#
  #RVT_MolluscaP.hmm is a simple database composed by two HMM profiles builded up from LTRs and LINEs-derived RVT
  awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$varSpecie"_ExtractedSequences_ORF.pep | grep -Ff LINEs/"$varSpecie"_LINES.txt - | tr "\t" "\n" > LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.fasta
  hmmscan -E 1e-05 --cpu 8  --domtblout LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.domtblout --tblout LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.tblout ../../../../data/RVT_MolluscaP.hmm LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.fasta
  #Print bet-hitted domain for each query sequences
  awk '!x[$4]++' LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.domtblout > LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.domtblout_BestHits
  cat LINEs/"$varSpecie"_Putative_LINES_RVT.ORF.domtblout_BestHits | grep "LINE"  | grep -v "#" | awk '{print$4}' > LINEs/"$varSpecie"_LINEs_RVT.ORF_Confirmed.txt
  sed -i 's/lcl|.* /lcl|/g' "$varSpecie"_ExtractedSequences_ORF.fasta
  #Recover both nucleotides and protein sequences of each confirmed LINE-derived ORF.
  awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$varSpecie"_ExtractedSequences_ORF.fasta | grep -Ff LINEs/"$varSpecie"_LINEs_RVT.ORF_Confirmed.txt - | tr "\t" "\n" > ../../../ALL_LINES/"$varSpecie"_LINEs_RVT.ORF.fasta;
  awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$varSpecie"_ExtractedSequences_ORF.pep | grep -Ff LINEs/"$varSpecie"_LINEs_RVT.ORF_Confirmed.txt - | tr "\t" "\n" > ../../../ALL_LINES/"$varSpecie"_LINEs_RVT.ORF.pep;
  cd ../../../
done;
