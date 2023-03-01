RVT_FASTA=$1 #Fasta files with LINEs-derived RVT. It must end with .fa

#--------------------------------------------Align LINEs-derived RVT---------------------------------------------------#
ginsi --thread 20 "$RVT_FASTA" > "${RVT_FASTA/.fa/.ginsi}"
trimal -resoverlap 0.75 -seqoverlap 80 -in "${RVT_FASTA/.fa/.ginsi}" -out "${RVT_FASTA/.fa/_Nosp.ginsi}"

#------------------------------------------NJ and ML analyses----------------------------------------------------------#
clearcut --alignment -P --in "${RVT_FASTA/.fa/_Nosp.ginsi}" --shuffle --neighbor --out=LINE_NJ.Tree.clearcut
iqtree2 -s "${RVT_FASTA/.fa/_Nosp.ginsi}" -nt AUTO -wbt -bb 1000 -m MFP --prefix Unconstrained_ML --runs 5

#--------------------------------------------Constrained ML analyses---------------------------------------------------#
#First constrained ML tree search using the topology obtained with Clearcut (basically we are only claculating branch lengths and the likelihood of the tree)
iqtree2 -s "${RVT_FASTA/.fa/_Nosp.ginsi}" -m MFP -te LINE_NJ.Tree.clearcut -nt AUTO --prefix NJ-Full_Costrain
#Second constrained ML tree search constraining only monophyly of LINE superfamilies, as recovered by the NJ tree search
iqtree2 -s "${RVT_FASTA/.fa/_Nosp.ginsi}" -m MFP -g NJ-SupFam.newick -nt 8 -bb 1000 --runs 5 --prefix NJ-SupFam_Constrain

#-----------------------------------------Tree topology test-----------------------------------------------------------#
mkdir Topology_Test
cat Unconstrained_ML.runtrees NJ-Full_Costrain.treefile NJ-SupFam_Constrain.runtrees > Topology_Test/ALL_Trees.newick
cd Topology_Test
#Here we will use the Q.yeast+F+F+R10 substitution model because it was the best-fit in the uncostrained ML tree searches
iqtree2 -s ../"${RVT_FASTA/.fa/_Nosp.ginsi}" -z ALL_Trees.newick -m Q.yeast+F+F+R10 -n 0 -zb 1000 -au -nt AUTO -zw
