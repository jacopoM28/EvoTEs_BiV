#!bin/bash
#Script to run all de-novo TE mining softwares (RepeatModeler2, MITE Tracker and HelitronScanner)

Path_Genome="$1"
Genome_database_name="$2"
N_of_cores="$3"
Path_MITEs_tracker="$4" #Path MITE_tracker python installation
Path_HelitronScanner="$5" #Path Helitronscanner installation

eval "$(conda shell.bash hook)"
varGenome=$( echo "$1" | rev | cut -d "/" -f1  | rev )

mkdir RepeatModeler
mkdir MITEs
mkdir Helitrons

###################################RepeatModeler##################################################
pathMaster=$(pwd)
conda activate Repeats

cd ./RepeatModeler
ln -s "$1"

#Path to NINJA installation. Necessary for LTR pipeline of RepeatModeler
NINJA_DIR=/home/jacopomartelossi/.conda/envs/Repeats/bin
BuildDatabase -name "$2" -engine ncbi "$varGenome"
RepeatModeler -pa "$3" -engine ncbi -database "$2" -LTRStruct -debug 2>&1 | tee LTR_struct.log

cd ..
###################################MITEs Tracker###################################################
conda deactivate
conda activate MITE_tracker

cd ./MITEs
varMITEs=$(pwd)
cd "$4"
ln -s "$1"

python3 MITETracker.py -g "$varGenome"  -j "$2"_MITEs -w "$3"
mv ./Results/"$2"_MITEs "$varMITEs"

cd "$pathMaster"

####################################Helitron Scanner##############################################
cd Helitrons
ln -s "$1"

java -jar "$5"/HelitronScanner.jar scanHead -lf "$5"/TrainingSet/head.lcvs -g "$varGenome" -bs 1000000 -o "$2".head
java -jar "$5"/HelitronScanner.jar scanTail -lf "$5"/TrainingSet/tail.lcvs -g "$varGenome" -bs 1000000 -o "$2".tail
java -jar "$5"/HelitronScanner.jar pairends -hs "$2".head -ts "$2".tail -o "$2".pairends -ht 10 -tt 10
java -jar "$5"/HelitronScanner.jar draw -g "$varGenome" -p "$2".pairends -o "$2".fasta --pure
