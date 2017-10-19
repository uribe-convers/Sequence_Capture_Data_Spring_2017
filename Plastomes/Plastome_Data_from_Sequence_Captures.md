# Getting Plastome Data from Sequence Capture

### Simon Uribe-Convers - June 05th, 2017
---

I used Michael McKain's pipeline [Fast-Plast](https://github.com/mrmckain/Fast-Plast) to try to assemble complete plastomes from the off-target reads from our sequence capture experiment. These are the commands I used on different machines.


## Pancho
```
while read name;
do perl /Applications/Fast-Plast-master/fast-plast.pl -1 /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/$name"*PE1.fastq" -2 /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/$name"*PE2.fastq" --single /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/$name"*SE.fastq" --name $name --user_bowtie /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/CBS_Bowtie_Index/CBS --adapters TruSeq --coverage_analysis
done < Name_List.txt

```
## Pancho Combined Data on Molaua

```
while read name
do perl /Applications/Fast-Plast-master/fast-plast.pl -1 /Volumes/Neobartsia/Projects/Burmeistera/Genomics/Data/Capture_Probe_Raw_Data_Spring_2017/Combined/Combined_Cleaned/$name"*PE1.fastq" -2 /Volumes/Neobartsia/Projects/Burmeistera/Genomics/Data/Capture_Probe_Raw_Data_Spring_2017/Combined/Combined_Cleaned/$name"*PE2.fastq" --single /Volumes/Neobartsia/Projects/Burmeistera/Genomics/Data/Capture_Probe_Raw_Data_Spring_2017/Combined/Combined_Cleaned/$name"*SE.fastq" --name $name --user_bowtie /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/CBS_Bowtie_Index/CBS --adapters TruSeq --coverage_analysis --skip trim --clean light
done < Name_List.txt
```

```
perl /Applications/Fast-Plast-master/fast-plast.pl -1 /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/37_C_erythraeus_R1.fastq.gz_Cleaned_PE1.fastq -2 /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/37_C_erythraeus_R1.fastq.gz_Cleaned_PE2.fastq --single /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/Data/37_C_erythraeus_R1.fastq.gz_Cleaned_SE.fastq --name 37_C_erythraeus --bowtie_index Asterales --adapters TruSeq --coverage_analysis
```

## Cluster
```
while read name;
do perl ~/bin/Fast-Plast/fast-plast.pl -1 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE1.fastq" -2 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE2.fastq" --single /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*SE.fastq" --name $name --bowtie_index Asterales --adapters TruSeq --coverage_analysis
done < Name_List.txt
```

## Using my own CBS Index

```
while read name;
do gunzip ../*gz
perl ~/bin/Fast-Plast/fast-plast.pl -1 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE1.fastq" -2 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE2.fastq" --single /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*SE.fastq" --name $name --user_bowtie /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/Plastomes/CBS_Bowtie_Index/CBS --adapters TruSeq --coverage_analysis
done < Name_List_B_ceratocarpa_100bp.txt

```

```
while read name;
do perl ~/bin/Fast-Plast/fast-plast.pl -1 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE1.fastq" -2 /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*PE2.fastq" --single /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/$name"*SE.fastq" --name $name --user_bowtie /mnt/lfs2/dtank/HybPiper_Simon/Combined_Data/Plastomes/CBS_Bowtie_Index/CBS --adapters TruSeq --coverage_analysis
done < Name_List_37.txt
```

## Move log files to each species directory
```
for i in *log; do echo ${i%%_Comb*}_Combined; done
for i in *log; do mv $i ${i%%_Comb*}_Combined; done
```

## Move final Scaffolds or contigs to one directory

```
for i in */Final_Assembly/*fasta; do cp $i /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/-Final_Contigs_and_Scoffolds;done
```
---
---
# Mapping contigs or scaffolds obtained from a Fast-Plast assembly back to a reference plastome



After successfully assembling large contigs of the plastome from my sequence capture data using Michael McKain's awesome pipeline [Fast-Plast](https://github.com/mrmckain/Fast-Plast), I needed a way to map the results back to a reference. This was needed because none of my samples produced a complete plastome; this is not surprising because I'm using the off-target reads to accomplish this. 


I first tried to do it on Sequencher, and although the program works, it requires a lot of manual work to do it. 

After some hard core Googling, I found a program called [CONTIGuator](http://contiguator.sourceforge.net) that was design to do exactly this task on microbial genomes. The program is no longer being maintained and they recommend using their other program [Medusa](http://combo.dbe.unifi.it/medusa). Both programs can be run online but they also provide the source code to run locally, yay! 

Easy steps:

1. Make a file with the names of every sample's final contig/scaffold file  
	`ls > Name_list.txt`
2. Run the code below

#### Running Medusa

```
mkdir Final_PseudoContigs
while read i
do
java -jar medusa.jar -f Reference_CBS_Plastomes/ -i /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/-Final_Contigs_and_Scoffolds/$i -o ${i%%.fa*}"Final_PseudoContigs.fasta" -v 1>> Medusa.log
done < Name_list.txt

mv *Final_PseudoContigs.fasta Final_PseudoContigs

```

#### Running CONTIGuator

```
while read i;
do
mkdir ${i%%.fa*}
cd ${i%%.fa*}
python /Applications/CONTIGuator_v2.7/CONTIGuator.py -r /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/-Final_Contigs_Mapped_Against_Reference/Reference_CBS_Plastomes/Burmeistera_borjensis_FINALgenome.fasta -c /Users/SAI/Documents/-Projects/--Burmeistera/Genomics/Sequence_Capture_Data_Spring_2017/Plastomes/-Final_Contigs_and_Scoffolds/$i
mv PseudoContig.fsa $i"_PseudoContigs.fasta"
cd ../
done < Name_list.txt

mkdir Final_PseudoContigs
for i in */*PseudoContigs.fasta; do cp $i Final_PseudoContigs/; done

```

I aligned all the sequences, but then I wanted to only work with species that have a certain amount of nuclear genes, what we are calling the "Good Data". To separate these species, I did the following grep command on the alignment. 


```
grep -A 1 -e \>Burmeistera_almedae_2016_001_Combined  -e \>Burmeistera_anderssonii_2016_113_Combined  -e \>Burmeistera_aspera_2016_089_Combined  -e \>Burmeistera_auriculata_2016_199_Combined  -e \>Burmeistera_borjensis_2016_219_Combined  -e \>Burmeistera_brachyandra_2016_114_Combined  -e \>Burmeistera_brighamioides_2016_121_Combined  -e \>Burmeistera_bullatifolia_2016_129_Combined  -e \>Burmeistera_ceratocarpa_2016_233_Combined  -e \>Burmeistera_cf_aeribacca_2016_087_Combined  -e \>Burmeistera_chiriquiensis_2016_161_Combined  -e \>Burmeistera_crassifolia_2016_200_Combined  -e \>Burmeistera_crispiloba_2016_201_Combined  -e \>Burmeistera_cyclostigmata_2016_104_Combined  -e \>Burmeistera_cylindrocarpa_2016_235_Combined  -e \>Burmeistera_darienensis_2016_149_Combined  -e \>Burmeistera_domingensis_2016_192_Combined  -e \>Burmeistera_draconis_2016_133_Combined  -e \>Burmeistera_dukei_2016_153_Combined  -e \>Burmeistera_fuscoapicata_2016_128_Combined  -e \>Burmeistera_glabrata_2016_227_Combined  -e \>Burmeistera_holm-nielsenii_2016_210_Combined  -e \>Burmeistera_huacamayensis_2016_115_Combined  -e \>Burmeistera_litensis_2016_214_Combined  -e \>Burmeistera_loejtnantii_2016_006_Combined  -e \>Burmeistera_lutosa_2016_209_Combined  -e \>Burmeistera_mcvaughii_2016_155_Combined  -e \>Burmeistera_multiflora_2016_196_Combined  -e \>Burmeistera_obtusifolia_2016_108_Combined  -e \>Burmeistera_oyacachensis_2016_230_Combined  -e \>Burmeistera_panamensis_2016_151_Combined  -e \>Burmeistera_parviflora_LL_20_Combined  -e \>Burmeistera_pirrensis_2016_148_Combined  -e \>Burmeistera_quercifolia_2016_009_Combined  -e \>Burmeistera_ramosa_2016_007_Combined  -e \>Burmeistera_refracta_2016_195_Combined  -e \>Burmeistera_resupinata_2016_240_Combined  -e \>Burmeistera_resupinata_var_heilbornii_2016_116_Combined  -e \>Burmeistera_rubrosepala_2016_239_Combined  -e \>Burmeistera_smaragdi_2016_236_Combined  -e \>Burmeistera_smooth_2016_135_Combined  -e \>Burmeistera_sodiroana_2016_197_Combined  -e \>Burmeistera_succulenta_2016_143_Combined  -e \>Burmeistera_succulenta_2016_188_Combined  -e \>Burmeistera_tenuiflora_2016_158_Combined  -e \>Burmeistera_truncata_2016_144_Combined  -e \>Burmeistera_utleyi_2016_156_Combined  -e \>Burmeistera_vulgaris_2016_100_Combined  -e \>Burmeistera_xerampelina_2016_086_Combined  -e \>Burmeistera_zurquiensis_2016_098_Combined  -e \>LL6_C_smithii_Combined  -e \>LL49_C_incanus_Combined  -e \>LL61_C_asclepideus_Combined  -e \>LL69_C_nigricans_Combined  -e \>LL83_C_brittonianus_Combined  -e \>LL86_C_mandonis_Combined  -e \>LL88_S_ayersiae_Combined  -e \>LL159_S_jelskii_Combined  -e \>LL334_S_aureus_Combined  -e \>LL363_S_krauseanus_Combined All_Plastomes_Burmeistera.fasta_non-interleaved.aln > Species_with_Good_Data/All_Plastomes_Burmeistera_Good_Data.aln 

```

## Phylogenetics

I did my usual phylogenetic approach to these data, align, clean, phylo.

Copy the code below and run the script with `./Align_Clean_Phylo.sh`

#### Code for phylogenetics
```
#!/bin/bash

# Load modules if running on a cluster
module load mafft
module load phyutility
module load fasttree

# Delete any empty files
find . -size 0 -delete

# Create directroies

mkdir Alignments Fasta Phylo

# Align with Mafft, clean with Phyutility, and analyze with FastTree

for i in *.fasta
do
# Align
echo ""
echo "~~~Aligning with Mafft~~~"
echo ""
time mafft --auto --thread 6 --preservecase $i > $i.aln

# Cleaninig Alignment
echo ""
echo "~~~Cleaninig alignment with Phyutility at 0.5~~~"
echo ""
time phyutility -clean 0.5 -in $i".aln" -out $i"_cleaned_05.aln"

# Phylogenetics
echo ""
echo "~~~Building phylogeny with FastTre~~~"
echo ""
time FastTree -nt -gtr < $i"_cleaned_05.aln" > $i".tre"

# House keeping
mv $i".aln" $i"_cleaned_05.aln" Alignments
mv $i".tre" Phylo
mv $i Fasta
done

```
