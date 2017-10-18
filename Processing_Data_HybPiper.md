# Installing and running HybPiper to process sequence capture data  
_Simon Uribe-Convers – April 05th, 2017 – www.simonuribe.com_
---

Sequence capture experiments have become very common among evolutionary biologists because of their high throughput, affordable prices, and the large amounts of data they produce. Because of their popularity, there has been a high demand for approaches to process the data in an efficient way. One of these methods is **[HybPiper](https://github.com/mossmatters/HybPiper)**, a [recently published](http://www.bioone.org/doi/full/10.3732/apps.1600016) pipeline designed and developed by [Matthew Johnson](http://mossmatters.net) and [Norman Wickett](https://www.chicagobotanic.org/research/staff/wickett) from the Botanic Chicago Garden.

The pipeline is very thorough and comprehensive, and Matt and Norm thought of everything your would need to process your data. Starting with demultiplexed clean reads (no barcodes or adapters), the pipeline 1) sorts the reads by mapping them to target sequences, using BLASTx (protein targets) or BWA (nucleotide targets), 2) assembles contigs for each gene separately, and 3) aligns contigs to target sequences and extracts exon (coding) sequence. It also produces very informative statistics for every step, identifies non-coding flanking regions (introns), and it even warns you of putative paralogous sequences, and gives you methods to help distinguish ancient from recent paralogs! What else can you ask for?

I recently used the pipeline with the latest data I received from 152 samples of *Burmesitera, Centropogon*, and *Siphocampylus*, and I could not be happier or more impressed with HybPiper. Plus the customer service is great, thanks Matt! :)

Here is my experience with HybPiper ran in a parallel computing cluster. You might need to modify some commands and **definitively** paths to files and/or scripts!


## Cleaning your data
Before running HybPiper, make sure your data are demultiplexed and clean of barcodes and adaptors. You can check if the data are clean with [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

`fastqc file`  

In case the data need to be cleaned, do this with SeqyClean or Trimmomatic.  
*Actually,* it's a good idea to clean the data to make sure that they don't contain adapters and to discard reads with low qualities.  

**[In SeqyClean](https://github.com/ibest/seqyclean)**

```{bash}
seqyclean -1 forward_read -2 reverse_read -o basename -qual
```

###For many files in a directory, code from [here](http://stackoverflow.com/questions/29306245/using-trimmomatic-on-multiple-illumina-paired-end-read-files), and Information on how this works [here:](https://debian-administration.org/article/150/Easily_renaming_multiple_files)


```{bash}
#Directory for clean reads
mkdir Cleaned_Data_SeqyClean_May_2017

#A for loop to run the program on many reads changing the variable mid loop to include the reverse read.
for f1 in *_R1.fastq.gz
do
# The line below changes R1 for R2. The way it is doing it is by stripping part of the name
# stored in the variable (i.e., R1.fastq.gz) and replacing it by something else (i.e., R2.fastq.gz)
f2=${f1%%R1.fastq.gz}"R2.fastq.gz"
# This can also be done with
# f2=${f1/R1.fastq.gz/}"R2.fastq.gz"

seqyclean -1 $f1 -2 $f2 -o ./Cleaned_Data_SeqyClean_May_2017/$f1"_Cleaned" -qual
done


#Dealing with the additional files SeqyClean produces

mkdir Reports
cd Cleaned_Data_SeqyClean_May_2017
mv *.txt *.tsv ../Reports
cd ../Reports
#appends all summary files into one
cat *.tsv >> All_Stats_temp.tsv
#removes duplicated headers
awk '!seen[$0]++' All_Stats_temp.tsv > All_Stats.tsv
rm *temp*
sed -i '/^$/d' All_Stats.tsv

#appends all report files into one
cat *.txt >> All_Reports.txt

```

**[In Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)**

```{bash}
java -jar path-to-program/trimmomatic-0.36.jar PE -phred33 Forward_reads Reverse_reads output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:path-to-adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:10:20 MINLEN:40 [other options]
```

### For many files in a directory


```{bash}
for f1 in RAPiD*R1_001.fastq.gz
do
f2=${f1%%R1_001.fastq.gz}"R2_001.fastq.gz"
java -jar /Applications/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $f1 $f2 $f1"_output"_forward_paired.fq.gz $f1"_output"_forward_unpaired.fq.gz $f1"_output"_reverse_paired.fq.gz $f1"_output"_reverse_unpaired.fq.gz ILLUMINACLIP:/Applications/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:10:20 MINLEN:40
done
```

I went with **SeqyClean** because it produces only three files, PE1, PE2, and SE, which are accepted in HybPiper

## Renaming the reads
My raw reads had some very long and uninformative names. I chose to rename the files to match the sample name (species) and some text informing me that they have been cleaned and if it's the forward (PE1), reverse (PE2) or single (SE) read.

We made a "sample sheet" that we sent to RapidGenomics with the names of each sample and they appended the corresponding filename for each sample. The spreadsheet looked like this:

Filename|Species
---|---
filename1|species 1
filename2|species 2

Because I know that `SeqyClean` will result in three files (PE1, PE2, and SE), I needed to triplicate each row. This was easily done in TextWrangler, searching for each row (`^(.+)`) and replacing it with itself three times (`\1\r\1\r\1`). Once I had this file with three identical lines per row, I added a PE1, PE2, and SE to one of each lines per species and the file extension.

Filename|Species
---|---
filename1|species1_PE1.fastq
filename1|species1_PE2.fastq
filename1|species1_SE.fastq
filename2|species2_PE1.fastq
filename2|species2_PE2.fastq
filename2|species2_SE.fastq

Finally, I added the command to move (`mv`) to the beginning of each line and ran it on the terminal:

```{bash}
mv filename1 species1_PE1.fastq
mv filename1 species1_PE2.fastq
mv filename1 species1_SE.fastq
mv filename2 species2_PE1.fastq
mv filename2 species2_PE2.fastq
mv filename2 species2_SE.fastq
```
Most of this was done with search and replace in TextWrangler.

## Running HybPiper

HybPiper requires a few programs to run. The easiest way to get everything installed, at least for Macs, is with [Homebrew](https://brew.sh).

Install HomeBrew if it's not on your computer yet, it's the easiest way to deal with convoluted installations.

On the terminal:

```{bash}
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Install programs and dependencies

```{bash}
brew tap homebrew/science
brew install exonerate
brew install bwa
brew install samtools
brew install spades
brew install blast
brew install parallel
brew install gcc
```

The pipeline runs on python 2.7 or later (already available in Macs) but it also requires [Biopython](http://biopython.org//wiki/Biopython). The **easiest** way to get Biopython installed on your machine, both PCs and Macs, is with [Anaconda](https://www.continuum.io/downloads).


Once you have all the programs installed via Homebrew and Biopython installed via Anaconda, download and install HybPiper either by downloading the zip file from [GitHub](https://github.com/mossmatters/HybPiper) or from the command line with

```
git clone https://github.com/mossmatters/HybPiper.git
```
If you have any problems, there are very detailed instruction on how to install HybPiper and all its dependencies [here](https://github.com/mossmatters/HybPiper/wiki/Installation)

Once everything is downloaded and installed, **check that everything is working properly!**

```{bash}
python path_to_HybPiper-master/reads_first.py --check
```


## Additional suggested installs

There are two programs that I cannot live without while working with DNA sequences. Both of them deal with format conversion (e.g., fasta to nexus) or concatenation very smoothly. These programs are:

Download and install [Phyutility](https://code.google.com/archive/p/phyutility/downloads) to do file format conversion.

Download and install [NCLconverter](http://ncl.sourceforge.net) for further file format conversion.


## Get your data ready for HybPiper
### The Target File
You need a "Target File" containing the genes used to designed the capture probes. If there are multiple exons per gene, the exons need to be concatenated.  
Our genes come from three different analyses and have different names. The first thing I did was to standardize the names within each probe analysis. For example, a gene containing multiple exons called `Hydra\_AT1G07970\_exon1, Hydra\_AT1G07970\_exon2, Hydra\_AT1G07970\_exon3` was renamed to just `Hydra\_AT1G07970`. This way I was able to concatenate the multiple exons easier. Because my sequences for the genes were all in a single fasta file, I renamed them in [TextWrangler](http://www.barebones.com/products/textwrangler/) using search and replace. Also, make sure that there are no empty lines between the multiples genes.
For my example:

```
Search: _exon.+
Replace: "nothing"
Search: ^\s+
Replace: "nothing"
```

With those changes, I was able to split the large files with all the fasta sequences into multiple files containing a single fasta sequence. This was necessary to concatenate the multiple exons per genes into a single continuous sequence. I did this using `split`, using the option `-l 2` to include two lines per files and `-a 15` to increase the suffix length for the new files created, i.e., this adds a unique letter combination to each file after its name.


```{bash}
mkdir Single_Fasta
split -l 2 -a 15 Final_Probes_Muchhala_et_al_Named_Sets.fasta ./Single_Fasta/gene
cd Single_Fasta
```

The next step is to concatenate the different exons for each gene. I'm sure there are more elegant ways to do this, but I used [**Phyutility**](http://blackrim.org/programs/phyutility/) to treat each new fasta file as a gene and to make a single large Nexus file. **However**, After over four hours, this process hadn't finished so I canceled it.

```{bash}
**This didn't finish:**
phyutility -concat -in gene* -out All_Genes_Concatenated.nex
```

I proceeded to break the concatenation step into smaller number of files. In my case, there are 2872 files, so to just work with one fourth of them (718 files) at a time instead of all of them I did:

```{bash}
#How many batches (directories) do you want to work with?
DIRECTORIES=4

#Counts number of files in a directory and stores it in a variable
NUMFILES=$(ls | wc -l)
#Divides the variable by a number specified in a variable called `DIRECTORIES`, which is the number of files I wanted to work with. It stores this number and its decimals (files/$DIRECTORIES) into a variable
NUMFILESDIR=$(echo "$NUMFILES/$DIRECTORIES"|bc -l)
#Divides the variable by four, which is the number of files I wanted to work with. It stores this number with NO decimals into a variable
NUMFILESDIR_NO_FLOAT=$(($NUMFILES / $DIRECTORIES))
#Informs how many files need to go to each directory
echo "<< The number of files per directory is ***$NUMFILESDIR*** >>"

eval mkdir {01..$DIRECTORIES} Used_Gene_Files

# In a case working with only 4 batches!
#Moves the right number of files into each directory (see below for a caveat)
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./01; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./02; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./03; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./04; done

```

***IMPORTANT***  
**_If_** the number of files per directory is divisible by 4, giving you an entire number (see the prompt on the terminal "The number of files per directory is X"), you don't have to do anything else. **_If NOT_**, only a portion of the files has been moved to the directories you just created because you can only move files in entire numbers! You'll need to manually move the remaining files to the directory **number 1** before continuing with the code below, which should take 10 to 15 minutes to run.

```{bash}
cd 1; phyutility -concat -in gene* -out temp1.nex; mv gene* ../Used_Gene_Files; mv temp* ..; cd ../2
phyutility -concat -in gene* -out temp2.nex; mv gene* ../Used_Gene_Files; mv temp* ..; cd ../3
phyutility -concat -in gene* -out temp3.nex; mv gene* ../Used_Gene_Files; mv temp* ..; cd ../4
phyutility -concat -in gene* -out temp4.nex; mv gene* ../Used_Gene_Files; mv temp* ..; cd ..
rm -r 1 2 3 4

```

Now that you have your sequences in multiple (4 in this case) NEXUS files, it's time to concatenate them into a single one to make sure that you don't have alleles from the same locus as different sequences.

These files are big (~130-200 sequences and ~50Mb) and concatenating them using Phyutility takes **way** too long! To solve this, I tried two approaches. The first one was to import the NEXUS files into Geneious and concatenate them there. This was easy and took only a minute. I then exported the concatenated file in NEXUS and fasta format.

An alternative approach for non-Geneious users can be done in **biopython** with the code below, which I found [here](http://biopython.org/wiki/Concatenate_nexus). I wrote the code in a file an executed from the command line: `python Name_of_file.py`

```{python}
from Bio.Nexus import Nexus

file_list = ['temp1.nex', 'temp2.nex', 'temp3.nex', 'temp4.nex']
nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]

combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('All_Sequences_Concatenated.nex', 'w'))
```

### After you have a single (or multiple) Nexus file with all the exons concatenated

These two methods create a matrix as a NEXUS files (~1.1Gb!) where every sequence has the same length and, thus, adds a dash "-" to sites/regions that are missing (i.e., it treats those as missing data in a phylogenetic analysis). The file has to be striped of those dashes and modified to only have the fasta information. First, convert the NEXUS file to a fasta file using [NCLconverter](http://ncl.sourceforge.net), you can skip this step if you exported directly to fasta from Geneious. The code below will also compress the nexus file (from 1Gb to 1Mb) and **_delete_** the original file. **Don't leave huge files laying around in your computer, they eat up disk space quickly!**

```{bash}
for i in *.nex; do NCLconverter $i -efasta -o$i; done
for i in *nex; do tar -czf $i.tar.gz $i; rm $i; done

```
Now that you have a fasta, get rid of the dashes (or question marks, depending on the conversion method used) in the sequences, delete any spaces (" ") in the species names that NCLconverter added and replace them for underscores ("_"), and delete any possible empty lines. This will reduce the file size from ~1Gb to ~1Mb! The code below **will overwrite the original file**. If you don't want this, remove the "-i" option and specify an output file.

```{bash}
for i in *.fasta; do sed -i "" -e 's/-//g' -e 's/?//g' -e 's/ /_/g' -e '/^$/d' $i ; done

```

**Totally Optional**

These two methods produce interleaved fasta files, which work in every program downstream but I don't like. If you, like me, want to change these files to non-interleaved (one line per sequence) file, use the code below:

```{bash}
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file_in.fasta > file_out.fasta
```


OK, we are done. Now rename your concatenated, processed fasta file to `Targets_Final.fasta`. This file now contains all of your genes with the exons concatenated in a single sequences. This is the file that you'll use as the "Target File" (`targets.fasta` in the tutorial) in HypPiper. The last step to finish this file is to be sure that it matches the format that HybPiper requires, which is: `>SpeciesName-GeneName`. It shouldn't have any spaces and, apparently, not underscores. This works `>Text-GeneSomeNumber` but this doesn't `>Text_Second_Text_Numbers`. Do this with a few search and replace steps in TextWrangler or other text editor.

### Name List

In order to run HybPiper on many samples with a single command, it is necessary to have a file with the names of every sample (one per line). This was easy because this information is already in the Sample Sheet file we sent to RapidGenomics. Simply copy the column with the sample names (just one per sample) and paste it in a text file. Called the file `Name_List.txt`. In my case, because I have three clean read files per sample, the names of the sample include a `_Cleaned_` at the end. This allows me to call either the paired or the unpaired reads with the correct command in HybPiper (more on this later).

## Combine the 100bp and the 150bp cleaned reads
I created a file with the names of all the 100bp and 150bp files. To do this I did the following within the directory with the cleaned reads: `ls > 100bp_names.txt`. The same with the 150bp and then copied paste both lists into BBEdit and added some commands and modified the text with regular expression. A few lines of the file look like this:

```
#!/bin/bash
cat 37_C_erythraeus_Cleaned_PE1.fastq 37_C_erythraeus_R1.fastq.gz_Cleaned_PE1.fastq > 37_C_erythraeus_Combined_Cleaned_PE1.fastq
cat 37_C_erythraeus_Cleaned_PE2.fastq 37_C_erythraeus_R1.fastq.gz_Cleaned_PE2.fastq > 37_C_erythraeus_Combined_Cleaned_PE2.fastq
cat 37_C_erythraeus_Cleaned_SE.fastq 37_C_erythraeus_R1.fastq.gz_Cleaned_SE.fastq > 37_C_erythraeus_Combined_Cleaned_SE.fastq
cat 99_S_brevicalyx_Cleaned_PE1.fastq 99_S_brevicalyx_R1.fastq.gz_Cleaned_PE1.fastq > 99_S_brevicalyx_Combined_Cleaned_PE1.fastq
cat 99_S_brevicalyx_Cleaned_PE2.fastq 99_S_brevicalyx_R1.fastq.gz_Cleaned_PE2.fastq > 99_S_brevicalyx_Combined_Cleaned_PE2.fastq
cat 99_S_brevicalyx_Cleaned_SE.fastq 99_S_brevicalyx_R1.fastq.gz_Cleaned_SE.fastq > 99_S_brevicalyx_Combined_Cleaned_SE.fastq
cat 122_S_lycioides_Cleaned_PE1.fastq 122_S_lycioides_R1.fastq.gz_Cleaned_PE1.fastq > 122_S_lycioides_Combined_Cleaned_PE1.fastq

```

Now, in a directory with both the 100bp and 150bp reads and the file with the concatenation instructions:

```
mkdir Combined_Data 100bp
./Combining_Reads.sh

mv *Combined_Cleaned*.fastq ./Combined_Data
mv *R1.fastq.gz_Cleaned* ./150bp/Cleaned_Data_SeqyClean_May_2017/
mv *_Cleaned_* ./100bp

```



## Running HybPiper

The [tutorial](https://github.com/mossmatters/HybPiper/wiki/Tutorial) is very thorough and all the details can be found there. The only thing worth mentioning is that if you are using three files for your reads (PE1, PE2, and SE [single reads without a sister]), you have to use the option `--unpaired File_with_SE`

So the command for running HybPiper would look something like this:

```{bash}
/path_to_HybPiper_directory/reads_first.py -b /path_to_Targets_Final.fasta -r /path_to_Paired_reads*.fastq --prefix SampleName --bwa --unpaired /path_to_Unpaired_reads.fastq
```
(--cpu N limits the number of processors, the default is all).

### Running multiple samples with one command

To run multiple samples with a single command, you'll need the `Name_List.txt` and the following command:

```{bash}
while read name;
do /path_to_HybPiper_directory/reads_first.py -b /path_to_Targets_Final.fasta -r $name"P*.fastq" --prefix $name --bwa --unpaired $name"S*.fastq"
done < Name_List.txt
```

The code above reads the first line of the `Name_List.txt` and uses it to start the run for that sample, each line is read and store in the variable `name`. That's why it was important for me to add the `_Cleaned_` to the names so that the variable `name` becomes `Species_1_Cleaned_` and I could add either P for paired reads or S for unpaired reads. The only problem/limitation with this code is that it runs sequentially and not in parallel.

### To run multiple samples in parallel in a cluster using SLURM
This code will run an array of 152 jobs in the cluster. The indices of the array are used to read one line of the Name_List.txt at a time and submit a job for that sample/line. It will submit jobs in batches of 30 using 8 processors for each job. Invaluable information to read the names of the files with the indices came from [here](https://research.csc.fi/taito-array-jobs).  

**IMPORTANT**  

If you run the code below in a cluster with parallel computing, you might experience some problems if the cluster can't deal with the large number of files that are read and written. This will result in the recovery of *very very* few loci. For example, I was recovering only **12** loci out of 745 possible genes for a specific sample after running it on the cluster. I ran the same sample on my laptop and recovered **675!**. It's worth it to test one or a few samples on a standalone computer and then move to the cluster.

```{bash}
#!/bin/bash

#This will run an array of 152 jobs. The indices of the array are used to read one line of the Name_List.txt at a time and submit a job for that sample.

#SBATCH -J HybPiper
#SBATCH -o info.output.o%j
#SBATCH -e info.error.e%j
#SBATCH --array=1-152%30


#Set the long queue
#SBATCH -A long

#Set amount of time to run
#SBATCH -t 7-00:00:00

#Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cpus
#SBATCH -n 8

module load hybpiper
module load openmpi

# Run the MPI program

#The option below doesn't work
#name=$(tail -n $PBS_ARRAYID Name_File.txt | head -1)

#Use the array's index to read the file with names into the $name variable
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p Name_List.txt)

#Run HybPiper
reads_first.py -b Targets_Final.fasta -r $name"P*.fastq" --prefix $name --bwa --unpaired $name"S*.fastq" --cpu 8
```

### Get the length of your sequences with this code in the cluster

```
#!/bin/bash

#SBATCH -J HybPiper
#SBATCH -o info.output.Visualize
#SBATCH -e info.error.Visualize

#Set the long queue
#SBATCH -A long

#Set amount of time to run
#SBATCH -t 7-00:00:00

#Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cpus
#SBATCH -n 8

module load hybpiper
module load openmpi
module load python

# Run the MPI program

#Visualize results
#Get sequence length to visualize HybPiper's success at recovering genes at each locus, across all the samples.
python /opt/modules/biology/hybpiper/1.2/bin/get_seq_lengths.py Targets_Final.fasta Name_List.txt dna > Seq_Length.txt
```

### Compute statistics about your experiment
```
#!/bin/bash

#SBATCH -J HybPiper
#SBATCH -o info.output.Stats
#SBATCH -e info.error.Stats

#Set the long queue
#SBATCH -A long

#Set amount of time to run
#SBATCH -t 7-00:00:00

#Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cpus
#SBATCH -n 8

module load hybpiper
module load openmpi
module load python

# Run the MPI program

#Produce Stats
#This script will summarize target enrichment and gene recovery efficiency for a set of samples.
python /opt/modules/biology/hybpiper/1.2/bin/hybpiper_stats.py Seq_Length.txt Name_List_Taxa_without_Genes_Deleted.txt > Stats.txt
```

### Teasing *introns* and *exons* apart and getting their sequences

```
#!/bin/bash

#This will run an array of 152 jobs. The indices of the array are used to read one line of the Name_List.txt at a time and submit a job for that sample.

#SBATCH -J HybPiper
#SBATCH -o info.output.o%j
#SBATCH -e info.error.e%j
#SBATCH --array=1-152%30


#Set the long queue
#SBATCH -A long

#Set amount of time to run
#SBATCH -t 7-00:00:00

#Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cpus
#SBATCH -n 8

module load hybpiper
module load openmpi

# Run the MPI program

#Use the array's index to read the file with names into the $name variable
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p Name_List.txt)

#Run HybPiper
python /opt/modules/biology/hybpiper/1.2/bin/intronerate.py --prefix $name

#Or run this line in a sequential way:
while read name;
do python /opt/modules/biology/hybpiper/1.2/bin/intronerate.py --prefix $name
done < Name_List.txt

```

### Retrieve all the genes that you recovered, either nucleotides (`dna`), proteins (`aa`), introns and exons (`supercontig`), or just introns (`intron`)
```
#!/bin/bash

#SBATCH -J HybPiper
#SBATCH -o info.output.Retrieving
#SBATCH -e info.error.Retrieving

#Set the long queue
#SBATCH -A long

#Set amount of time to run
#SBATCH -t 7-00:00:00

#Set the number of nodes
#SBATCH -N 1

# Set the number of tasks/cpus
#SBATCH -n 8

module load hybpiper
module load openmpi
module load python/2.7.3

# Run the MPI program

#Retrieving Sequences
#This script fetches the sequences recovered from the same gene for many samples and generates an unaligned multi-FASTA file for each gene.
#Change dna below for supercontig, contig, or aa to retrieve different type of data.

python /opt/modules/biology/hybpiper/1.2/bin/retrieve_sequences.py Targets_Final.fasta . dna > exons.log
python /opt/modules/biology/hybpiper/1.2/bin/retrieve_sequences.py Targets_Final.fasta . aa > proteins.log
python /opt/modules/biology/hybpiper/1.2/bin/retrieve_sequences.py Targets_Final.fasta . intron > introns.log
python /opt/modules/biology/hybpiper/1.2/bin/retrieve_sequences.py Targets_Final.fasta . supercontig > supercontigs.log

mkdir Final_Genes
mkdir ./Final_Genes/Exons ./Final_Genes/Proteins ./Final_Genes/Supercontig ./Final_Genes/Introns
mv *FNA exons.log ./Final_Genes/Exons
mv *introns* ./Final_Genes/Introns
mv *FAA proteins.log ./Final_Genes/Proteins
mv *supercontig* ./Final_Genes/Supercontig

```

### Analyzing potential paralog loci

Matt has a really good explanation on what you can expect and do with the possible paralog gene. Go read it [here](https://github.com/mossmatters/HybPiper/wiki/Paralogs)

To get a list of the number of copies per gene in each sample to stdout, type the following:

```
while read i
do
echo $i
python /opt/modules/biology/hybpiper/1.2/bin/paralog_investigator.py $i
done < Name_List.txt
```

You'll need to select all and copy/paste the list to a text file.

Now you have to create a file with the genes that have potential paralogs. Open the file you created in the previous step (number of copies per gene per sample) and remove all the information from it **except** for the gene name. Do this in a text editor with search and replace.

Example: this line of text `2 paralogs written for Gene000002284419` changes to `Gene000002284419`. After you've done this, delete any duplicated lines. In BBEdit go to "Text->Process Duplicated Lines". Finally, you need to delete all returns `\r` and replace them for spaces ` ` so that all the genes are in a single line.

Now you have a list of genes that have multiple copies.

Finally run the following command:


```
# Using a file with one gene per line:
mkdir Final_Genes/Paralogs
parallel "python /opt/modules/biology/hybpiper/1.2/bin/paralog_retriever.py Name_List.txt {} > {}.paralogs.fasta" :::: ./Final_Genes/Paralog_Genes.txt 2> ./Final_Genes/Paralog_table.txt
mv *paralogs* Final_Genes/Paralogs
```


```
# Using the gene names
mkdir Final_Genes/Paralogs
parallel "python /opt/modules/biology/hybpiper/1.2/bin/paralog_retriever.py Name_List.txt {} > {}.paralogs.fasta" ::: Gene000000000064 Gene000000000099 Gene000000000106 Gene000000000107 Gene000000000108 Gene000000000111 Gene000000000119 Gene000000013958 Gene000000026916 Gene000000034941 Gene000000072925 Gene000000088452 Gene000000093415 Gene000000095243 Gene000000105734 Gene000000109936 Gene000000117000 Gene000000139714 Gene000000140381 Gene000000148794 Gene000000151236 Gene000000155182 Gene000000155339 Gene000000166541 Gene000000168219 Gene000000179355 Gene000000185393 Gene000000185460 Gene000000198698 Gene000000210298 Gene000000213937 Gene000000216209 Gene000000226910 Gene000000252250 Gene000000261154 Gene000000272966 Gene000000284062 Gene000000314285 Gene000000370726 Gene000000372933 Gene000000408989 Gene000000456987 Gene000000467115 Gene000000480359 Gene000000510272 Gene000000520052 Gene000000560456 Gene000000630862 Gene000000665080 Gene000001291309 Gene000001298132 Gene000002089674 Gene000002090980 Gene000002093127 Gene000002162473 Gene000002165004 Gene000002165930 Gene000002168236 Gene000002182899 Gene000002213682 Gene000002219132 Gene000002219408 Gene000002221889 Gene000002240027 Gene000002281776 Gene000002284419 Gene000002305250 Gene000002341791 Gene000002346357 Gene000002348025 Gene000002349895 Gene000002410184 Gene000002412774 Gene000002413865 Gene000002420086 Gene000002420898 Gene000002437988 GeneAT1G01180 GeneAT1G04110 GeneAT1G04640 GeneAT1G07740 GeneAT1G10460 GeneAT1G15510 GeneAT1G21840 GeneAT1G26900 GeneAT1G32520 GeneAT1G53600 GeneAT1G65380 GeneAT2G16880 GeneAT2G28790 GeneAT2G29900 GeneAT2G42700 GeneAT2G43760 GeneAT3G05340 GeneAT3G06920 GeneAT3G21470 GeneAT3G48150 GeneAT3G58520 GeneAT4G21300 GeneAT4G23890 GeneAT4G31850 GeneAT4G33990 GeneAT5G06370 GeneAT5G18390 GeneAT5G18475 GeneAT5G21970 GeneAT5G43280 GeneAT5G50290 2> ./Final_Genes/Paralog_table.txt
mv *paralogs* Final_Genes/Paralogs
```

### Separate by Clade

I'm going to separate the final sequences (the supercontigs in this case) in each fasta file by taxonomic group.

Start with the **non-interleaved** fasta files that HybPiper created:

```{bash}
cd Supercontig
mkdir Non-Interleaved
for i in *fasta; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > $i"_non_interleaved.fasta"; done
mv *non-interleaved.fasta Non-Interleaved/

```

For the *Burmeistera* species that have good names and enough loci, plus some outgroups:

```
for i in *.aln; do grep -A 1 -e \>Burmeistera_almedae_2016_001_Combined -e \>Burmeistera_anderssonii_2016_113_Combined -e \>Burmeistera_aspera_2016_089_Combined -e \>Burmeistera_auriculata_2016_199_Combined -e \>Burmeistera_borjensis_2016_219_Combined -e \>Burmeistera_brachyandra_2016_114_Combined -e \>Burmeistera_brighamioides_2016_121_Combined -e \>Burmeistera_bullatifolia_2016_129_Combined -e \>Burmeistera_ceratocarpa_2016_233_150bp -e \>Burmeistera_cf_aeribacca_2016_087_Combined -e \>Burmeistera_chiriquiensis_2016_161_Combined -e \>Burmeistera_crassifolia_2016_200_Combined -e \>Burmeistera_crispiloba_2016_201_Combined -e \>Burmeistera_cyclostigmata_2016_104_Combined -e \>Burmeistera_cylindrocarpa_2016_235_Combined -e \>Burmeistera_darienensis_2016_149_Combined -e \>Burmeistera_domingensis_2016_192_Combined -e \>Burmeistera_draconis_2016_133_Combined -e \>Burmeistera_dukei_2016_153_Combined -e \>Burmeistera_fuscoapicata_2016_128_Combined -e \>Burmeistera_glabrata_2016_227_Combined -e \>Burmeistera_holm-nielsenii_2016_210_Combined -e \>Burmeistera_huacamayensis_2016_115_Combined -e \>Burmeistera_litensis_2016_214_Combined -e \>Burmeistera_loejtnantii_2016_006_Combined -e \>Burmeistera_lutosa_2016_209_Combined -e \>Burmeistera_mcvaughii_2016_155_Combined -e \>Burmeistera_multiflora_2016_196_Combined -e \>Burmeistera_obtusifolia_2016_108_Combined -e \>Burmeistera_oyacachensis_2016_230_Combined -e \>Burmeistera_panamensis_2016_151_Combined -e \>Burmeistera_parviflora_LL_20_Combined -e \>Burmeistera_pirrensis_2016_148_Combined -e \>Burmeistera_quercifolia_2016_009_Combined -e \>Burmeistera_ramosa_2016_007_Combined -e \>Burmeistera_refracta_2016_195_Combined -e \>Burmeistera_resupinata_2016_240_Combined -e \>Burmeistera_resupinata_var_heilbornii_2016_116_Combined -e \>Burmeistera_rubrosepala_2016_239_Combined -e \>Burmeistera_smaragdi_2016_236_Combined -e \>Burmeistera_smooth_2016_135_Combined -e \>Burmeistera_sodiroana_2016_197_Combined -e \>Burmeistera_succulenta_2016_143_Combined -e \>Burmeistera_succulenta_2016_188_Combined -e \>Burmeistera_tenuiflora_2016_158_Combined -e \>Burmeistera_truncata_2016_144_Combined -e \>Burmeistera_utleyi_2016_156_Combined -e \>Burmeistera_vulgaris_2016_100_Combined -e \>Burmeistera_xerampelina_2016_086_Combined -e \>Burmeistera_zurquiensis_2016_098_Combined -e \>LL6_C_smithii_Combined -e \>LL49_C_incanus_Combined -e \>LL61_C_asclepideus_Combined -e \>LL69_C_nigricans_Combined -e \>LL83_C_brittonianus_Combined -e \>LL86_C_mandonis_Combined -e \>LL88_S_ayersiae_Combined -e \>LL159_S_jelskii_Combined -e \>LL334_S_aureus_Combined -e \>LL363_S_krauseanus_Combined  $i > Species_with_Good_Data/$i"_Good_Data.aln"; done

#Clean the two dashes in extra lines

cd Species_with_Good_Data
for i in *.aln; do sed -i "" 's/^--$//g' $i; done


```


For **all** *Burmeistera*, including some outgroup species:

```
cd Non-Interleaved
mkdir Burmeistera_and_Outgroups

for i in *non_interleaved.fasta; do grep -A 1 -e \>Burm -e \>LL69_C_nigricans -e \>LL363_S_krauseanus -e \>LL6_C_smithii -e \>LL159_S_jelskii -e \>LL61_C_asclepideus -e \>LL334_S_aureus -e \>LL49_C_incanus -e \>LL83_C_brittonianus -e \>LL86_C_mandonis -e \>LL88_S_ayersiae $i > Burmeistera_and_Outgroups/$i"_Burmeistera_and_Outgroups.fasta"; done

#Clean the two dashes in extra lines

cd Burmeistera_and_Outgroups
for i in *.fasta; do sed -i "" 's/^--$//g' $i; done

```

### Alignment, Clean up, and Phylogenetics

**Before** running this code, read how to "parallelize" the work load below!

If you don't have access to a cluster with many nodes and cores but have access to a few machines with many core, you can still work on parallel by assigning some work to each computer.

You need to figure out how many files you have, how many "batches" you are going to divide the work in, move the right amount of files into each batch, and run the analyses in each batch.

Using the same code I used during concatenation step if the target genes:

```
#Delete any empty files
find . -size 0 -delete

#How many batches (directories) do you want to work with?
DIRECTORIES=30

#Counts number of files in a directory and stores it in a variable
NUMFILES=$(ls | wc -l)
#Divides the variable by a number specified in a variable called `DIRECTORIES`, which is the number of files I wanted to work with. It stores this number and its decimals (files/$DIRECTORIES) into a variable
NUMFILESDIR=$(echo "$NUMFILES/$DIRECTORIES"|bc -l)
#Divides the variable by four, which is the number of files I wanted to work with. It stores this number with NO decimals into a variable
NUMFILESDIR_NO_FLOAT=$(($NUMFILES / $DIRECTORIES))
#Informs how many files need to go to each directory
echo "<< The number of files per directory is ***$NUMFILESDIR*** >>"

eval mkdir {01..$DIRECTORIES}

#Move the corresponding number of files into each directory (number of lines must match $DIRECTORIES)

for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./01; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./02; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./03; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./04; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./05; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./06; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./07; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./08; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./09; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./10; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./11; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./12; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./13; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./14; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./15; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./16; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./17; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./18; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./19; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./20; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./21; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./22; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./23; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./24; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./25; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./26; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./27; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./28; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./29; done
for file in $(ls -p | grep -v / | tail -$NUMFILESDIR_NO_FLOAT); do mv $file ./30; done

## Check that you moved all the files!!!!!
## If not, move the remainder manually to the last (largest number) directory.

```

Within the directory with the cleaned fasta files, either execute the code below or, even better, copy paste it into a file and run it as a script. I called my script `Align_Clean_Phylo.sh`, and I copied the script into **in each** directory/batch of files with the code below:

```
echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 | xargs -n 1 cp Align_Clean_Phylo.sh

```

The run the script with `./Align_Clean_Phylo.sh`

#### Code for phylogenetics
```
#!/bin/bash

# This script will align, clean, and analyze multiple fasta files in a single
# directory. It assumes that you have Mafft, Phyutility, and FastTree installed
# and in your path. The first part of the script is meant to be used in a cluster
# where modules have to be loaded, so delete it or comment it out if you are working
# on a single machine. The fasta files, alignments, and trees will be placed in 
# separate directories.

# Simon Uribe-Convers - http://simonuribe.com - June 15th, 2017

# Load modules if running on a cluster

module load mafft
module load phyutility
module load fasttree

# Delete any empty files

find . -size 0 -delete

# Create directories

mkdir Alignments Fasta Phylo

## Align with Mafft, clean with Phyutility, and analyze with FastTree

# Alignment
for i in *.fasta
do
# Align
echo ""
echo "~~~Aligning with Mafft~~~"
echo ""
time mafft --auto --thread 6 --preservecase $i > $i.aln

# Cleaninig Alignment at 50% occupancy per site
echo ""
echo "~~~Cleaninig alignment with Phyutility at at 50% occupancy per site~~~"
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

### If you need to standardize names within the trees

Maybe I'll need to make the tips match for species tree analyses or concatenation. The code below will get rid off the appended gene information:

```
perl -pe 's/_Gene\w+//g' InFile > OutFile
```

## Species Tree Inference

After all the gene trees are built, it's time to make a species tree under the multispecies coalescent. I'm going to use [ASTRAL II](https://github.com/smirarab/ASTRAL) and [SVDquartets](http://www.stat.osu.edu/~lkubatko/software/SVDquartets/) for this purpose.

### ASTRAL II  

Starting in a directory with all the gene trees (the best tree from RAxML for example), make a file with all the gene trees in it:

```
cat genetree* > All_gene_trees.tre
```

Make sure that the names of the species in each gene tree are the **same!!** Otherwise the species tree program won't be able to know they are the same. This can be skipped if you standardized the names with the perl step above.

For example, "species1_gene1" for gene tree 1 needs to be changed to "species1". This is easily done in a text editor.

Now that they all match, run ASTRAL II

```
java -jar /Applications/ASTRAL-master/Astral/astral.4.11.1.jar -i All_gene_trees.tre -o Species_tree_Astral.tre

```

### SVDquartets

Start with a file with all loci (NSPs) concatenated into a NEXUS file, you will use this in Paup.

SVDquartets allows for multiple individuals from the same species to be included in the analysis, and those individual, and their information, will be combined into the taxon they belong to. For this to work however, you need to include a taxon block on your NEXUS file specifying which samples belong to which species. It should start with the species name, followed by a colon, the lines the individuals are located (or the range of lines) and a comma. Here is an example:

```
begin sets;
taxpartition species =
sp1: 1-4,
sp2: 5,
sp3: 6-7,
;
END;

``` 
This is easy to do for a few samples but if you have hundreds of individuals it becomes tedious very quickly. To make it a bit easier, use the R code below. **Note:** this code works on a *Phylip* file and not a *NEXUS*, so convert the NEXUS to Phylip using NCLconverter: `NCLconverter infile.nex -erelaxedphylip -ofileout`.


```{R}
### This script will format the settings file for SVDquartets
### It will parse a phylip file and output the line number in which each taxon is located. Then it will write how many
### occurrences a specific species has and the lines of each.
### This works best with species names separated by an underscore, and it will assume that there are no subspecies,
### i.e., only the fisrt two parts of the name will be used.
### by Matthew Pennell, July 23 2014 - http://mwpennell.github.io/

## Read in and parse phylip file
phy.tab <- read.table(phylip.filename)
phy.tab <- as.character(phy.tab[-1,1])

## specific to my naming scheme using underscores
get.species.name <- function(x){
  tmp <- strsplit(x, split="_")
  paste(tmp[[1]][1], tmp[[1]][2], sep="_")
}

## Get species names 
sp.names <- as.character(sapply(phy.tab, function(x) get.species.name(x)))

## get identity of matches
sp.lab <- lapply(unique(sp.names), function(i) which(sp.names == i))
names(sp.lab) <- unique(sp.names)

## output to input file
sink("SVDquartets_settings.txt")
for (i in 1:length(sp.lab)){
	cat(paste(names(sp.lab)[i], length(sp.lab[[i]]), sep=" : "))
	cat("\n")
	cat("\t")
	cat(as.numeric(sp.lab[[i]]))
	cat("\n")

}
sink()
```
One you have the location of every sample from the code above, modified the text slightly to match the correct format, i.e., put every occurrence in one line, add commas, etc. Finally, copy paste your sample-to-species information at the end of your concatenated NEXUS file—don;t forget to include the few lines that the code above doesn't generate, see the format example above!

Now that we have the file ready, get the latest command-line version of [Paup](https://people.sc.fsu.edu/~dswofford/paup_test/) and type the following:


```
paup

#Within Paup

exe All_Loci_Concatenated_Good_Data_SVDquartets.nex

#No Bootstrap
SVDQuartets nthreads=25 nquartets=10000000 partition=species speciesTree=yes

#Save the tree
SaveTrees file = Burmeistera_SVDquartet.tre format = Newick brLens = yes supportValues = Both trees = all

#With Bootstrap
SVDQuartets nthreads=25 nquartets=10000000 partition=species bootstrap=standard speciesTree=yes nreps=50 treeFile= Burmeistera_SVDquartet_Bootstrap3_Good_Data.trees


```


### Cleaning up

run the code below to delete unnecessary files, mostly the results from SPAdes. This will reduce to size of the directories by 75%.

```
python /opt/modules/biology/hybpiper/1.2/bin/cleanup.py
```

### Final Commands and Notes

- **Remember** to compress your reads after you are done!

```{bash}
gzip ./Cleaned_Data_SeqyClean_May_2017/*fastq
```

- Be in the directory that contains the `Name_List.txt` and `Targets_Final.fasta`. The reads/genes will be processed here.
- Make sure that the paths on your commands point to the directory of HybPiper and the appropriate reads. **Important** If the cluster/computer you are using gives you an error that the read files don't exist, even though you triple checked that the path is correct and that there are not typos in the names, you will need to move the reads to the directory you are processing the genes in.
- Use the flag `--cpu N` to specify the number of processors available to HybPiper. This is especially important if you are using a cluster or if you want to limit the available resources of your computer.
- **Have fun!**
