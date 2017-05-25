# Processing Sequence Capture Data, version 2

This time we got the right data for our Burmeistera capture sequencing experiment, paired-end 150bp reads. This are the steps I took to process the data.

## Renaming of the reads
Starting from the spreadsheet with the sample names, copy and paste the columns with the samples and the corresponding sequencing file name on a text editor. Copy paste these two columns again in a new file. The first file will be for the forward reads and the second for the reverse reads. 

Do some text parsing with search and replace.

**For the forward reads:**

```
Search:
^(.+)\t(.+)
Replace:
mv \2.gz \1_R1.fastq.gz

```

**For the reverse reads:**  

```
^(.+)\t(.+)
Replace:
mv \2.gz \1_R1.fastq.gz

Search:
_R1\.
Replace:
_R2.
```

## Cleaning the data


Cleaning the data with SeqyClean:

```
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

The next instruction will be brief because they are the same as in the last post where I processed the 100bp reads.

Once you have the necessary files to run HybPiper, i.e., the targets file and the names file, the basic command for one sample is:

```
reads_first.py -b Targets_Final.fasta -r Burmeistera_brachyandra_2016_114_Combined_Cleaned_PE* --prefix Burmeistera_brachyandra_2016_114_Combined --bwa --unpaired Burmeistera_brachyandra_2016_114_Combined_Cleaned_SE.fastq --cpu 8
```
