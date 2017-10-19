# Generate an inclusive matrix from multiple alignments

Say you have hundreds of alignments and you want to figure out which genes are shared by all of the species. Here are the steps how to do it.

## Species in each gene

First you need to know which species are present in each gene. The Python code below will do this:

```{python}

"""
Within a directory with alignment files in fasta format, this script will read 
every file and return the ID for that sequence (the species) that is present in 
that file, along with the file it belongs to. In other words, it will tell you 
what species have that genes and in which file they are. It will create a file 
for each gene with that information.

After this is done, use the file produced here with Matt Pennell's R script "Get_Inclusive_Matrix_for_Genes.R"
"""


import glob
from Bio import SeqIO
results = ""
for filename in glob.glob('*.aln'):
    # Uncomment the line below to separate each gene's information with a header
    #results = results + "\r" + "Species_Name" + "\t" + "Gene_Name = " + filename + "\r"
    for record in list(SeqIO.parse(filename, "fasta")):
        results = results + record.id + "," + filename +"\r"
    file = open("Genes_For_Each_Sample.csv","w")
    file.write(results)
    file.close()

```

## Generate a table with the inclusive matrix

Now that you know which species are where, you need to compare this information and generate a list with the inclusive data. The R code below does this:


```{R}
## By Matthew Pennell - June 22nd, 2017
# The code below will constrcut a matrix for species that have the same genes from a file with
# information on which species has which gene. 
# It creates an inclusive list with all species having the exact number of genes
# You can then use one of the rows to generate a list of species and grep those
# out of a set of alignments or fasta.



library(dplyr)
t <- read.csv("Genes_For_Each_Sample.csv", as.is=TRUE)
colnames(t) <- c("species", "gene")

## set threshold
thresh <- 0.8

n_gene <- t %>% count(gene) %>% nrow()
spp    <- t %>% count(species) %>% filter(n >= n_gene * thresh)
t_thr  <- t %>% filter(species %in% spp$species)

## most inclusive set of genes (no missing data)
most_incl <- t_thr %>% count(gene) 
## let's just say (arbitrarily) that you want at least 76 spp
## and we'll refilter. Modify depending on how many species
## there are in the dataset.
min_spp <- 60
gene_min <- most_incl %>% filter(n >= min_spp)
t_thr_min <- t_thr %>% filter(gene %in% gene_min$gene)

g_by_spp <- lapply(unique(t_thr_min$gene), function(x) {
  x <- t_thr_min[which(t_thr_min$gene == x), "species"]
  return(tbl_df(x))
})

## length
ll <- length(g_by_spp)

for (i in 2:ll){
  g_by_spp[[i]] <- dplyr::intersect(g_by_spp[[i]], g_by_spp[[i-1]])
}

out <- g_by_spp[[ll]]
nrow(out)

final <- t_thr_min %>% filter(species %in% out$value)
gene_count <- final %>% count(gene)
spp_count <- final %>% count(species)

f <- data.frame(matrix(nrow=nrow(gene_count), ncol=nrow(spp_count)))
rownames(f) <- unique(final$gene)
colnames(f) <- sapply(c(1:ncol(f)), function(x) paste0("spp", x))
for (i in 1:nrow(f))
  f[i,] <- unique(final$species)

write.csv(f, "fullmatrix_80.csv")

```

The code above generates a table with gene in rows and species in columns. Because it is inclusive (same data across the table), you only need one row to create a list of species that you can then grep out of the original fasta or alignment files. 

Copy the first row and make it a grep command in a text editor. For *Burmeistera* and outgroups, this is the command.

```{bash}
mkdir All_Species_All_Genes

for i in *aln
do
grep -A 1 -e \>Burmeistera_almedae_2016_001_Combined -e \>Burmeistera_anderssonii_2016_113_Combined -e \>Burmeistera_aspera_2016_089_Combined -e \>Burmeistera_auriculata_2016_199_Combined -e \>Burmeistera_borjensis_2016_219_Combined -e \>Burmeistera_brachyandra_2016_114_Combined -e \>Burmeistera_brighamioides_2016_121_Combined -e \>Burmeistera_bullatifolia_2016_129_Combined -e \>Burmeistera_cf_aeribacca_2016_087_Combined -e \>Burmeistera_chiriquiensis_2016_161_Combined -e \>Burmeistera_crassifolia_2016_200_Combined -e \>Burmeistera_crispiloba_2016_201_Combined -e \>Burmeistera_cyclostigmata_2016_104_Combined -e \>Burmeistera_cylindrocarpa_2016_235_Combined -e \>Burmeistera_darienensis_2016_149_Combined -e \>Burmeistera_domingensis_2016_192_Combined -e \>Burmeistera_draconis_2016_133_Combined -e \>Burmeistera_dukei_2016_153_Combined -e \>Burmeistera_fuscoapicata_2016_128_Combined -e \>Burmeistera_glabrata_2016_227_Combined -e \>Burmeistera_holm_nielsenii_2016_210_Combined -e \>Burmeistera_huacamayensis_2016_115_Combined -e \>Burmeistera_litensis_2016_214_Combined -e \>Burmeistera_loejtnantii_2016_006_Combined -e \>Burmeistera_lutosa_2016_209_Combined -e \>Burmeistera_mcvaughii_2016_155_Combined -e \>Burmeistera_multiflora_2016_196_Combined -e \>Burmeistera_obtusifolia_2016_108_Combined -e \>Burmeistera_oyacachensis_2016_230_Combined -e \>Burmeistera_panamensis_2016_151_Combined -e \>Burmeistera_parviflora_LL_20_Combined -e \>Burmeistera_pirrensis_2016_148_Combined -e \>Burmeistera_quercifolia_2016_009_Combined -e \>Burmeistera_ramosa_2016_007_Combined -e \>Burmeistera_refracta_2016_195_Combined -e \>Burmeistera_resupinata_2016_240_Combined -e \>Burmeistera_resupinata_var_heilbornii_2016_116_Combined -e \>Burmeistera_rubrosepala_2016_239_Combined -e \>Burmeistera_smaragdi_2016_236_Combined -e \>Burmeistera_smooth_2016_135_Combined -e \>Burmeistera_sodiroana_2016_197_Combined -e \>Burmeistera_succulenta_2016_143_Combined -e \>Burmeistera_succulenta_2016_188_Combined -e \>Burmeistera_tenuiflora_2016_158_Combined -e \>Burmeistera_truncata_2016_144_Combined -e \>Burmeistera_utleyi_2016_156_Combined -e \>Burmeistera_vulgaris_2016_100_Combined -e \>Burmeistera_xerampelina_2016_086_Combined -e \>Burmeistera_zurquiensis_2016_098_Combined -e \>LL6_C_smithii_Combined -e \>LL49_C_incanus_Combined -e \>LL61_C_asclepideus_Combined -e \>LL69_C_nigricans_Combined -e \>LL83_C_brittonianus_Combined -e \>LL86_C_mandonis_Combined -e \>LL88_S_ayersiae_Combined -e \>LL159_S_jelskii_Combined -e \>LL334_S_aureus_Combined -e \>LL363_S_krauseanus_Combined $i > All_Species_All_Genes/$i"_inclusive_data.aln"

# The line below deletes the weird dashes that appear after after each sequence
sed -i "" 's/^--$//g' All_Species_All_Genes/$i"_inclusive_data.aln"

done 

```

Now you can build trees! 