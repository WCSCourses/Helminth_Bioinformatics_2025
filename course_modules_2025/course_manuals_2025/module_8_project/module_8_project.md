# Project 

Congratulation!! You have reached the final day of the course. This is where you will gather the skills you have learnt and practice them on a new problem. Exchange the knowledge, tips and solutions with your coursemates, and ask your instructors if you get stuck.  


There are two projects to choose from:
1. Transcriptomics of the filarial nematode *Brugia malayi* and its wolbachia endosymbiont
2. Analysis of cryptic genetic variation in a family of *Haemonchus contortus*


--



## Project 1: Transcriptomics of the filarial nematode *Brugia malayi* and its wolbachia endosymbiont

There are three directories with files that are relevant to this project. **Don't get scared with the amount of files**

**_Brugia_**
- The only directory you really need to get started is the **Brugia**. This directory contains everything you need to look at gene expression in the _B. malayi_ worms. 
```
./RNAseqData_Brugia/Bm_counts:
F120a.counts  F30a.counts  F42a.counts  L4a.counts  M120a.counts  M30a.counts  M42a.counts
F120b.counts  F30b.counts  F42b.counts  L4b.counts  M120b.counts  M30b.counts  M42b.counts
```

**_Wol_**
- The **Wol** directory is only needed if you want to do dual-RNA seq and to look at bacteria gene expression as well. 
```
./RNAseqData_Wol:
Wol_F120a_htseqcount.txt  Wol_F30a_htseqcount.txt Wol_M120a_htseqcount.txt  Wol_M30a_htseqcount.txt
Wol_F120b_htseqcount.txt  Wol_F30b_htseqcount.txt Wol_M120b_htseqcount.txt  Wol_M30b_htseqcount.txt
```

**Note: There is no need to download the FASTQ files for this project session because**
- 1)	These FASTQ files are HUGE!! To process these files, we normally work on a computer cluster or a high-performance computing server, not on individual laptop
- 2)	We have already downloaded and processed the FASTQ files for you, up to the point that the data can be analysed on a laptop. These are the read count files you see in your downloaded files. 

- There are also read count files and a text file containing the accession IDs of the original data.
- The accession IDs are for the original FASTQ files, and the IDs are there for your reference and for completeness only. 


### The steps that were done for you are
1)	Download the FASTQ files from a public repository (see appendix [Finding and download sequence data from public repository](https://wcscourses.github.io/HelminthBioinformatics_2023/manuals/other_information/Sequence_data_on_public_repo.html) if you want to learn how to do this)

2)	Download the reference genome and annotation file (GTF) of the worms from WormBaseParaSite. Additionally, download the reference genome of Wolbachia of _B. malayi_ from [Ensembl Bacteria](http://bacteria.ensembl.org/index.html). The usage is very similar to WormBaseParaSite but it contains bacterial reference genomes. 

3)	Index the reference genome using `hisat2` (see transcriptomics module)

4)	Mapped the FASTQ to the reference genome using `hisat2` (see transcriptomics module)

5)	Convert SAM to BAM using `samtools view` (see transcriptomics module)

6)	Sort BAM files by read names using `samtools sort` (see transcriptomics module)

7)	Count reads per features (in this case, reads per gene) using `htseq-count` (an alternative to featureCounts tool)


### The steps that you will perform by yourself in this project are (most of them refer to the transcriptomics module):
1)	Starting in RStudio

2)	Import the read counts data into RStudio

3)	Explore the data using various data visualisation

4)	Perform pair-wise comparison of any comparisons that you are interested in and produce relevant plots

5)	Look at your results, compare your results, follow the lead of your results, explore genes using public databases and tools e.g. WormBaseParaSite, InterProScan, AlphaFold etc. 

6)	Perform functional analysis using GO term enrichment

### GO enrichment tips!!!! 
Notice that we have provided a _GO term annotation reference_ for you in the [transcriptomics module], but not in this project module. You cannot use the same GO annotation reference here because...

We are using a completely different organism here. Therefore, we need a new GO annotation for this one. 

You will need to download a new GO annotation by yourself using WormBaseParaSite BioMart. 🤓

- See appendix [Downloading and formatting GO annotation](https://wcscourses.github.io/HelminthBioinformatics_2023/manuals/other_information/GO_ref_download_and_formatting.html) on how to do this
- You are also provided with the R script `GO_formatting.R` to reformat the downloaded GO reference into the format the `topGO` require. **However, this is not a complete-and-ready-to-use script**, and it will need some editing _(a real life scenario when you borrow script from someone else)_ 
- Have a look through the script `GO_formatting.R` and try to figure out which parts need to be edited. A simple way to do this would be to copy the whole code from `GO_formatting.R` and place it into your current R script and test your editing there. 




---

[↥ **Back to top**](#top)
<br>
<br>

## Project 2: Analysis of cryptic genetic variation in a family of *Haemonchus contortus*
### Background
Using unpublished data from Steve's lab, we would like you describe genetic variation among a closely related group of *H. contortus* parasites. Whole genome sequencing data was generated from hundreds of individual worms from the 5th generation of a genetic cross between susceptible female worms and resistant male worms. These crosses were described in two papers [here](https://doi.org/10.1093/gbe/evx269) and [here](https://doi.org/10.1016/j.celrep.2022.111522) (which might be useful for background and for your presentation).

We have extracted the reads from the mitochondrial genome for you for mapping. 

In theory, mitochondrial DNA is inherited from the mother (however, this is not exactly known in nematodes) - even though this is the 5th generation of the genetic cross, we might expect that mitochondrial variation in the progeny should look very similar the the original susceptible female parent. Is this true? We want you to find out.



We have provided sequencing data for the following sample sets. You will see the following codes in the sample names, which will help you in your analysis:
- ```MHCO18``` is the resistant parent (n=1)
- ```MHCO3``` is the susceptible parent (n=1)
- ```XQTL``` is the progeny. Within the progeny, we have two main groups:
    - ```SUS``` - which are highly susceptible to ivermectin (n=180)
    - ```RES``` - which are highly resistant to ivermectin (n=28)




### The steps that you will perform by yourself in this project are
1. download the refererence from WormBase ParaSite, and extract the mitochondrial genome
2. download the raw sequencing data from Steve's FTP site. To help you, you can use the following command:

```bash
wget ftp://ngs.sanger.ac.uk//production/pathogens/sd21/helminth_bioinformatics_2023_variation_project/*gz

```

3. map the reads to the mitochondrial genome. Hint: modify the ```for loop``` in the genetic variaiton module. 
3. Call variants
4. perform a PCA analysis of SNP variants
5. Colour the PCA by different grouping factors - you will need to make a metadata table. Hint: look closely at the sample names
6. find specific variants that discriminate different groups

### Working in a group
- see if you can split the tasks to save time
- for example, you could have a "mapping team" and a "metadata" team



---

[↥ **Back to top**](#top)
<br>
<br>
