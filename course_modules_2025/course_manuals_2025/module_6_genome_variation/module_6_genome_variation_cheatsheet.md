# Genetic variation module - cheat sheet

#### 2.2. Questions:
- **How does the "all sample" report compared to the "Australian-only sample" report?**
        - The all sample report will obviously contain a lot more samples, and should be possible to see some variaiton between samples, including variaiton between groups of samples per country.


- **From the "General statistics" section, can we see any sample groups that look different and that might be problematic for our analyses?**
    - some of the samples from Benin, Brazil, and Cape Verde have high duplication rates

- **Are all of the read lengths the same? How can you tell?**
    - no, there are some differences in read lengths
    - comparing read lengths, esp if trimming of reads is invovled is useful

- **Can you think of any other quality control checks on the raw data that you might want to do before analysing your data? Hint: think about where these worms are collected from.**
    - contamination screen using a tool like kraken would be useful. Helminths often contain host contamination, however, they may be heavily contaminated with bacteria. If most of the reads are contamination, could be worth stopping
    - we could screen and remove adapters and poor quality sequences. Using a tool like trimmomatic or similar to trim the reads can be useful. Poorly constructed libraries, esp those with very small insert sizes, will have adapator contamination which will affect how reads map.



#### 3.1. Questions:
- **Looking at the WormBase ParaSite website for *Haemonchus contortus* - how big is the genome? how many genes are present?**
    - genome is 283 Mb
    - annotation has 	19,623 genes


- **There is a second *Haemonchus contortus* genome resource also present - how do the two genomes compare?**
    - the second Haemonchus genome is from the University of Melbourne, under accession 	PRJNA205202
    - it is a little bigger (319 Mb) and has more genes (23,610)
    - it is also much more fragmented and has lower BUSCO scores, suggesting it is incomplete / poorer quality


#### 5.1. Question: 
- **can you think of a reason why variant calling in a haploid sample is less complicated than compared to a diploid (or polyploid) sample?**
    - in a haploid genome, at a given position, all reads will either support a reference base or a variant base. If there is a mix of reference and variant, there is likely to be a technical problem that generally can be filtered
    - in a diploid or polyploid sample, sufficient reads are needed to identify positions that contain both reference and variant bases - a simple fix is having higher sequencing converage. A general rule is 30X for a diploid genome. However, some technical errors are more difficult to distinguish.


#### 8.1. Questions:
- **how many variants were kept after filtering?**
    - should be about 1600 variants

- **it is possible that not all samples will contain all of the SNPs, ie. there is some degree of "missingness". Can you find a flag in the vcftools manual to test this? Are there samples with a lot of missing data?**
    - in this case, there shouldnt be any missingness. It is a bit unusual but it is simply due to having a test dataset
    - normally, some missingness would be expected. Looking at this could help identify samples will really high missingness, which may need to be excluded


#### 9.2. Questions:
- **what proportion of SNPs are in coding regions vs non-coding regions? Why would this happen?**
    - 

- **what proporiton of variants are a "synonymous_variant" and what proportion are "missense_variants"? What effect to these variants have on the coding sequence?**
    - 
    - a missense variant will change the amino acid sequence of a protein coding gene. Also called a non-synonymous variant


- **find a gene with a missense variant - what is the amino acid change, and is it likely to have an effect on the protein? (use the following table to help you: [table](https://en.wikipedia.org/wiki/File:ProteinogenicAminoAcids.svg))**
    - using the table, participants should look to see if the amino acid change results in a different amino acid property 
        - some missense variants will still be within the same functional class, and so perhaps will be less likely to cause an overall functional change to the protein
        - however, a change in functional class may have a greater effect - differnces in charge or hydrophobility, for example, may change the way the protein folds or interacts with other proteins or molecules 

- **can you think of other ways to determine if this variant might impact the function of this protein?**
    - hopefully will link back to the WBPStutorial on using alphafold
    - could perform an experiment, knock-in, knockout and complement, or express in heterologous system


#### 11.1. Questions: 
- **How do these plots compare?**
    - 

- **What is the relative contribution of variance in the PC3/PC4 plot compared to the PC1/PC2 plot?**
    - 


#### 11.2. Questions:
- **Looking at the ellipses specifically, can you see any countries that have a different distribution than the others, and describe this difference?**


#### 12.1. Questions:
- **How do the samples cluster on the tree?**
    - some groups of samples coming from the same country clearly cluster
    - others less so, and are more distributed
    - patterns in both the PCA and tree reflect some deeper ancestral structure in the samples, and also a lot of admixture in the samples due to large scale sheep transport aorund the world in the last few hundred years.

- **Are the PCAs easier or harder to interpret than the tree? Why?**
    - the tree is probably easier to interpret - all the samples are laid out clearly and are labelled
    - however, the tree can also be misleading - if the branch supports are not very high, then some samples might be equally be found on other branches of the tree
    - 