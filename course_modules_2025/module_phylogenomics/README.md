# Helminth Bioinformatics
18 - 24 May 2025 

## Topic:  Mitogenomics-Phylogenomics
**Instructor: Dr. Andrés Iriarte**


## Introduction to the Hands-on exercises
**DATASET: Sequences and Table**

Trematodes and cestodes are parasitic flatworms, and their evolutionary relationships are key to understanding how parasitism evolved in Platyhelminthes. Phylogenetic reconstruction is essential not only to trace evolutionary history but also to identify species, study host associations, and explore flatworm biodiversity.

In this practical activity, we will reconstruct phylogenetic relationships based on multiple sequences (multi-gene phylogenetic approach) of selected organisms, mainly trematodes, and identify the main lineages within this group.

For this purpose, we have built a dataset of six mitochondrial protein-coding genes, each consisting of 66 amino acid sequences in fasta format: 50 from trematodes and 16 from cestodes (atp6.faa, cox1.faa, cox2.faa, nad1.faa, nad2.faa, and nad3.faa). 

A file containing the four-letter code, complete name, and class (trematode/cestode) of each species is provided and required for the analysis (Codes-Names.tab). 

## Activity I: Multiple Sequence Alignment (MSA)

**Specific objectives of this practice:**

- To become familiar with the MAFFT and visualizing alignments.
- To obtain an adequate sequence alignment for later use in phylogeny.

To carry out an alignment, the sequences need to be available in a file that can be read by the programs. In general, the FASTA format is accepted by most sequence alignment and edition programs. The time required for the analysis will depend on the computer's processing power as well as the length and number of sequences to be analyzed. In general, it can be estimated that the computation time will increase linearly with the length of the sequences, and exponentially with the number of sequences to be aligned.

The dataset (unaligned set of sequences: **.faa files** and table) is in:

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

### A.	Alignment with MAFFT

MAFFT is an advanced tool that can align using different alignment algorithms for different applications such as L-INS-i (accurate; recommended for <200 sequences), FFT-NS-2 (fast; recommended for >2,000 sequences), etc. It can be run locally or on online servers. To understand the algorithms and their use cases, please refer to https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html

To use it on the VM, type mafft on the command-line, mafft --help will give you information about the proper syntax. 

>*Please note that the procedure below is for all six proteins. As you need to align the sequences of each protein independently.*


1. Open a Terminal and go into the directory	that contains the dataset to align.

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

2. Type: 

```
mafft --auto atp6.faa > atp6.alg.faa
mafft --auto cox1.faa > cox1.alg.faa
mafft --auto cox2.faa > cox2.alg.faa
mafft --auto nad1.faa > nad1.alg.faa
mafft --auto nad2.faa > nad2.alg.faa
mafft --auto nad3.faa > nad3.alg.faa
```

>**Usage:**
>mafft [options] input > output
>
>--auto: automatically switches algorithms according to data size.

-------------------------
[OPTIONAL] Visualizing the alignment online in https://alignmentviewer.org/).
Open your preferred web browser and upload a file containing one of your recent sequence alignments.

![alignmentviewer](https://github.com/user-attachments/assets/6aa9ed93-c215-45bd-855f-78632c8fb706)


-----------------

### B.	Trimming the alignment with TrimAl

TrimAl is a command-line tool used to automatically trim multiple sequence alignments by removing poorly aligned or overly divergent regions. These regions can negatively affect downstream analyses such as phylogenetic reconstruction. TrimAl supports various trimming modes, from conservative to aggressive, and is particularly useful for preparing alignments for accurate phylogenetic inference. To understand options and their use cases, please refer to https://trimal.readthedocs.io/en/latest/.

To use it on the VM, type trimal on the command-line, trimal -h will give you information about the proper syntax. 

>*"Again, please note that the procedure below is for all six proteins. As you need to trim the alignment of each protein independently.*

1. Open a Terminal and go into the directory	that contains the alignments.

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

2. Type: 

```
trimal -in atp6.alg.faa -out atp6.alg.trim.faa -gt 0.3 -cons 0.7
trimal -in cox1.alg.faa -out cox1.alg.trim.faa -gt 0.3 -cons 0.7
trimal -in cox2.alg.faa -out cox2.alg.trim.faa -gt 0.3 -cons 0.7
trimal -in nad1.alg.faa -out nad1.alg.trim.faa -gt 0.3 -cons 0.7
trimal -in nad2.alg.faa -out nad2.alg.trim.faa -gt 0.3 -cons 0.7
trimal -in nad3.alg.faa -out nad3.alg.trim.faa -gt 0.3 -cons 0.7
```

The option -gt 0.3 removes columns (alignment positions) with gaps in more than 30% of the sequences, while -cons 0.7 retains only columns where the residue conservation score is at least 0.7. 

### C.	Convert the alignment to Clustal format to facilitate visualization using the **seqret** program from the EMBOSS package.

seqret is a command-line tool from the EMBOSS (European Molecular Biology Open Software Suite) package used to convert between different sequence formats. For example, it can convert a FASTA file to Clustal, GenBank, Phylip, and many other formats. It works with both DNA and protein sequences, and is often used to prepare files for input into other bioinformatics tools. Please refer to https://www.bioinformatics.nl/cgi-bin/emboss/help/seqret for further details and options.

Clustal format is ideal for viewing alignments because it is human-readable. Sequences are displayed in blocks with aligned residues stacked vertically, making them easy to compare. 

To use it on the VM, type man seqret on the command line to see the manual information. Press letter “q” to exit the manual page.

```
man seqret
```

1. Open a Terminal and go into the directory	that contains the trimmed files.

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

2. Type: 

```
seqret atp6.alg.trim.faa -osformat clustal atp6.alg.trim.aln
seqret cox1.alg.trim.faa -osformat clustal cox1.alg.trim.aln
seqret cox2.alg.trim.faa -osformat clustal cox2.alg.trim.aln
seqret nad1.alg.trim.faa -osformat clustal nad1.alg.trim.aln
seqret nad2.alg.trim.faa -osformat clustal nad2.alg.trim.aln
seqret nad3.alg.trim.faa -osformat clustal nad3.alg.trim.aln
```
3. Use *less* command to visualize the clustal files generated and press letter “q” to exit the less command.

```
less atp6.alg.trim.aln
```
>*Do you think this is a conserved protein sequence?*

### D.	Concatenate the alignments of the six alignment and trimmed sequences using a script named *Concatenator.scp*

Concatenator.scp  is a locally developed script written in bash that we use in our lab to concatenate sequences from different alignments to build a “super protein” or “super gene”. This script takes a *list of organism codes and a list of alignment files*. For each organism, it extracts matching sequences from all input alignment files, appends them in order, and produces a final concatenated FASTA file (concatenated.fas).


1. Open a Terminal and go into the directory	that contains the trimmed files.

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

2. To generate a file with the list of alignment files to concatenate type: 

```
ls *.alg.trim.faa > alignment_to_concatenate.list
```

3. To extract codes (first column) from the table of codes and names and generate a file type:

```
awk '{print $1}' Codes-Names.tab > codes.list
```
4. Run Concatenator.scp:

```
chmod +x Concatenator.scp
./Concatenator.scp alignment_to_concatenate.list codes.list
```
>*Remember to provide the first and second arguments to the script: a file containing the list of alignment files, and a file containing the list of codes (one per species), respectively.*

5. Convert the alignment to Clustal format to facilitate visualization using the **seqret** tool.

```
seqret concatenated.fas -osformat clustal concatenated.aln
less concatenated.aln
```
-------------------------
[OPTIONAL] Get some info from alignment using infoalign tool from emboss package. 

The output file from infoalign provides summary statistics for each sequence in a multiple sequence alignment. It includes the sequence name, its aligned length, the number and percentage of residues identical to the consensus, the number of differences, and the number of insertions or deletions (indels). Depending on the options used, it may also include a similarity score and a consensus sequence. The level of detail in the output can be customized using various command-line options.


```
infoalign concatenated.fas concatenated.infoalign
less concatenated.infoalign
```
-------------------------

### E.	Add the species name to the headers using the script named *Replacetator.scp*.

This script requires three input files: (i) the FASTA file whose headers will be modified, (ii) a table with at least two columns—the first containing the codes and the second the full species names, and (iii) the name you choose for the output file (in this case, “concatenated.rn.fas”).

1. Run Replacetator.scp:

```
chmod +x Replacetator.scp
./Replacetator.scp concatenated.fas Codes-Names.tab concatenated.rn.fas
```
2. Use *less* command to visualize the generated file and press letter “q” to exit:

```
less atp6.alg.trim.aln
```

## Activity 2: Phylogenetic analysis

Specific objectives of this practice:
- To become familiar with the IQ-TREE and Figtree programs.
- To build a Maximum Likelihood tree to identify main lineages within the Trematoda class.

As the input file, we will use the FASTA file generated in the previous activity: concatenated.rn.fas.

**Introduction to the IQ-TREE program:**

This program allows you to perform phylogenetic analysis by Maximum Likelihood. It uses highly efficient algorithms to explore the tree space, allowing very large matrices to be analyzed with reliable results (hundreds or thousands of sequences). It allows estimating the evolutionary model (ModelFinder module) followed by the phylogenetic inference and implements support measures to evaluate the reliability of the groupings or branches (Bootstrap, Ultrafast Bootstrap Approximation and probabilistic contrasts). The program can be downloaded and run locally (http://www.iqtree.org/), or on online servers such as http://iqtree.cibiv.univie.ac.at/.
You can find many basic and advanced tutorials at http://www.iqtree.org/doc/

>*Please note that the following procedure applies to the six concatenated protein sequences of mitochondrial-encoded genes. If you are analyzing a different dataset, be sure to replace the file names in the instructions and adjust the parameters accordingly.*

### A.	Phylogenetic Inference by Maximum Likelihood with IQ-TREE

**Phylogenetic inference + support (Ultrafast Bootstrap Approximation + SH-aLRT)**

1.	Open a Terminal and go into the folder that contains the concatenated alignment file previously generated:

```
/home/manager/course_github/course_modules_2025/module_phylogenomics/Data/
```

2.	Type: 
```
iqtree2 -h 
```
>*This command allows you to see all available options, check those that you will use in the next step.*

3.	Type: 
```
iqtree2 -s concatenated.rn.fas -m mtZOA+F+I+R6 -B 1000 -nt 6
```

>**Usage:**

>**-s:** Specifies the input sequence alignment file.

>**-m:** specify a model selection strategy (if no option is specified, -m MFP is used by default). mtZOA+F+I+R6: This model is designed for mitochondrial protein-coding genes from invertebrates (e.g., cnidarians, arthropods) and captures their unique substitution patterns, including empirical amino acid frequencies (+F), invariant sites (+I), and rate variation across sites with 6 categories (+R6). You can use “-m MFP” if you want iqtree2 to run automatic model selection using the Model Finder Plus (MFP) approach.

>**-B:** specify the number of replicates for Ultrafast Bootstrap Approximation in IQ-TREE

>*(OPTIONAL *-alrt:* Performs the number of SH-aLRT (Shimodaira-Hasegawa-like approximate likelihood ratio test) replicates for additional support.)*

>**-nt:** specify the number of threads for the analysis (you can use “-nt AUTO” if you want iqtree2 to automatically detect the number of CPU cores available to parallelize the analysis).


Once the process is finished, the output files will be found in the folder, including:
**.treefile:** the ML tree in NEWICK format, which can be visualized by any supported tree viewer programs like FigTree. **.iqtree:** the main report file that is self-readable. You should look at this file to see the computational results. It also contains a textual representation of the final tree. **.log:** log file of the entire run (also printed on the screen).

### B.	Tree visualization

1. Open a Terminal and type: 
```
figtree 
```

2. Open the ML tree: 
```
File -> Open -> select the file concatenated.rn.fas.treefile
```

3. Select a name for annotated values: **“UFB”**

4. In In this exercise, root the tree at the branch that separates cestodes from trematodes:
```
Select the branch -> Reroot
```

5. It is recommended to order the tree to improve its visualization:
```
Tree -> Increasing Node Order
```

6. Display support values: 
```
Branch labels -> Display "UFB"
```

7. Modify tree presentation according to your preference:

```
e.g.  can modify the size of the fonts (in Tip Labels, Legend, etc).
```

8. Check whether the taxonomy aligns with the phylogenetic pattern we observed:

### Trematode Taxonomy of Species Included in the Phylogenetic Analysis

| **Order**         | **Family**           | **Species**                                                                                                                                 |
|-------------------|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| **Plagiorchiida** | Allocreadiidae       | *Allocreadium lobatum*                                                                                                                        |
|                   | Azygiidae            | *Azygia hwangtsiyui*, *Azygia longa*, *Azygia susquehannae*                                                                                         |
|                   | Brachycladiidae      | *Brachycladium goliath*                                                                                                                       |
|                   | Bucephalidae         | *Carassotrema koreanum*                                                                                                                        |
|                   | Clinostomidae        | *Clinostomum complanatum*, *Clinostomum piscidium*                                                                                               |
|                   | Opisthorchiidae      | *Clonorchis sinensis*, *Metorchis bilis*, *Metorchis orientalis*, *Metorchis xanthosomus*, *Opisthorchis sudarikovi*, *Opisthorchis viverrini*           |
|                   | Cryptogonimidae      | *Creptotrematina aguirrepequenoi*                                                                                                             |
|                   | Heterophyidae        | *Cryptocotyle lingua*                                                                                                                         |
|                   | Cyathocotylidae      | *Cyathocotyle prussica*                                                                                                                       |
|                   | Echinostomatidae     | *Echinoparyphium aconiatum*, *Echinostoma caproni*, *Echinostoma hortense*, *Echinostoma miyagawai*, *Echinostoma paraensei*, *Echinostoma revolutum* |
|                   | Fasciolidae          | *Fasciola hepatica*, *Fasciola jacksoni*, *Fascioloides magna*, *Fasciolopsis buskii*                                                               |
|                   | Morishitiidae        | *Morishitium polonicum*     
|                   | Diplodiscidae        | *Diplodiscus japonicus*, *Diplodiscus mehrai*, *Diplodiscus nigromaculati*  |
|                   | Notocotylidae        | *Ogmocotyle ailuri*                                                                                                                           |
|                   | Paragonimidae        | *Paragonimus skrjabini miyazakii*, *Paragonimus westermani*                                                                                     |
|                   | Plagiorchiidae       | *Plagiorchis elegans*, *Proterometra macrostoma*                                                                                                 |
|                   | Allocreadiidae       | *Pseudoparacreptotrema yaguezani*                                                                                                             |
| **Diplostomida**  | Schistosomatidae     | *Schistosoma bovis*, *Schistosoma curassoni*, *Schistosoma guineensis*, *Schistosoma haematobium*, *Schistosoma indicum*, *Schistosoma japonicum*, *Schistosoma mansoni*, *Schistosoma margrebowiei*, *Schistosoma mattheei*, *Schistosoma mekongi*, *Schistosoma spindalis*, *Trichobilharzia regenti*, *Trichobilharzia szidati* |

-------------------------
[OPTIONAL] Check Olson PD, Cribb TH, Tkach VV, Bray RA, Littlewood DTJ. (2003). *Phylogeny and classification of the Digenea (Platyhelminthes: Trematoda).* International Journal for Parasitology, 33(7), 733–755. for furhter

