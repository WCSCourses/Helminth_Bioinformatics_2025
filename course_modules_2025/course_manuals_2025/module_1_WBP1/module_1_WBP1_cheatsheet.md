## 1. Genome assembly metrics exercise <a name="genomes_exercise"></a>

1. Find the two other genome assemblies from different _Brugia_ species in WormBase ParaSite, which are of lower quality than _Brugia malayi_.
   1. From the WormBase ParaSite homepage, click either the ”Genome List” tab in the tools bar, or the “Genomes” icon.
   2. Scroll to find the other Brugia genomes: B. pahangi_ and _B. timori_
2. According to their scaffold statistics and BUSCO scores, which of these two assemblies is more contiguous and complete?
   1. Contiguity: Their N50 values are in the last column of the table.
   2. Completeness: The BUSCO scores are summarised in the pie charts. To view them in detail, hover the mouse over each pie chart.
   
   Summary:

    | Genome | N50 | BUSCO Assembly Complete (%) | BUSCO Annotation Complete (%) |
    | ------ | --- | -------------- | -------------------- |
    | B. pahangi | 65,530 | 89.7 | 89.8 |
    | B. timori | 4,903 | 53.3 | 57.3 |
   
   *B. pahangi* is more complete and more contiguous.

---
## 2. Gene page exercise <a name="gene_page_exercise"></a>

Go to the gene page for the _Trichuris muris_ gene TMUE_2000008757 and retrieve the following information:
1. What is the summary description of the gene? Do you have any idea what the gene might be doing from this description?
   1. Go to [WormBase ParaSite](https://parasite.wormbase.org/index.html).
   2. Search for TMUE_2000008757 in the search box at the top right of the page.
   3. On the results page, click the WBGene00286032 link (there should be only 1 search result).
   4. The summary description is located under the Gene Name at the top next to the "Description" header: "Mothers against decapentaplegic homolog \[Source:UniProtKB/TrEMBL;Acc:A0A5S6QPG2\]"

2. How many transcripts of the gene are annotated?
   1. This can be found on the transcript table at the centre of the page.
   2. There is only one transcript.

3. Which strand is the gene on? What is the name of the 5’ neighbouring gene?
   1. Click the ‘Region in Detail’ link in the “Genomic context” section.
   2. The gene is on the reverse (negative) strand - you can see this from the ‘<’ symbol located next to the gene name in the protein coding genes track.
   3. Since the gene is on the reverse (negative) strand 5' would be on the right side: 5' neighbouring gene is TMUE_2000008758.

4. Download the 3’UTR sequence.
   1. Navigate back to the gene page by clicking the "Gene: TMUE_2000008757" tab at the top of the page. 
   2. Click "Sequence" in the navigation menu on the left hand side of the page.
   3. Click "Download sequence" at the center of the page.
   4. Select "Sequences to export" -> "3' UTRs" on the pop-up dialog box.
   5. Click "Download"
5. What identifier would you use to search for the gene in Uniprot?
   1. We could use the Gene Name (WBGene00286032)
   2. We could use its UniProt ID:
      1. On the gene page, click the transcript ID (TMUE_2000008757.1) in the transcript table to navigate to the transcript page.
      2. Click "External references" in the navigation menu on the left hand side of the page.
      3. The UniProt ID is the one next to the "UniProtKB/TrEMBL" external database on the "Extrnal References" table: A0A5S6QPG2.11
   3. You can try both at UniProt's website.
6. Where is this gene’s protein predicted to localise to?
   1. Navigate back to the gene page by clicking the "Gene: TMUE_2000008757" tab at the top of the page.
   2. Click "Cellular component" under "Gene Ontology" in the navigation menu on the left hand side of the page.
   3. According to the Cellular component gene ontologies the protein is predicted to localise to: nucleus, cytoplasm
7. Which Pfam domains is the protein predicted to have? Which of these is responsible for its DNA binding activity?
   1. On the gene page, click the transcript ID (TMUE_2000008757.1) in the transcript table to navigate to the transcript page.
   2. Click the “Domains & features” menu option in the navigation menu on the transcript page.
   3. According to the table there are 2 PFAM domains annotated on the protein: "MAD homology 1, Dwarfin-type" and  "SMAD domain, Dwarfin-type".
   4. You can navigate to the InterPro pages by clicking on the links for these domains in the "Accession" column (PF03165, PF03166). According to the InterPro website pages, the "MAD homology 1, Dwarfin-type" domain is the one responsible for the protein's DNA binding activity.
8. Download the protein alignment of TMUE_2000008757 and its _C. elegans_ orthologue. Is there any published literature on the _C. elegans_ orthologue?
   1. Protein Alignment:
      1. Navigate back to the gene page by clicking the "Gene: TMUE_2000008757" tab at the top of the page.
      2. Click the “Orthologues” menu option under "Comparative Genomics" in the navigation menu on the gene page.
      3. Find the C. elegans orthologue for this gene (hint: you can use the search box at the top right of the table to type "elegans" and filter the results).
      4. Under the "Compare" column of the Orthologues table, click Alignment (protein) to visualise the orthologue alignment between the two proteins.
      5. While on the comparison page, Click "Download Homology" to download the alignment on various formats. 
   2. Literature:
      1. Click the orthologue's gene ID (WBGene00004858) to navigate to the orthologue's gene page.
      2. Click "Literature" in the navigation menu on the left hand side of the C. elegans gene page to find a table with published literature on this gene.
9. Are there any phenotypes associated with this _T. muris_ gene according to the gene page? Which one(s)? Where are these gene-phenotype associations inferred from?
* TMUE_2000008757 gene page -> Click "Phenotypes" in the left navigation menu:
   * Phenotypes associated with this gene: None Found
   * Phenotype, disease and trait annotations associated orthologues of this gene in other species: List of phenotypes associated with the C. elegans ortholog of TMUE_2000008757
* The second table shows phenotype, disease and trait annotations associated orthologues of this gene in other species. Usually, there are not many phenotypic data available for parasitic worms, for this reason exploring phenotypes associated with other species like C. elegans is really useful for scientists. You can use the links in the table to navigate to the orthologous gene in the other species and find more information there. 

---
## 3. BioMart exercise <a name="biomart_exercise"></a>

Use the following _S. ratti_  gene **names** (note: names, not stable IDs) and use BioMart to answer questions 1-5:

Use the list of genes above and generate an output with:
1. their WormBase gene IDs and UniProtKB/TrEMBL IDs. 
2. the InterPro domains that they have been annotated with (InterPro short description). [Q: why do some of the output rows appear multiple times?]
3. the gene stable IDs of their _Strongyloides stercoralis_ orthologues. [Q: which gene has more than one _S. stercoralis_ orthologue?]. 
4. the names of any GO terms associated with the genes.
5. FASTA file of their peptide sequences.

Solution:
For all the above questions, participants should set the same "Query Filters":
- From the WormBase ParaSite homepage, select BioMart from the tool bar, or the BioMart icon. 
- Set the "Query Filters":
   * Select “SPECIES”, tick the “genome” checkbox and scroll down to select “Strongyloides ratti (PRJEB125) [WS285]”. 
   * Select "GENE", tick the "ID list limit" check box, select "Gene Name(s)" on the dialogue box above the text box and copy-paste the given gene list into the text box underneath.

For each question the "Output Attributes" change
- Set the "Output Attributes" by clicking "Output Attributes" on the left menu.
- For question 1:
  - Select "SPECIES AND GENOME INFORMATION": Untick the by-default selected "Genome project". 
  - Select "GENE": Untick the by-default selected "Gene stable ID". 
  - Select "EXTERNAL DATABASE REFERENCES AND ID CONVERSION": Tick "WormBase gene ID" and "UniProtKB/TrEMBL ID" boxes.
  - Click "Results" at top left.
- For question 2:
  - Click "Output Attributes" on the left menu to change the previously set output attributes. Untick all the ones you ticked for the previous question.
  - Select "GENE": Tick the by "Gene stable ID" box.
  - Select "INTERPRO PROTEIN DOMAINS": Tick the "InterPro short description" tick box.
  - Click "Results" at top left.
  - Why do some of the output rows appear multiple times? Some genes will have multiple annotated InterPro domains. These will appear as duplicate rows.
- For question 3:
  - Tick "Output Attributes" on the left menu to change the previously set output attributes. Untick all the ones you ticked for the previous question.
  - Select "GENE": Tick the by "Gene stable ID" box.
  - Select "ORTHOLOGUES": Scroll until you find the "Strongyloides stercoralis (PRJEB528) Orthologues" tab. Then tick the "Strongyloides stercoralis (PRJEB528) gene stable ID" tick box.
  - Click "Results" at top left.
  - which gene has more than one _S. stercoralis_ orthologue? The one with duplicated rows (WBGene00256613 -> SSTP_0001203000 and SSTP_0001203100).
- For question 4:
  - Tick "Output Attributes" on the left menu to change the previously set output attributes. Untick all the ones you ticked for the previous question.
  - Select "GENE ONTOLOGY (GO)": Tick the "GO term name" tick box.
  - Click "Results" at top left.
- For question 5:
  - Tick "Output Attributes" on the left menu to change the previously set output attributes. Untick all the ones you ticked for the previous question.
  - Select the "Retrive sequences" round button at the top of the page.
  - Select "SEQUENCES": Tick the "Peptide" round button.
  - Click "Results" at top left.

---
Use the following _S. mansoni_ gene stable IDs to answer questions 6-9:

* Set the "Query Filters":
   * Select “SPECIES”, tick the “genome” checkbox and scroll down to select “Schistosoma mansoni (PRJEA36577) [WS285]”. 
   * Select "GENE", tick the "ID list limit" check box, select "Gene Name(s)" on the dialogue box above the text box and copy-paste the given gene list into the text box underneath.

6. How many of these genes have orthologues in _S. haematobium_?
* Expand the "Query Filters" selection:
   *  HOMOLOGY -> Restict results to genes with orthologues in... (tick box) -> Schistosoma haematobium (PRJNA78265)
* Click Count at the top left of the page: 16 genes.

7. Generate a table listing the genes in question 6. The table should also has the gene stable ID for the homologue in both species, the homology type (1-1, 1-many, etc), and the % identity between the two orthologues.
* Keep the same "Query Filters"
* Set the "Output Attributes"
   * Untick Genome Project from the SPECIES AND ... tab
   * Tick "Gene stable ID" from the GENE tab
   * Go to the "ORTHOLOGUES" tab and scroll until you find "Schistosoma haematobium (PRJNA78265) Orthologues". Tick "gene stable ID", "Homology type", "%identity"
   * Click Results

8. Of these genes, how many also do not have a human orthologue?
* Just add to the "Query Filters": HOMOLOGY -> Restrict results to genes without orthologues in... (tick box) -> Human
* Click Count at the top left of the page: 11 gene.

9. Retrieve (a) a FASTA file with the CDS sequence of each transcript encoded by the genes from question 6. Make sure that the transcript stable ID is in the header; and (b) a FASTA file containing the CDS sequence plus 100 nt downstream of the stop codon of each of those transcripts. In the header, include the transcript stable ID and the name of the scaffold that the transcript is on.
Use the "Query Filters" from Question 6

* (a)
   * "Output Attributes" -> "Retrieve sequences" at the top of the menu:
   * "SEQUENCES" -> Coding Sequence
   * "HEADER INFORMATION" -> Trascript Attributes -> Transcript stable ID -> Result
* (b)
   * "Output Attributes" -> "Retrieve sequences" at the top of the menu:
   * "SEQUENCES" -> Coding Sequence & "SEQUENCES" -> Downstream flank (tick box): 100
   * "HEADER INFORMATION" -> Trascript Attributes -> Transcript stable ID & Chromosome/scaffold name -> Result

---