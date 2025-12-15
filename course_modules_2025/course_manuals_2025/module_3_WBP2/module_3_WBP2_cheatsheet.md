
### BLAST exercise <a name="blast_exercise"></a>

Use WormBase ParaSite BLAST to find out the identity of this sequence, and which species it belongs to. Does it have any close hits in other genomes? Try BLASTing against both cDNA and a genomic DNA databases. What kind of sequence is this?

- Use the provided sequence as input in WBPS Blast tool and run two different jobs one against genomic DNA (Search Against: All species, DNA database: Genomic Sequence) and one agaist cDNA (Search Against: All species, DNA database: cDNAs).

- Examine the results:
    - Alignment against Genomic Sequence: First top result - Only 339 from the query's 1453 bases are aligned against the genomic sequence ("Length" column). Click the Sequence button in the "length" column -> only a part of the query is aligned against the genome perfectly (%ID=100%, click the "Alignment" button to examine the alignment)
    - Alignment against cDNA: First top result - The whole query sequence aligns against the cDNA of daf-38 sequence apart from its first 8 bases which are probably part of the UTR that is not annotated. These 8 bases are before the ATG codon. Click the "Sequence" and "Alignment" buttons of the first result to examine this further. 


#### VEP exercise <a name="vep_exercise"></a>

Download the VEP results from the example above as a “VEP file”. Use this file and the original VCF file to answer the following questions:

- Click the "VEP" button located within the facet, which is positioned at the center of the page on the VEP results page.
- Place the downloaded *.vep.txt file inside the module's directory and rename it:
```
mv ~/Downloads/*.vep.txt ~/Module_3_WormBaseParaSite_2/sratti.vep.txt
cd ~/Module_3_WormBaseParaSite_2/
```
1. How many variants were there in the original dataset?
```
# grep to remove the header. wc -l to count the variants of the original vcf file
grep -v "^#" sratti_*.vcf | wc -l| wc -l
```

2. What are the different types of consequence that are found in the file, and how often does each occur?

Lots of ways of doing it...
e.g
```
# remove first line with grep.  extract column 7 with unix cut command, then sort, collapse to unique instances and count
% grep -v '^#' sratti.vep | cut -f 7 | sort | uniq -c
```
Can also be done with awk:
```
# Print a unique list of all the values under the "Consequence" column of the file:
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Consequence") {col=i; break}} NR>1 {print $col}' sratti.vep.txt | sort | uniq

# Count each consequence's occurence and sort them in ascending order
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "Consequence") {col=i; break}} NR>1 {print $col}' sratti.vep.txt | sort | uniq -c | sort
```
3. List all of the variants found in SRAE_2000005500.1.  Which variant or variants show the greatest impact?
```
# One way to list them is:
grep "\tSRAE_2000005500.1\t" sratti.vep.txt

#OR
awk '$5=="SRAE_2000005500.1"' sratti.vep | grep HIGH 

# Another more specific solution by looking for "SRAE_2000005500.1" in the Feature column:
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i == "Feature") {col=i; break} } $col=="SRAE_2000005500.1"' sratti.vep.txt

# The impact of the variants is noted in the IMPACT field of the last column (e.g. IMPACT=MODERATE)
# First list all the possible IMPACT values:
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i == "Feature") {col=i; break} } $col=="SRAE_2000005500.1"' sratti.vep.txt | grep -o 'IMPACT=[^;]*' | sort | uniq

# It looks like the IMPACT=HIGH is the tag for the greatest impact
# List the variants with IMPACT=HIGH for SRAE_2000005500.1
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i == "Feature") {col=i; break} } $col=="SRAE_2000005500.1"' sratti.vep.txt | grep 'IMPACT=HIGH;'
```


4. Create a list of genes where a missense variant is found.
```
# First list all variants that are missense variants (Consequence=missense_variant)
# Simple method:
grep "\tmissense_variant\t" sratti.vep.txt
# More advanced specific method:
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) if ($i == "Consequence") {col=i; break} } $col=="missense_variant"' sratti.vep.txt 

# Then, from the above output get a unique list of all the values of the "Gene" column
grep "\tmissense_variant\t" sratti.vep.txt | cut -f4 | sort | uniq

# Another more advanced way
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) { if ($i == "Consequence") consequenceCol=i; if ($i == "Gene") geneCol=i } } $consequenceCol=="missense_variant" { print $geneCol }' sratti.vep.txt | sort | uniq
```  

5. Find out which genes has the highest number of missense mutations.  View the distribution of variants along the coding sequence in Jbrowse.
```
# Just add a sort | uniq -c in the previous command
awk -F'\t' 'NR==1 { for (i=1; i<=NF; i++) { if ($i == "Consequence") consequenceCol=i; if ($i == "Gene") geneCol=i } } $consequenceCol=="missense_variant" { print $geneCol }' sratti.vep.txt | sort | uniq -c | sort
```
From the previous command we get that the answer is the "WBGene00260238" gene.

To view the VCF in JBrowse:
1. compress and index it.
```
bgzip file.vcf && tabix -p vcf file.vcf.gz
```
2. From the WormBase ParaSite homepage, click either the ”Genome List” tab in the tools bar, or the “Genomes” icon.
3. Scroll down the page to find _Strongyloides ratti_ and click the "Jbrowse" link.
4. Use the Jbrowse search bar and search for "WBGene00260238".
5. From the top menu bar select Track -> Open track file or URL -> Select Files (select your files both .gz.vcf and .gz.vcf.tbi) -> Open


### API exercises



1. List the _Meloidogyne sp._ assemblies by size, smallest to largest.
```
curl -sL 'https://parasite.wormbase.org/rest/info/genomes/taxonomy/Meloidogyne?' -H'Content-type:application/json' | jq -r '.[] | "\(.name)\t\(.assembly.base_count)"' | sort -k 2,2
```
2. Retrieve the protein sequence of the guinea worm transcript DME_0000938001-mRNA-1.
```
curl -Ls 'https://parasite.wormbase.org/rest/sequence/id/DME_0000938001-mRNA-1?' -H 'Content-type:text/plain'
```
3. Write a small program, `get_sequence_for_transcript.sh`, that takes any transcript ID as an argument and returns its protein sequence. For example, running

```
./get_sequence_for_transcript.sh DME_0000938001-mRNA-1
```    
should print:
```
MAKHNAVGIDLGTTYSC...
```
(Hint: shell scripts put arguments from the command line into special variables, named $1, $2 etc )
```
tr_id=$1

curl -Ls "https://parasite.wormbase.org/rest/sequence/id/${tr_id}?" -H 'Content-type:text/plain'
```

4. Retrieve a GFF file of all of the genes located on the AgB01 scaffold of the Ascaris suum PRJNA62057 assembly, between the following coordinates: 5284000 to 5836000.
```
curl -sL 'https://parasite.wormbase.org/rest/overlap/region/ascaris_suum_prjna62057/AgB01:5284000-5836000?feature=gene' -H 'Content-type:text/x-gff3'
```

5. Write a program, `retrieve_genes_in_region.sh` which takes species, scaffold, start and end coordinates as arguments and can return the above for any given region. For example, calling

```
./retrieve_genes_in_region.sh ascaris_suum_prjna62057 AgB01 5284000 5836000
```
should print the same result as question 4.

```
genome=$1
scaffold=$2
start=$3
end=$4

curl -sL "https://parasite.wormbase.org/rest/overlap/region/${genome}/${scaffold}:${start}-${end}?feature=gene" -H 'Content-type:text/x-gff3'
```

### Gene-set enrichment analysis exercise

Use the 24-hour-schistosomule-vs-cercariae.tsv from the previous section and print a list of genes with an adjusted p-value that is less than 0.05, which are most strongly upregulated in the cercariae v the 24h schistosomules.

1. Use gProfiler and perform a Gene-set enrichment analysis for these 40 genes from the "Schistosoma mansoni (PRJEA36577)" organism.
- Extract genes
```
grep -v "^#" 24-hour-schistosomule-vs-cercariae.tsv | grep -v "^gene_id" | awk -F'\t' '$7 != "NA" && $7 < 0.05 && $3 < 0' | sort -g -k 3,3 | head -n 40 | cut -f1
```
-  Paste the list of gene IDs into the central text box. Select "Schistosoma mansoni (PRJEA36577)" under WormBase ParaSite using the "Organism" drop-down menu and then click on "Run Query"

- When results appear, scroll down and hover over the points in the graph to explore gene ontologies which are over-represented in your list of genes. You can also click on "Detailed Results" tab to see a table with all the enriched Gene ontology terms.

2. Which are the 3 most significantly enriched Cellular Component terms? Are they relevant to this developmental stage comparion we're performing?
- When results appear, click on the Detailed Results tab.
- These are the top 3 enriched CC: respirasome, mitochondrion, cytochrome complex

3. Expand the stats by using the ">>" in the header of the GO:CC result table. Try to interpret the T, Q, TnQ and U metrics. What do they represent?
You can read more here: https://biit.cs.ut.ee/gprofiler/page/docs
    * T - Term Size: How many S. mansoni genes are in general associated with this term.
    * Q - Query Size: The number of genes in our gene list (the one we run the analysis for). In our case this number should theoretically be 40, why it's 14?
    * Q - Overlap Size: How many genes from our list are associated with this term.
    * U - Total number of S. mansoni genes.
