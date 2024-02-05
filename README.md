# Compute Gene to Anatomy Derived From Text-Mined Data
Computes associations between genes and anatomies as derived from mentions in biomedical research articles (MEDLINE/PubMed). Note that the genes are likely mixed with proteins.

## Picking Genes
How to use it: From the command line, you can use the go_terms flag to pick a Gene Ontology term of interest, selecting a cellular component, molecular function, or biological process you wish to study. This will automatically curate the genes of interest. (See example command in the Computing Associations section.)

Alternatively, you can also write out the gene list of NCBI/Entrez gene IDs as integers and put it in a "gene_ids.json" file within the input folder. 
Example file:
```
[1, 1324, 142, 790, 341, 927, 817, 814, 42]
```
How it works: PubTator3 has attempted to extract genes from all PubMed articles by a NER + dictionary-based(?) method that normalizes the NER-identified genes to the actual gene IDs. My script here uses their information to find the PubMed article IDs in which the users' genes of interest are found. 

## Picking Anatomy
How to use it: Create a file in the input/ folder. The file should be a 'anatomy_tree_terms.json' which is a dictionary of keys representing the names you provide and values being a list of a list of MeSH tree numbers defining that category.
Example file:
```
{"heart": [["A07.541"]], "brain": [["A08.186.211"]], "eye": [["A09.371"]]}
```

How it works: These MeSH trees are used to query PubMed for the PubMed article IDs of the articles annotated to be studying these topics. Once the anatomy-studying PMIDs are obtained, they are compared with the gene-studying PMIDs above. The intersection between these two lists are used to compute association scores.

## Computing Associations
How to use it: The main script to run is map_gene_to_anatomy.py. This will download what you need (e.g., gets anatomy annotations by calling get_pubmed_docs.py).
Example script (finds all mitochondrial genes in your anatomies defined above):
```
python map_gene_to_anatomy.py --gene_group mitochondrial_genes --go_terms "GO:0005739"
```

How it works: The intersection between the two aforementioned lists are used for the scores. There are three scores: mentions, popularity, and distinctiveness. The mentions are the the number of PMIDs in which a gene was found. The popularity is a parallel metric normalized for the total number of PMIDs in which the genes were found (i.e., gene popularity = # PMIDs of that gene / # PMIDs of all your genes). Distinctivness attempts to measure how specific a gene is to a particular anatomy relative to all other anatomies (i.e., gene distinctiveness = # PMIDs of that gene in your anatomy / # PMIDs of that gene across PubMed) [Note that the score may be slightly off due to the denominator being larger than it should, as some studies mentioning the gene may not be annotated as belonging to an anatomy or that they may actually be unannotated and about your anatomy). 


Feel free to reach out if you have any questions about this. 

Current bugs:
I need to fix the automatic download of protein-to-gene mapping. It doesn't download the correct file sometimes.
