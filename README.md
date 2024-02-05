# Compute Gene to Anatomy Derived From Text-Mined Data
Computes associations between genes and anatomies as derived from mentions in biomedical research articles (MEDLINE/PubMed). 

## Picking Genes
From the command line, you can use the go_terms flag to pick a Gene Ontology term of interest, selecting a cellular component, molecular function, or biological process you wish to study. This will automatically curate the genes of interest. You can also write out the gene list of NCBI/Entrez gene IDs as integers and put it in a "gene_ids.json" file within the input folder. 
Example file:
```
[1, 1324, 142, 790, 341, 927, 817, 814, 42]
```

## Picking Anatomy
Create a file in the input/ folder. The file should be a 'anatomy_tree_terms.json' which is a dictionary of keys representing the names you provide and values being a list of a list of MeSH tree numbers defining that category.
Example file:
```
{"heart": [["A07.541"]], "brain": [["A08.186.211"]], "eye": [["A09.371"]]}
```

## Computing Associations
The main script to run is map_gene_to_anatomy.py. This will download what you need (e.g., gets anatomy annotations by calling get_pubmed_docs.py).
Example script (finds all mitochondrial genes in your anatomies defined above):
```
python map_gene_to_anatomy.py --gene_group mitochondrial_genes --go_terms "GO:0005739"
```



Feel free to reach out if you have any questions about this. 

Current bugs:
I need to fix the automatic download of protein-to-gene mapping. It doesn't download the correct file sometimes.
