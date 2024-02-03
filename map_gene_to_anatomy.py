import argparse
import os
import subprocess
import json
from bisect import bisect_left
import pickle
from sys import exit
from collections import defaultdict

class AssociateGeneToAnatomy():
        
    def __init__(self, download=True, gene_to_pmid_filename='input/gene2pubtator3', genes=[]):
        self.gene_to_pmid_filename = gene_to_pmid_filename
        self.gene_id_to_gene_name = {}
        self.gene_ids_of_interest = sorted(genes)
        self.gene_to_anatomy_freq = {}
        self.gene_to_freq_in_same_anatomy = {}
        self.gene_uniqueness_to_anatomy = {}
        self.pmid_to_genes_filt = {}
        self.gene_to_pmids_filt = {}
        self.topic = ''
        
        # Download mapping file
        self.download = download
        if self.download:
            self.download_gene_annotations()
        
        # Rows in gene2pubtator
        self.number_of_gene_entries = self.get_number_of_gene_entries()
        
        
    def download_gene_annotations(self):
        '''
        Download PubTator3's file mapping genes to PubMed IDs
        '''
        pubtator3_download = f'https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/{self.gene_to_pmid_filename}.gz'
        os.system(f'wget -NP input/ {pubtator3_download}')
        try:
            os.system(f'gunzip {self.gene_to_pmid_filename}.gz')
            os.system(f'mv {self.gene_to_pmid_filename} input/{self.gene_to_pmid_filename}')
        except:
            pass
    
    def get_number_of_gene_entries(self):
        '''
        Get the total number of rows in the PubTator3 gene file
        '''
        result = subprocess.run(['wc', self.gene_to_pmid_filename, '-l'], 
                                stdout=subprocess.PIPE, 
                                text=True)
        output = result.stdout.strip()
        number_of_gene_entries = int(output.split()[0])
        
        return number_of_gene_entries
        
        
    def extract_gene_to_pmids(self):
        '''
        Based on the PubTator3 file, map the genes to the PubMed IDs of the
        articles they were mentioned in.
        '''
        self.pmid_to_genes = defaultdict(list)
        self.gene_to_pmids = defaultdict(list)

        with open(self.gene_to_pmid_filename) as fin:
            for idx, line in enumerate(fin):
                pmid, _, gene_ids, gene_names, _ = map(str.strip, line.split('\t'))
                gene_ids = map(int, gene_ids.split(';'))
                pmid = int(pmid)
                for gene_id in gene_ids:
                    if gene_names != '':
                        self.gene_id_to_gene_name[gene_id] = gene_names
                    if binary_search_bool(self.gene_ids_of_interest, gene_id):
                    #if gene_id in self.gene_ids_of_interest:
                        self.pmid_to_genes[pmid].append(gene_id)
                        self.gene_to_pmids[gene_id].append(pmid)
                if idx%10000 == 0:
                    print(f'Parsing PubTator3 Table: {idx}/{self.number_of_gene_entries}', end='\r')
        print('\nDone mapping gene to PMID')
        

    def set_anatomy_pmids(self, anatomy):
        pmid_to_category_path = f'output/{anatomy}/pmid_to_category_{anatomy}.json'
        if not os.path.exists(pmid_to_category_path):
            print(f'Downloading the PMIDs for {anatomy}')
            os.system(f'python get_pubmed_docs.py --topic {anatomy} --get_pmids_via_mesh')
        anatomy_pmids = list(json.load(open(pmid_to_category_path)).keys())
        self.anatomy_pmids = sorted([int(pmid) for pmid in anatomy_pmids])
        self.topic = anatomy
        
        
    def set_pmid_to_genes(self, pmid_to_genes):
        '''
        In case you wanted to manually override the pmid-to-genes 
        '''
        self.pmid_to_genes = pmid_to_genes
        
        
    def set_gene_to_pmids(self, genes_to_pmid):
        '''
        In case you wanted to manually override the pmid-to-genes 
        '''
        self.gene_to_pmids = genes_to_pmid
    
    
    def filter_for_gene_anatomy_pmids(self):
        '''
        Start with the PMID-to-genes mapping. Remove all entries whose PMIDs do not
        also study your anatomy of interest.
        '''
        ##print('len(self.pmid_to_genes))', len(self.pmid_to_genes))
        pmid_to_genes_filt = {pmid:genes for pmid, genes in self.pmid_to_genes.items() 
                              if binary_search_bool(self.anatomy_pmids, pmid)}        
        self.pmid_to_genes_filt = pmid_to_genes_filt
        self.gene_to_pmids_filt = invert_a_dict_with_list_values(pmid_to_genes_filt)        
        ##print('len(self.gene_to_pmids_filt)', len(self.gene_to_pmids_filt))
        ##print('len(self.pmid_to_genes_filt)', len(self.pmid_to_genes_filt))
        
        
    def calculate_frequency(self):
        # Gene mentions in the anatomy of interest (Anatomy-Specific Distribution)
        '''
        Calculate how frequent a gene is studied/mentioned in "your anatomy" studies.
        Not normalized.
        '''
        gene_to_anatomy_freq = {}
        for gene, anatomy_pmids in self.gene_to_pmids_filt.items():
            gene_to_anatomy_freq[gene] = len(anatomy_pmids)
            
        anat_genes_not_found = set(self.gene_ids_of_interest).difference(gene_to_anatomy_freq.keys())
        for gene in anat_genes_not_found:
            gene_to_anatomy_freq[gene] = 0    
        self.gene_to_anatomy_freq = sort_dict_by_value(gene_to_anatomy_freq)
    
    
    def calculate_frequency_relative_to_anatomy(self):
        # Gene relative to all genes of interest in the same anatomy of interest
        '''
        Calculate how frequent a gene is studied in a "your anatomy" study 
        relative to the other genes of interest.
        Frequency of mentions in "your anatomy" studies.
        '''
        total_mentions = 0
        for gene, anatomy_pmids in self.gene_to_pmids_filt.items():
            total_mentions += len(anatomy_pmids)

        gene_to_freq_in_same_anatomy = {}
        for gene, anatomy_pmids in self.gene_to_pmids_filt.items():
            mentions = len(anatomy_pmids)
            relative_mentions = mentions/total_mentions
            gene_to_freq_in_same_anatomy[gene] = relative_mentions
        self.gene_to_freq_in_same_anatomy = sort_dict_by_value(gene_to_freq_in_same_anatomy)

        
    def calculate_frequency_relative_to_all_pubmed(self):
        # Gene relative to all anatomies
        '''
        Calculate how unique a gene is to "your anatomy" studies.
        Relative frequency in of mentions in "your anatomy" studies 
        over all studies. Counts mentions across studies, not studies.
        '''
        gene_to_freq_in_all_studies = {}
        for gene, pmids in self.gene_to_pmids.items():
            gene_to_freq_in_all_studies[gene] = len(pmids)

        gene_uniqueness_to_anatomy = {}
        for gene, all_freq in gene_to_freq_in_all_studies.items():
            anatomy_freq = self.gene_to_anatomy_freq[gene]
            unique_score = anatomy_freq/all_freq
            gene_uniqueness_to_anatomy[gene] = unique_score
        self.gene_uniqueness_to_anatomy = sort_dict_by_value(gene_uniqueness_to_anatomy)
    
    
    def calculate_scores(self):
        '''
        Calculate all the scores of genes in anatomy.
        Uses the genes and anatomy you chose. 
        '''
        self.calculate_frequency()
        self.calculate_frequency_relative_to_anatomy()
        self.calculate_frequency_relative_to_all_pubmed()
        
        
    def display_score(self):
        '''
        Prints each gene-to-anatomy score
        '''
        # Shows the most popular gene in your anatomy (# mentions)
        print(self.gene_to_anatomy_freq, '\n')

        # Shows the most popular gene in your anatomy (relative to other genes in the anatomy)
        print(self.gene_to_freq_in_same_anatomy, '\n')

        # Shows the genes most unique to / enriched in your anatomy 
        # 0 indicates a gene is never mentioned related to your anatomy
        # 1 indicates a gene is only mentined related to studies of your anatomy, not another. 
        # Note that 1 could still occur for e.g., heart and be non-zero for e.g., brain if there
        # is a study that is of heart and brain
        print(self.gene_uniqueness_to_anatomy, '\n')

        
    def switch_gene_id_keys_to_gene_names(self, gene_to_score):
        gene_name_to_score = {}
        for gene_id, score in gene_to_score.items():
            try:
                gene_name = self.gene_id_to_gene_name[gene_id]
            except:
                gene_name = gene_id
            gene_name_to_score[gene_name] = score
        return gene_name_to_score
    
    
    def save_scores(self):
        '''
        Export the score dictionaries
        '''
        topic = self.topic
        num_genes = len(self.gene_ids_of_interest)
            
        # Popularity (Gene ID)
        outfile = f'output/{topic}/mentions_{num_genes}_gene_ids.json'
        with open(outfile,'w') as fout:
            json.dump(self.gene_to_anatomy_freq, fout)

        # Popularity (Gene Name)
        gene_name_to_anatomy_freq = self.switch_gene_id_keys_to_gene_names(self.gene_to_anatomy_freq)
        outfile = f'output/{topic}/mentions_{num_genes}_gene_names.json'
        with open(outfile,'w') as fout:
            json.dump(gene_name_to_anatomy_freq, fout)
            
        # Popularity Relative to Other Genes in the anatomy (Gene ID)
        outfile = f'output/{topic}/popularity_{topic}_{num_genes}_gene_ids.json'
        with open(outfile,'w') as fout:
            json.dump(self.gene_to_freq_in_same_anatomy, fout)
            
        # Popularity Relative to Other Genes in the anatomy (Gene Name)
        gene_name_to_freq_in_same_anatomy = self.switch_gene_id_keys_to_gene_names(self.gene_to_freq_in_same_anatomy)
        outfile = f'output/{topic}/popularity_{topic}_{num_genes}_gene_names.json'
        with open(outfile,'w') as fout:
            json.dump(gene_name_to_freq_in_same_anatomy, fout)
            
        # Popularity Relative to Other Anatomies (Gene ID)
        outfile = f'output/{topic}/distinctiveness_{topic}_to_anatomies_{num_genes}_genes.json'
        with open(outfile,'w') as fout:
            json.dump(self.gene_uniqueness_to_anatomy, fout)
            
         # Popularity Relative to Other Anatomies (Gene Name)
        gene_name_uniqueness_to_anatomy = self.switch_gene_id_keys_to_gene_names(self.gene_uniqueness_to_anatomy)
        outfile = f'output/{topic}/distinctiveness_{topic}_to_anatomies_{num_genes}_gene_ids.json'
        with open(outfile,'w') as fout:
            json.dump(gene_name_uniqueness_to_anatomy, fout)    
        
        # Gene ID to Name 
        with open(f'output/{topic}/gene_id_to_names_{num_genes}_gene.json','w') as fout:
            json.dump(self.gene_id_to_gene_name, fout)
        
        
        
        
def initialize_gene_to_anatomy_class(gene_group, gene_ids):
    class_path = f'input/GeneToAnatomyClass_{gene_group}.pkl'
    if os.path.exists(class_path): 
        with open(class_path, 'rb') as fin:
            GTA = pickle.load(fin)
    else:
        GTA = AssociateGeneToAnatomy(genes=gene_ids)
        GTA.extract_gene_to_pmids()
        with open(class_path, 'wb') as fout: 
            pickle.dump(GTA, fout)

    return GTA

        
def binary_search_bool(sorted_list, target):
    '''
    Uses binary search to look for an integer item,
    "targeT" in a sorted list
    '''
    index = bisect_left(sorted_list, target)
    if index != len(sorted_list) and sorted_list[index] == target:
        return True  # Found the target
    return False  # Target not found


def invert_a_dict_with_list_values(the_dict):
    '''
    Used for taking the filtered pmid-to-genes dictionary
    and making a filtered gene-to-pmids dictionary
    '''
    inverted_dict = {}
    for key, values in the_dict.items():
        for value in values:
            inverted_dict.setdefault(value, []).append(key)
    return inverted_dict


def sort_dict_by_value(my_dict):
    sorted_dict_desc = dict(sorted(my_dict.items(), key=lambda item: item[1], reverse=True))
    return sorted_dict_desc


def prepare_anatomy_file(anatomy, tree_numbers):
    anatomy_tree_numbers_path = f'input/{anatomy}_tree_numbers.json'
    if not os.path.exists(anatomy_tree_numbers_path):
        with open(anatomy_tree_numbers_path, 'w') as fout:
            json.dump(tree_numbers, fout)
            

def switch_dict_set_to_dict_list(dict_set):
    dict_set_copy = dict_set.copy()
    for k,v in dict_set_copy.items():
        dict_set[k] = list(v)
            
    
'''Gene Ontology Functions'''
def download_go_to_protein():
    if not os.path.exists('input/goa_human.gaf'):
        go_gaf_url = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
        os.system(f'wget -N -P input/ {go_gaf_url}')
        os.system(f'gunzip input/goa_human.gaf.gz')
        

def map_protein_to_go():
    protein_to_go = defaultdict(set)
    go_to_protein = defaultdict(set)

    with open('input/goa_human.gaf') as fin:
        for i, line in enumerate(fin):
            if line.startswith('!'):
                continue
            line = line.split('\t')
            protein_id = line[1]
            go_term = line[4]
            protein_to_go[protein_id].add(go_term)
            go_to_protein[go_term].add(protein_id)
    
    switch_dict_set_to_dict_list(go_to_protein)
    
    with open('input/go_to_protein.json','w') as fout:
        json.dump(go_to_protein, fout)


def download_gene_to_protein():
    gene_to_protein_file = 'input/gene_to_protein_hgnc.tsv'
    if not os.path.exists(gene_to_protein_file):
        base_url = 'https://www.genenames.org/cgi-bin/download/'
        download_link =  "custom?col=md_eg_id&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=md_prot_id&format=text&submit=submit"
        os.system(f'wget -NP input/ {base_url+download_link}')
        os.system(f'mv input/custom?col=md_eg_id {gene_to_protein_file}')
        #os.system(f'mv input/{download_link} {gene_to_protein_file}')

    
def map_gene_to_protein():
    protein_to_gene = defaultdict(set)

    with open('input/gene_to_protein_hgnc.tsv') as fin:
        for idx, line in enumerate(fin):
            if idx == 0:
                continue
            try:
                gene, protein = [item.strip() for item in line.split('\t')]
                if gene != '' and protein != '':
                    gene = int(gene)
                    protein_to_gene[protein].add(gene)
            except:
                pass
            
    switch_dict_set_to_dict_list(protein_to_gene)
    with open('input/protein_to_gene.json','w') as fout:
        json.dump(protein_to_gene, fout)
    
    
def map_go_to_genes(go_terms):
    go_terms = go_terms.split(',')
    go_to_protein = json.load(open('input/go_to_protein.json'))
    protein_to_gene = json.load(open('input/protein_to_gene.json'))
    ##print(protein_to_gene)
    gene_list = []
    ##print(go_terms)
    y, n = 0,0
    for go_term in go_terms:
        proteins = go_to_protein[go_term]
        print(len(proteins), 'proteins')
        for protein in proteins:
            try:
                genes = protein_to_gene[protein]
                y += 1
            except:
                n += 1
                continue
            gene_list.extend(genes)
        ##print('genes', len(gene_list))
        
    ##print(y, n)
    gene_list = list(set(gene_list))
    ##print(gene_list)
    
    return gene_list
        
    
def get_genes_for_go_terms(go_terms):
    # Get mapping files
    if True:
    #if not os.path.exists('input/protein_to_gene.json') or not os.path.exists('input/go_to_protein.json'):
        download_go_to_protein()
        map_protein_to_go()
        download_gene_to_protein()
        map_gene_to_protein()
    
    print(go_terms)
    gene_ids = map_go_to_genes(go_terms)
    
    return gene_ids
            

    
    
    
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Map genes to anatomy')
    parser.add_argument('--anatomy_file', '-a', type=str, default='input/anatomy_tree_numbers.json',
                        help='File name of dictionary mapping the name of anatomies to their MeSH numbers')
    parser.add_argument('--gene_id_file', '-g', type=str, default='input/gene_ids.json',
                        help='File name of the gene lists')
    parser.add_argument('--gene_ids', type=str, default='',
                        help='Comma-separated list of NCBI gene IDs')
    parser.add_argument('--gene_group', type=str, default='',
                        help='Name of the gene group')
    parser.add_argument('--go_term_file', '-gf', type=str, default='',
                        help='File name of a list of GO terms. Will get the associated genes')
    parser.add_argument('--go_terms', '-go', type=str, default='',
                        help='Comma-separated list of GO terms')
    
    args = parser.parse_args()
    anatomies_file = args.anatomy_file
    gene_id_file = args.gene_id_file
    gene_ids = args.gene_ids
    go_term_file = args.go_term_file
    go_terms = args.go_terms
    gene_group = args.gene_group
    
    if os.exists.path('input'):
        os.mkdir('input')
    if os.exists.path('output'):
        os.mkdir('output')
        
        
    ''' Map genes to PMIDs '''
    if gene_id_file != '':
        gene_ids = json.load(open(gene_id_file))
        
    if go_terms != '':
        gene_ids = get_genes_for_go_terms(go_terms)
        
    if gene_group == '':
        gene_group = f'{len(gene_ids)}_genes'
        
    GTA = initialize_gene_to_anatomy_class(gene_group, gene_ids)
    
        
    ''' Map genes to anatomy '''
    anatomies = json.load(open(anatomies_file)) 
    for anatomy, tree_numbers in anatomies.items():
        print(f'\nScoring {anatomy}.')
        
        ''' Map Anatomy to PMID '''
        prepare_anatomy_file(anatomy, tree_numbers)
        GTA.set_anatomy_pmids(anatomy)

        ''' Compute associations between the Genes and Anatomy '''
        GTA.filter_for_gene_anatomy_pmids()
        GTA.calculate_scores()
        #GTA.display_score()
        GTA.save_scores()
    


        
# change gene name file to not include the big one

# provide the PMIDs
# do this by sving the filtere gene_to_

# provide the text 

# automate gene-protein mapping