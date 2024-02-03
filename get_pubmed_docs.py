from datetime import datetime
import time
import urllib.request
import xml.etree.ElementTree as ET
import json
import os
from itertools import chain
try:
    from pubmed_api.pubmed import PubMedAPI
except:
    os.system('pip3 install pubmed-api')
    from pubmed_api.pubmed import PubMedAPI
import requests
import re
import html
import codecs
import csv
import argparse

'''
Generic utility function
'''
def switch_dictset_to_dictlist(the_dict):
    '''
    FUNCTION:
    - Make a new dictionary with values as lists 
      instead of values as sets
      
    PARAMS:
    - the_dict: The initial dict with values of sets
    '''
    dictlist = dict()
    for k in the_dict.copy():
        dictlist[k] = list(the_dict[k])
        
    return dictlist


'''
Pick MeSH IDs
'''
def download_mesh_xml(year=-1):
    print('Downloading MeSH XML...')
    if year == -1:
        year = datetime.now().year
    url = f'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc{year}.xml'
    dest = f'input/desc{year}.xml'
    if os.path.exists(dest):
        print('Already downloaded')
    else:
        urllib.request.urlretrieve(url, dest);
        print('Download complete!')

def parse_mesh_xml(year=-1):
    print('Parsing MeSH XML...')
    if year == -1:
        year = datetime.now().year 
    tree = ET.parse(f'input/desc{year}.xml')
    root = tree.getroot()   
    print('Parsing complete!')
    return root


def align_mesh_trees_with_terms(root, mesh_roots = ('A','C','F03','G')):
    print('Aligning MeSH tree numbers with MeSH terms...')
    name_to_id, id_to_name, id_to_tree, tree_to_id, name_to_tree, tree_to_name = {},{},{},{},{},{}
    all_tree_numbers = list()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(mesh_roots): # exclude animal diseases C22?

                    tree_number = tree_number.text
                    all_tree_numbers.append(tree_number)
                    
                    # Name to Tree
                    try:
                        name = ele.find('DescriptorName').find('String').text
                        name_to_tree.setdefault(name, set()).add(tree_number)
                        tree_to_name.setdefault(tree_number, set()).add(name)
                    except:
                        pass
        except:
            continue        

    all_tree_numbers = sorted(all_tree_numbers)

    name_to_tree = switch_dictset_to_dictlist(name_to_tree)
    with open('output/name_to_tree.json','w') as fout:
        json.dump(name_to_tree, fout)
        
    tree_to_name = switch_dictset_to_dictlist(tree_to_name)
    with open('output/tree_to_name.json','w') as fout:
        json.dump(tree_to_name, fout)
        
    with open('output/tree_numbers.json','w') as fout:
        json.dump(all_tree_numbers, fout)
        
    print('Alignment complete!')
    return name_to_tree, tree_to_name, all_tree_numbers
     

def map_seed_tree_numbers_to_seed_terms(seed_tree_numbers, tree_to_term):
    '''Map a list of lists of MeSH tree numbers to all the MeSH terms that
       correspond to the tree number and that are lower level terms under that
       part of the tree'''
    all_seeds_terms = []
    for seeds in seed_tree_numbers:
        seed_terms = []
        for seed in seeds:
            seed_terms += [term for tree, term in tree_to_term.items() if tree.startswith(seed)]
        seed_terms = list(chain(*seed_terms))
        all_seeds_terms.append(seed_terms)
    return all_seeds_terms
    

def retrieve_ontopic_pmids(categories_of_terms, 
                                  cat_of_pmid_path,
                                  pmid_to_cat_path,
                                  max_num_pmids=99999999999):
    ### Category to PMIDs
    pa = PubMedAPI()
    categories_of_pmids = []
    # For each category...
    for cat_num, category_of_terms in enumerate(categories_of_terms):
        print(cat_num+1, '/', len(categories_of_terms), end='\r')
        category_pmids = []
        category_of_terms = set(category_of_terms)
        # ...identify the MeSH terms defining that category...
        for term in category_of_terms:
            # ...and find the PubMed IDs of documents studying those MeSH topics...
            try:
                response = pa.extract(f'{term}[MeSH Terms]')
                time.sleep(1)
            except:
                print('Error with response for', term)#, 'Trying normal API')
                continue            
            # ...save the PMIDs...
            pmids = [str(pmid) for pmid in response.pmids]
            category_pmids += pmids
            
            # ...but if there is a limit, only save some PMIDs
            if len(category_pmids) > max_num_pmids:
                category_pmids = category_pmids[:max_num_pmids]
                print('Cutting off PMIDs at',max_num_pmids,'\n')
                break
        
        print(f'Category #{cat_num}: {len(set(category_pmids))} PMIDs')
        categories_of_pmids.append(list(set(category_pmids)))
            
    with open(cat_of_pmid_path, 'w') as fout:
        json.dump(categories_of_pmids, fout)

        
    ### PMIDs to Categories
    pmid_to_categories = {}
    for category_num, pmids in enumerate(categories_of_pmids):
        for pmid in pmids:
            pmid_to_categories.setdefault(pmid, []).append(category_num)
        
    with open(pmid_to_cat_path,'w') as fout:
        json.dump(pmid_to_categories, fout)

        
def get_all_pmids(pmid_to_categories_path):
    with open(pmid_to_categories_path,'r') as fin:
        pmid_to_categories = json.load(fin)
    all_pmids = list(set(pmid_to_categories.keys()))
    return all_pmids
    

def submit_pmids_get_article_xml(pmids, topic, retmax=10000):
    '''Submits a batch of PMIDs to PubMed, retrieves the document text,
       outputs it to a file, and repeats the process for a new batch to be
       saved to a new file. These are then iterated through in the next 
       function to extract the PMID, title, and abstract'''
    print('Total', len(pmids), 'submit pmids')
    for batch_num, batch_idx in enumerate(range(0, len(pmids), retmax)):
        pmid_batch = pmids[batch_idx:batch_idx+retmax]
        print(len(pmid_batch), 'in batch ', batch_num)
        pmid_batch = ','.join(pmid_batch)
        
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        params = {
            'db': 'pubmed',
            'rettype': 'abstract',
            'retmode': 'xml',
            'retmax': retmax,
            #'restart': (batch_idx+1)*retmax
        }
        data = {
            'id': pmid_batch,
        }
        response = requests.post(url=base_url, params=params, data=data)
        with open(f'output/{topic}/response_text_{batch_num}.json','w') as fout:
            json.dump(response.text, fout)
        time.sleep(1)
        
            

        
# ALT idea: Iterate through queried PMIDs and extract that part from the XML    
def extract_pmid_title_abstract(topic,
                                queried_pmids,
                                max_num_docs,
                                pmid_to_categories_path, 
                                feature_matrix_outpath,
                                retmax,
                                get_docs_on_pubmed,
                                get_pmids_via_mesh,
                                get_unlabeled_docs,
                                get_labeled_docs=True,
                                verbose=True):
    num_pmids, num_titles, num_abstracts = 0, 0, 0
    
    # Load documents' topic labels
    with open(pmid_to_categories_path,'r') as fin:
        pmid_to_categories = json.load(fin)
    
    # Open final text-to-label output file 
    with open(feature_matrix_outpath,'w') as fout:
        writer = csv.writer(fout)
        writer.writerow(['pmid','title','abstract','labels'])

        # Load batch of documents' xml text
        for batch_num, batch_idx in enumerate(range(0, len(queried_pmids), retmax)):
            pmid_batch = queried_pmids[batch_idx:batch_idx+retmax]
            print('Loading', len(pmid_batch), 'in batch ', batch_num)
            with open(f'output/{topic}/response_text_{batch_num}.json','r') as fin:
                articles_xml_text = json.load(fin)

            # Preprocess the document xml text
            split_on = r'<PMID Version="\d+">'
            pmid_split_xmls = re.split(split_on, articles_xml_text)
            
            # Iterate through each document
            for pmid_split_xml in pmid_split_xmls:

                ### Check whether the document has been MeSH-labeled
                if get_docs_on_pubmed and not get_pmids_via_mesh:
                    if 'MeshHeadingList' in pmid_split_xml:
                        # it's a labeled document
                        if get_labeled_docs:
                            pass
                        elif get_unlabeled_docs:
                            continue                            
                    else:
                        # it's an unlabeled document
                        if get_labeled_docs:
                            continue
                        elif get_unlabeled_docs:
                            pass
                
                ### Extract PMID
                pmid = pmid_split_xml.split('</PMID>')[0] 

                ### Extract title
                try:
                    title_pattern = r'<ArticleTitle>(.*?)<\/ArticleTitle>'
                    title = re.findall(title_pattern, pmid_split_xml)[0].strip('[]')
                    if pmid not in pmid_batch:
                        continue
                    num_titles += 1
                    num_pmids += 1
                    no_title = False
                except:
                    no_title = True
                    pass

                ### Extract abstract
                abstract_pattern = r'<AbstractText.*?>'
                abstract = ''.join(pmid_split_xml.split('</AbstractText>')[:-1])
                # Don't include textless documents
                if len(abstract) == 0:
                    if no_title:
                        continue
                    else:
                        abstract = ''
                else:
                    try:
                        # Extract the abstract from the xml
                        abstract = re.sub(abstract_pattern, '', abstract)
                        abstract_delim = '<Abstract>'
                        other_delim = '<OtherAbstract'
                        if abstract_delim in abstract:
                            abstract = abstract.split(abstract_delim)[1]  
                        elif other_delim in abstract:
                            abstract = abstract.split(other_delim)[1]
                   
                        # Decode abstract
                        html_decoded = html.unescape(abstract)
                        codesc_decoded = codecs.decode(html_decoded.encode('latin1', 'ignore'), 'unicode_escape')
                        abstract = codesc_decoded.replace('\xa0',' ')
                        num_abstracts += 1
                    except:
                        print('ERROR with ', abstract)
                        if no_title:
                            continue
                   
                ### MeSH Topic Labels
                if get_unlabeled_docs: # when getting documents of unknown topics
                    topic_labels = ['nan']
                else:  # when getting labeled on-topic or labeled off-topic
                    try:
                        categories = list(set(pmid_to_categories[pmid]))
                        topic_labels = [str(cat_num) for cat_num in categories]
                    except:
                        if pmid in pmid_batch:
                            print(pmid in pmid_batch,' Alarm: Problem with pmid_to_categories. PMID in queried PMIDs:')
                
                    
                ### Final (Features || Labels)
                writer.writerow([pmid, title, abstract, ','.join(topic_labels)])
            os.remove(f'output/{topic}/response_text_{batch_num}.json')

        if get_unlabeled_docs:
            with open(f'output/{topic}/{topic}_unlabeled_{max_num_docs}_docs_feature_matrix_path.txt', 'w') as fout:
                fout.write(feature_matrix_outpath)
        if verbose:
            print('PMIDs:', num_pmids)
            print('Titles:', num_titles)
            print('Abstracts:', num_abstracts)
    
    

if __name__ == '__main__':
    batch_size = 10000
    
    parser = argparse.ArgumentParser(description='PubMed document API')
    parser.add_argument('--topic', '-t',
                        type=str, default='hf', 
                        help='the name of the topic you are studying (you should have run the first pubmed API with this topic name, and the file should be available)')
    parser.add_argument('--get_pmids_via_mesh', '-gp',
                       action='store_true', default=False)
    parser.add_argument('--get_docs_on_pubmed', '-gd',
                       action='store_true', default=False)
    parser.add_argument('--max_num_docs', '-max',
                        type=int, default=99999999999999)
    parser.add_argument('--get_offtopic_docs', '-off',
                   action='store_true', default=False)    
    parser.add_argument('--get_unlabeled_docs', '-unlab',
                   action='store_true', default=False)
    args = parser.parse_args()
    topic = args.topic
    get_docs_on_pubmed = args.get_docs_on_pubmed
    get_pmids_via_mesh = args.get_pmids_via_mesh 
    ft_mtrx_pth = f'output/{topic}/feature_matrix_{topic}.csv'  # specify if you are retrieving abstracts and titles
    categories_of_pmids_path = f'output/{topic}/category_of_pmids_{topic}.csv'# specify if retrieving PMIDs via MeSH
    pmid_to_categories_path = f'output/{topic}/pmid_to_category_{topic}.json'# specify if not retreiving PMIDs via MeSH or to save a new one
    max_num_docs = args.max_num_docs              # maximum number of documents
    get_labeled_docs = args.get_offtopic_docs     # specify if you want known offtopic docs
    get_unlabeled_docs = args.get_unlabeled_docs  # specify if you want unknown topic docs     
    categories_path = f'input/{topic}_tree_numbers.json'# specify if retreiving PMIDs via MeSH, (lists of lists of tree numbers, descendants are auto included)
    if not os.path.exists(f'output/{topic}'):
        os.makedirs(f'output/{topic}')

    if get_pmids_via_mesh:        
        ### Initial download of all of MeSH (Run once)
        try:
            download_mesh_xml()
            root = parse_mesh_xml()
            print('Done trying to download MeSH tree')
        except: 
            next_year = datetime.now().year+1
            print(f"This year's MeSH file didn't work. Trying next year {next_year}")
            download_mesh_xml(year=next_year)
            root = parse_mesh_xml(year=next_year)
            print('Done trying to download MeSH tree')

        _, _, _ = align_mesh_trees_with_terms(root)

        # Identify categories of study
        tree_to_term = json.load(open('output/tree_to_name.json'))
        seed_tree_numbers = json.load(open(categories_path))  # Categories to study 
        seed_terms = map_seed_tree_numbers_to_seed_terms(seed_tree_numbers, tree_to_term)  # Collects descendant terms

        # Download the PubMed IDs of the categories via an API
        retrieve_ontopic_pmids(seed_terms, 
                               categories_of_pmids_path,
                               pmid_to_categories_path,
                               max_num_docs) 
    
    if get_docs_on_pubmed:
        print('pmid_to_categories_path', pmid_to_categories_path)
        # Extract sections (title, abstract) from downloaded PubMed documents
        all_pmids = get_all_pmids(pmid_to_categories_path)
        print(len(all_pmids), 'all pmids')
        submit_pmids_get_article_xml(all_pmids, 
                                     topic,
                                     batch_size,)
        extract_pmid_title_abstract(topic,
                                    all_pmids,
                                    max_num_docs,
                                    pmid_to_categories_path,
                                    feature_matrix_outpath,
                                    batch_size,
                                    get_docs_on_pubmed,
                                    get_pmids_via_mesh,
                                    get_unlabeled_docs,
                                    get_labeled_docs,
                                    verbose=True) 
        
        