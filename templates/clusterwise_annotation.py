#!/usr/bin/env python


### Inputs _____________________________________________________


STR_dataset = "${params.db}"
STR_basename = "${sample}"
STR_gtdb_release = "${VAL_gtdb_release}"

# STR_dataset = "gorg_dark"
# STR_basename = "100ksub"
# STR_gtdb_release = "R207"
# BASEDIR = "./input_" + STR_dataset + "/"


PATH_out = STR_basename + "_clusters_v_" + STR_dataset + "_GTDB_" + STR_gtdb_release + ".tsv.gz" 

### Auto-inputs
PATH_in = "${classified_reads}"
PATH_annots = "${annotations}"
PATH_lineages = "${lineages}"
PATH_clusters = "${clusters}"

##PATH_in=STR_basename+ "_reads_classified_v_"+ STR_dataset + "_v_"+ STR_gtdb_release + ".tsv.gz" 
# PATH_clusters = BASEDIR + STR_dataset + "_60minid_m80.tsv"
# PATH_annots = BASEDIR + "annotations.tsv"
# PATH_lineages = BASEDIR + "gtdb_R207_lineages.tsv"

NAME_lineage_column = "gtdb_" + STR_gtdb_release + "_lineage"

from collections import Counter
import pandas as pd
pd.set_option("display.max_rows", 6)
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
import numpy as np
import time

TAB = "\\t" ### change to "\\t" if running this code in nextflow

# Column names to use for the output file (based on the input file)
names=pd.read_csv(PATH_in,compression='gzip', nrows=0, sep=TAB)
header=names.columns.to_list()

def which_rank_encompasses_cluster(STR_num_taxa_per_rank):
    '''
    Given a colon-separated string like "1:1:1:1:1:2:5" (from domain on left to species on right)
    This function determines the lowest rank encompassing the diversity.
    
    Examples:
    'family' for "1:1:1:1:1:2:5" because a single family contains 2 genera and 5 species
    Likewise:
        "species for "1:1:1:1:1:1:1"
        "genus" for "1:1:1:1:1:1:4"
        "phylum" for "1:1:2:5:8:9"
    '''
    DICT_step_to_rank = {1:'prokaryotes',2:'domain',3:'phylum',4:'class',5:'order',6:'family',7:'genus',8:'species'}
    rank = 0
    LIST_num_taxa_per_rank = STR_num_taxa_per_rank.split(':') # Turn the string into a list
    LIST_num_taxa_per_rank.append("1") 
    for num_taxa in LIST_num_taxa_per_rank:
        rank +=1
        if num_taxa != '1': break  
    STR_encompassing_rank= DICT_step_to_rank[rank] # Determines the encompassing rank based on how many steps it took until num_taxa became more than 1
    return STR_encompassing_rank



#### I. #### Set up: (DF_annots ) _______________________________________________________________________________

print("Reading in lineage and annotation info"); print()

_DF = pd.read_csv(PATH_lineages, sep=TAB)
DICT_taxid_to_lineage = dict(zip(_DF['taxid'],_DF['lineage']))

DF_annots = pd.read_csv(PATH_annots, sep=TAB).drop(columns=["Unnamed: 0"])

DF_clust = pd.read_csv(PATH_clusters, sep=TAB, names=['gene_cluster_rep', 'gene'])


#### II. ##### Prepare protein cluster info --> DF2_clust ________________________________________________________

print("Contextualizing protein clusters with taxonomic info"); print()
DF_clust['TEMP_taxid']=[i.split("_")[-1] for i in DF_clust['gene']]                         # Extract taxid
DF_clust['TEMP_lineage'] = [DICT_taxid_to_lineage[int(i)] for i in DF_clust['TEMP_taxid']]  # Retrieve lineage

LIST_ranks = ['domain','phylum','class','order','family','genus','species']

### Split TEMP_lineage into cols "TEMP_domain" "TEMP_phylum" etc.
DF_clust[['TEMP_'+i for i in LIST_ranks]] = DF_clust["TEMP_lineage"].str.split(";",expand=True)

### Make columns called "TEMP_domain_count" "TEMP_phylum_count" etc.
for STR_rank in LIST_ranks:
    _DF_this_rank_per_cluster = DF_clust[DF_clust['TEMP_'+STR_rank] != 'Unclassified'].drop_duplicates(subset=['gene_cluster_rep','TEMP_'+STR_rank]).groupby('gene_cluster_rep', as_index=False)['gene'].count().rename(columns={'gene':'TEMP_'+STR_rank+'_count'})
    DF_clust = DF_clust.merge(_DF_this_rank_per_cluster[['gene_cluster_rep','TEMP_'+STR_rank+'_count']], how='left')  
    
### Concatenate those count values into a string like "1:1:1:1:1:1:1" (if it has one taxon each for domain, phylum, etc.)
DF_clust['TEMP_count_per_rank'] = DF_clust['TEMP_domain_count'].astype(str) +  ':' + DF_clust['TEMP_phylum_count'].astype(str) +  ':' + DF_clust['TEMP_class_count'].astype(str) +  ':' + DF_clust['TEMP_order_count'].astype(str) +  ':' + DF_clust['TEMP_family_count'].astype(str) +  ':' + DF_clust['TEMP_genus_count'].astype(str) +  ':' + DF_clust['TEMP_species_count'].astype(str)
    
### Based on that string, find the rank that encompasses the cluster
DF_clust['rank_encompassing_cluster'] = DF_clust['TEMP_count_per_rank'].apply(which_rank_encompasses_cluster)

### Report the taxonomic breadths of the protein clusters
_SER = DF_clust['rank_encompassing_cluster'].value_counts()
for i in list(_SER.keys()):
    PERC = str(round(_SER[i] / len(DF_clust) * 100, 2))
    print(PERC + "% of proteins cluster within a single "+i)
    
DF2_clust = DF_clust.loc[:, ~DF_clust.columns.str.startswith('TEMP')]  # Drop TEMP cols



#### III. ##### Cluster the Kaiju_hits+gtdb+prokka and reannotate cluster-wise --> DF_summary _________________________

start_time = time.time()

# specify columns that need to have whitespaces stripped
cols = [NAME_lineage_column, 'gene_cluster_rep']

chunksize = 100000
counter=1   ## counter to track progress
check=1

### Read the file in chunks and group the reads by each gene/tax lineage for each chunk
for chunk in pd.read_csv(PATH_in, compression='gzip', sep=TAB, chunksize=chunksize, skiprows=0, names=header):
    if counter % 10 == 0:
        print('The number of processed reads is',counter*chunksize, '. The total number of classified reads is', DF_summary.number_of_reads.sum(), DF_summary.number_of_reads.sum()/(counter*chunksize)*100, '% of reads recruiting.')
        
    ### Process first chunk -> DF_summary
    if check==1:
        # Extract classified reads only
        DF_summary = chunk[(chunk['status'] == 'C')]
        if len(DF_summary)==0:
            # For some of the samples there are no classified reads. This check makes sure that this first part is repeated
            # until there is at least one classified read and then it moves on to the second part that stacks the chunks.
            check=1
            continue
            
        else:
            check=2
            
            # Extract the first gene (contig;start;stop) that the read is assigned to (discarding other potential best hits)
            DF_summary['accessions_best_matches']=[i.split(",")[0] if "," in i else i for i in DF_summary['accessions_best_matches']]

            # Append the taxid to the gene. (contig;start;stop) --> (contig;start;stop_taxid)
            DF_summary['accessions_best_matches'] = DF_summary['accessions_best_matches'].str.split('_').str[0] + '_' + DF_summary['accessions_best_matches'].str.split('_').str[1] + '_' + DF_summary['taxid'].astype(str)
            
            # OLD: for these samples only I have to modify the gene name to add the '-' back in
            # OLD: DF_summary['accessions_best_matches']=DF_summary['accessions_best_matches'].str[:2] + '-' + DF_summary['accessions_best_matches'].str[2:5] + '-'+DF_summary['accessions_best_matches'].str[5:]
        
            # Link each gene with its cluster by merging DF_summary + DF2_clust (cluster is identified by its 'gene_cluster_rep')
            DF_summary=DF_summary.merge(DF2_clust, left_on='accessions_best_matches', right_on='gene', how='left')

            # Extract the gene annotation data from the summary file and drop duplicates
            genes_summary=DF_summary[['gene_cluster_rep', 'rank_encompassing_cluster','prokka_EC_number', 'prokka_product', 'swissprot_gene',
                                      'swissprot_EC_number', 'swissprot_eggNOG', 'swissprot_KO', 'swissprot_Pfam',
                                      'swissprot_CAZy', 'swissprot_TIGRFAMs']].copy()
            genes_summary.drop_duplicates(inplace=True)
            
            
            # Group the DF by taxonomic lineage, and gene cluster
            DF_summary=DF_summary.groupby([NAME_lineage_column, 'gene_cluster_rep'], as_index=False)['status'].count() #NAME_lineage_column
            DF_summary.rename(columns={'status':'number_of_reads'}, inplace=True)
                            
            
    # Read all other chunks into a temporary df and process them
    else:
        # Same as above but store each successive chunk as a tempory_df
        tdf_summary = chunk[(chunk['status'] == 'C')]
        if len(tdf_summary)==0:
            # Proceed to the next chunk if there are no classified reads
            continue
        
        else:               
            
            # Extract the first gene (contig;start;stop) that the read is assigned to (discarding other potential best hits)
            tdf_summary['accessions_best_matches']=[i.split(",")[0] if "," in i else i for i in tdf_summary['accessions_best_matches']]

            # Append the taxid to the gene. (contig;start;stop) --> (contig;start;stop_taxid)
            tdf_summary['accessions_best_matches'] = tdf_summary['accessions_best_matches'].str.split('_').str[0] + '_' + tdf_summary['accessions_best_matches'].str.split('_').str[1] + '_' + tdf_summary['taxid'].astype(str)
        
            tdf_summary=tdf_summary.merge(DF2_clust, left_on='accessions_best_matches', right_on='gene', how='left')
    
            # extract all the gene information from this chunk and drop duplicates
            tgenes_summary=tdf_summary[['gene_cluster_rep', 'rank_encompassing_cluster', 'prokka_EC_number', 'prokka_product', 'swissprot_gene',
                                      'swissprot_EC_number', 'swissprot_eggNOG', 'swissprot_KO', 'swissprot_Pfam', 
                                      'swissprot_CAZy', 'swissprot_TIGRFAMs']].copy()
            tgenes_summary.drop_duplicates(inplace=True)
            
            tdf_summary=tdf_summary.groupby([NAME_lineage_column, 'gene_cluster_rep'], as_index=False, dropna=False)['status'].count()
            tdf_summary.rename(columns={'status':'number_of_reads'}, inplace=True)
            
            # Merge the temporary df with the summary df and convert nan to 0 (this allows for adding together the columns)
            DF_summary=DF_summary.merge(tdf_summary, on=[NAME_lineage_column, 'gene_cluster_rep'], how='outer')
            DF_summary["number_of_reads_x"].fillna(0, inplace=True)
            DF_summary["number_of_reads_y"].fillna(0, inplace=True)
    
            # Add the read numbers, and drop the read number columns that were duplicated upon the merge
            DF_summary['number_of_reads']=DF_summary['number_of_reads_x']+DF_summary['number_of_reads_y']
            DF_summary.drop(columns=['number_of_reads_x', 'number_of_reads_y'], inplace=True)
            
            # Stack gene data with initial gene data and drop duplicates. This will be appended at the end
            genes_summary=pd.concat([genes_summary, tgenes_summary], ignore_index=True, axis=0)
            genes_summary.drop_duplicates(inplace=True)
            

    counter+=1
    
# Now re-annotate each gene cluster based on its first annotated product.
if check==2:
    
    genes=genes_summary['gene_cluster_rep'].unique()
    final_genes_summary=pd.DataFrame()
    for i in genes:
        tdf_genes=genes_summary[genes_summary['gene_cluster_rep']==i]
        ## If there is only one product use that
        if len(tdf_genes['prokka_product'].unique())==1:
            final_genes_summary=pd.concat([final_genes_summary, tdf_genes])
        ## If there are more than one product get rid of the hypothetical
        else:
            tdf_genes=tdf_genes[tdf_genes['prokka_product']!='hypothetical protein']
            if len(tdf_genes['prokka_product'].unique())>1:
                ## if there are still multiple products use only the first one
                tdf_genes=tdf_genes.iloc[0]
                final_genes_summary=pd.concat([final_genes_summary, tdf_genes])
                
            else:
                final_genes_summary=pd.concat([final_genes_summary, tdf_genes])
            
    
    DF_summary=DF_summary.merge(final_genes_summary, on='gene_cluster_rep', how='left')
    DF_summary.to_csv(PATH_out, compression="gzip")


print("It took "+str(round(time.time() - start_time,2))+" seconds to add secondhand annotations to the reads")