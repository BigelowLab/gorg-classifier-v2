#!/usr/bin/env python

import csv
import gzip
import pandas as pd
import numpy as np
import sys
import warnings
# Silence pandas warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
pd.options.mode.chained_assignment = None  # default='warn'

STR_dataset = "${params.db}"
STR_basename = "${sample}"
STR_gtdb_release = "${VAL_gtdb_release}"

PATH_output = STR_basename + "_reads_classified_v_" + STR_dataset + "_GTDB_" + STR_gtdb_release + ".tsv.gz" 

### Auto-inputs
PATH_kaiju = "${kaiju_hits}"
PATH_annots = "${annotations}"
PATH_lineages = "${lineages}"
PATH_clusters = "${clusters}"

# PATH_kaiju = "BASEDIR + '100ksub_hits.txt'"
# PATH_annots = BASEDIR + "annotations.tsv"
# PATH_lineages = BASEDIR + "gtdb_R207_lineages.tsv"
# PATH_clusters = BASEDIR + STR_dataset + "_60minid_m80.tsv"

PATH_temp = "TEMP_lineages_added.tsv.gz"

TAB = "\\t"  ### change to "\\t" if running this code in nextflow


### Read in lineage table -> DICT_taxid_to_lineage
DF = pd.read_csv(PATH_lineages, sep=TAB)
print("Reading " + str(len(DF)) + " GTDB lineages into dictionary")
DICT_taxid_to_lineage = dict(zip(DF['taxid'], DF['lineage']))

### Read in annotation table -> DICT_genepos_to_annots 
DICT_genepos_to_annots = {} 
'''
a dictionary with a format like:
        
        KEY ("genepos")      VALUE (list of 10 "annots")
 _________________
 Contig;Start;Stop      :     [prokka_gene,
                                 prokka_EC_number,
                                 prokka_product,
                                 swissprot_gene,
                                 swissprot_EC_number,
                                 swissprot_eggNOG,
                                 swissprot_KO,
                                 swissprot_Pfam,
                                 swissprot_CAZy,
                                 swissprot_TIGRFAMs']
                                 
    EXAMPLE:
 
 {'AH-770-A01_N1;4070;4606':  ['1.4.9.2',
                              'Aralkylamine dehydrogenase light chain',
                              'DHML_METFK',
                              '1.4.9.1',
                              'ENOG4108XDC,CENOG4111IC4',
                              'K15228',
                              'PF02975',
                                '',
                              'TIGR02659']}
'''
keep = []
with open(PATH_annots, 'r') as fh:
    header = fh.readline().strip().split(TAB)  
    fh.seek(0)
    keep = header[7:]
    reader = csv.DictReader(fh, delimiter=TAB) 
    for row in reader:
        DICT_genepos_to_annots[f"{row['contig_id']};{row['start']};{row['stop']}"] = [row[i] for i in keep]
print("Reading " + str(len(DICT_genepos_to_annots)) + " reference annotations into dictionary")


### I. Kaiju hits -- > DF_kaiju ______________________________________________________________________
LIST_kaiju_colnames = ['status','read_name','taxid','score_best_match','taxids_best_matches','accessions_best_matches','sequence_matching_fragment']
DF_kaiju = pd.read_csv(PATH_kaiju, sep=TAB, names=LIST_kaiju_colnames)


### II. Add lineages to kaiju hits -- > DF2_lineages_added ___________________________________________

NAME_lineage_column = "gtdb_" + STR_gtdb_release + "_lineage"
DF_unclassified = DF_kaiju.loc[DF_kaiju['taxid']==0]
DF_unclassified[NAME_lineage_column]=np.nan

# Only add lineages to classified reads
DF_classified = DF_kaiju.loc[DF_kaiju['taxid']!=0]
print("Analyzing "+str(len(DF_classified))+" hits in kaiju for "+str(len(DF_kaiju))+" reads")

print("Matching lineages to taxid of each hit")
DF_recognized_taxid = DF_classified.loc[DF_classified['taxid'].isin(list(DICT_taxid_to_lineage.keys()))]
DF_recognized_taxid[NAME_lineage_column] = [DICT_taxid_to_lineage[i] for i in DF_recognized_taxid['taxid']]

# Check if any of the taxids didn't link to a lineage
DF_unrecognized_taxid = DF_classified.loc[~DF_classified['taxid'].isin(list(DICT_taxid_to_lineage.keys()))]
if len(DF_unrecognized_taxid) > 1:
    print(str(len(DF_unrecognized_taxid)) + ' however, had unrecognized taxids and recieved no lineage.')
    print("These unrecognized taxids and counts were...")
    print(DF_unrecognized_taxid['taxid'].value_counts())
    sys.exit("Some taxids do not match. Check stdout message.")
    
# Put everything back together and save to TEMP intermediate file
DF2_lineages_added = pd.concat([DF_recognized_taxid, DF_unclassified], sort=False)
DF2_lineages_added = DF2_lineages_added.sort_index()
DF2_lineages_added.to_csv(PATH_temp,sep=TAB, header=False, index=False, compression="gzip")


### III. Add annotations --> PATH_output _________________________________________________________________________________

# Read in the TEMP file and add annotation columns via DICT_genepos_to_annots 
print(str("Matching annotations to gene name;start;stop of each hit"))

with gzip.open(PATH_temp, 'rt') as in_fh, gzip.open(PATH_output,"wt") as out_fh:
    output_header = [
    'status',
    'read_name',
    'taxid',
    'score_best_match',
    'taxids_best_matches',
    'accessions_best_matches',
    'sequence_matching_fragment',
    NAME_lineage_column,
    ]
    output_header.extend(keep)
    print(*output_header, sep=TAB, file=out_fh)

    for line in in_fh:
        toks = line.strip("\\r\\n").split(TAB)                ### gremlin
        if toks[0] != "U":
            contig_id = toks[5].partition(",")[0]
            toks.extend(DICT_genepos_to_annots[contig_id])
        print(*toks, sep=TAB, file=out_fh)

        
### IV. Report ___________________________________________________________________________________________________________
DF4 = pd.read_csv(PATH_output, sep=TAB, compression="gzip")
print(str(len(DF4.loc[~DF4['prokka_product'].isna()])) + " reads given annotations")
if 'hypothetical protein' in DF4['prokka_product'].unique(): print(str(DF4['prokka_product'].value_counts()['hypothetical protein']) + " of which are 'hypothetical protein'")
print("Done")
