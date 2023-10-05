#!/usr/bin/env python

STR_sample_id = "${sample}"
PATH_contigs="${contigs}"
PATH_proteins="${faa}"
PATH_gff="${gff}"
STR_classif="${classif}"
STR_taxid="${taxid}"
TAB = "\t"

import pandas as pd
import sys
DF_classif = pd.read_csv(PATH_classif)

### I. Read in GTDB-Tk classifications ######  -> DF_classif ________________________________________________________

DF_classif = pd.read_csv(PATH_classif)
DF_out = DF_classif

# Keep only classified SAGS
DF_classif = DF_classif.loc[DF_classif['classification']!='0']
DF_classif[['domain','phylum','class','order','family','genus','species']] = DF_classif['classification'].str.split(';', expand=True) # Break ranks into columns
DF_classif = DF_classif[['sag', 'species',  'genus', 'family', 'order', 'class', 'phylum', 'domain']] # reverse order


### II. Read in GTDB lineages, linked to taxIDS ##### --> DF_names _______________________________________________________
DF_names = pd.read_csv(PATH_lineage_names, sep=tab)
DF_names[['domain','phylum','class','order','family','genus','species']] = DF_names['lineage'].str.split(';', expand=True) # Break ranks into columns

### III. cross-reference -> DICT_sag2taxid ______________________________________________________
#(by matching first at the species, then genus if necessary, etc.)

verbose = False

### Start writing output file
#DF_classif = DF_classif.head()

DICT_sag2taxid = {}

species_found = 0; genera_found = 0; families_found = 0; orders_found = 0; classes_found = 0; phyla_found = 0; domains_found = 0; not_found = 0

LIST_gtdb_names_not_found_in_namesdmp = []

for index, row in DF_classif.iterrows():
    DF_this_row = row
    SAG = DF_this_row['sag']
    taxid = ''; mapped_to_rank = ''; TAXID_species = ''; TAXID_genus = ''; TAXID_family = ''

    ### I ### Can we find this SPECIES?
    rank = "species"
    STR_species = DF_this_row[rank]
    if STR_species in DF_names[rank].to_list():
        DF_name_hits = DF_names.loc[DF_names[rank] == STR_species]

        # print(f'{len(DF_name_hits)} match found')                       # How many species in Names.dmp matched this?
        TAXID_species = DF_name_hits.iloc[0]['taxid']                   # Take the taxid of first match (there are rare cases [e.g. 3 in a million] of two different taxids having the same 'name_txt')
        if verbose == True: print("species "+STR_species+" = "+TAXID_species)
        species_found += 1; taxid= TAXID_species; name=STR_species; mapped_to_rank = rank
    else:
        if STR_species != 's__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_species)
            
        ### II ### If not, can we find this GENUS?
        rank = "genus"
        STR_genus = DF_this_row[rank]
        if STR_genus in DF_names[rank].to_list():
            DF_name_hits = DF_names.loc[DF_names[rank] == STR_genus]
            TAXID_genus = DF_name_hits.iloc[0]['taxid']
            if verbose == True: print("genus "+STR_genus+" = "+TAXID_genus)
            genera_found += 1; taxid= TAXID_genus; name=STR_genus; mapped_to_rank = rank
        else:
            if STR_genus != 'g__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_genus)
                
            ### III ### If not, can we find the FAMILY?
            rank = "family"
            STR_family = DF_this_row[rank]
            if STR_family in DF_names[rank].to_list():
                DF_name_hits = DF_names.loc[DF_names[rank] == STR_family]
                TAXID_family = DF_name_hits.iloc[0]['taxid']
                if verbose: print("family "+STR_family+" = "+TAXID_family)
                families_found += 1; taxid = TAXID_family; name=STR_family; mapped_to_rank = rank
            else:
                if STR_family != 'f__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_family)
                    
                ### IV ### If not, can we find the ORDER?
                rank = "order"
                STR_order = DF_this_row[rank]
                if STR_order in DF_names[rank].to_list():
                    DF_name_hits = DF_names.loc[DF_names[rank] == STR_order]
                    TAXID_order = DF_name_hits.iloc[0]['taxid']
                    if verbose: print("order "+STR_order+" = "+TAXID_order)
                    orders_found +=1; taxid = TAXID_order; name=STR_order; mapped_to_rank = rank     
                else:
                    if STR_order != 'o__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_order)
                        
                    ### V ### If not, can we find the CLASS?
                    rank = "class"
                    STR_class = DF_this_row[rank]
                    if STR_class in DF_names[rank].to_list():
                        DF_name_hits = DF_names.loc[DF_names[rank] == STR_class]
                        TAXID_class = DF_name_hits.iloc[0]['taxid']
                        if verbose: print("class " + STR_class + " " + TAXID_class)
                        classes_found += 1; taxid = TAXID_class; name=STR_class; mapped_to_rank = rank    
                    else:
                        if STR_class != 'c__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_class)
                            
                        ### VI ### If not, can we find the PHYLUM?
                        rank = "phylum"
                        STR_phylum = DF_this_row[rank]
                        if STR_phylum in DF_names[rank].to_list():
                            DF_name_hits = DF_names.loc[DF_names[rank] == STR_phylum]
                            TAXID_phylum = DF_name_hits.iloc[0]['taxid']
                            if verbose: print("order "+STR_phylum+" = "+TAXID_phylum)
                            phyla_found += 1; taxid = TAXID_phylum; name=STR_phylum; mapped_to_rank = rank         
                        else:
                            if STR_phylum != 'p__': LIST_gtdb_names_not_found_in_namesdmp.append(STR_phylum)
                            print("SAG "+SAG+" not placed to phylum")
                            not_found += 1
                            
    DICT_sag2taxid[SAG] = str(taxid)
                            
##### Summarize
print('species found: '+str(species_found))
print('genera found: '+str(genera_found))
print('families found: '+str(families_found))
print('orders found: '+str(orders_found))
print('classes found: '+str(classes_found))
print('phyla found: '+str(phyla_found))

if len(LIST_gtdb_names_not_found_in_namesdmp) > 0:
    sys.exit("gtdb taxon names that were not found in "+PATH_lineage_names+": "+",".join(LIST_gtdb_names_not_found_in_names))
    
#### IV. ####### Save output ###### ______________________________________________________________________________________\

SER_taxids = DF_out['sag'].map(DICT_sag2taxid)
DF_out['taxid'] = SER_taxids ### Adds in a new 'taxid' column with the findings
DF_out.to_csv(PATH_out)
