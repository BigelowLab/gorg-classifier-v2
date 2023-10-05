#!/usr/bin/env python

PATH_out = "${sample}_contigs_proteins_annot.tsv"

STR_sample_id = "${sample}"
PATH_contigs="${contigs}"
PATH_proteins="${faa}"
PATH_gff="${gff}"
STR_classif="${classif}"
STR_taxid="${taxid}"
TAB = "\t"

import pandas as pd
from Bio import SeqIO
import sys

#### I. Convert prokka gff -> tsv ###### _______________________________________________________________________

with open(PATH_gff) as gff_fh, open(STR_sample_id + "_prokka_gff_tabulated.tsv","w") as output:
	## print header to output tsv
	print("contig_id","locus_tag","ftype","start","end","strand","prokka_gene","prokka_EC_number","prokka_product","swissprot_gene","swissprot_EC_number","swissprot_product","swissprot_eggNOG","swissprot_KO","swissprot_Pfam","swissprot_CAZy","swissprot_TIGRFAMs",sep=TAB,file=output)

	for line in gff_fh:

		### initialize column values
		locus_tag=''; Name=''; prokka_gene=''; prokka_EC_number=''; prokka_product=''; swissprot_gene=''; swissprot_EC_number=''; swissprot_product=''; swissprot_eggNOG=''; swissprot_KO=''; swissprot_Pfam=''; swissprot_CAZy='';  swissprot_TIGRFAMs=''

		if line.startswith("##FASTA"):
			break
		if line.startswith("#"):
			continue
		LIST_of_fields = line.strip().split(TAB)

		### Fill required columns
		#contig_id = LIST_of_fields[0]
		# change contig_id from "gnl|Prokka|AM-276-A01_1" -> "SCGC_AM-276-A01_contig1"
		contig_id = LIST_of_fields[0].replace('gnl|Prokka|','')
		# contig_id = LIST_of_fields[0].replace('gnl|Prokka|','SCGC!').replace('_','_contig').replace('!','_')
		ftype = LIST_of_fields[2]
		start = LIST_of_fields[3]
		end = LIST_of_fields[4]
		strand = LIST_of_fields[6]

		## Skip features that are just genes
		if ftype == "gene": continue

		## Read annots into a dictionary
		DICT_of_annots = {}
		LIST_of_annots = LIST_of_fields[-1].split(';')   #annots are always in last field output by prokka
		for i in LIST_of_annots:
			split_keyvalue = i.split('=') # key is on right of '=' value is onleft
			if len(split_keyvalue) != 2 and MODE_verbose == True:
				print(split_keyvalue)
			if len(split_keyvalue) == 2:
				key = split_keyvalue[0]
				value = split_keyvalue[1]
				DICT_of_annots[key] = value # load into dictionary

		## The 'product' annotation may contain additional annots (eggNOG, KO, Pfam, CAZy, & TIGRFAMs) all sep by '^^'
		if 'product' in DICT_of_annots.keys():
			STR_product = DICT_of_annots.get('product')
			if "^^" not in STR_product:
				prokka_product = STR_product

			else:
				## Load in those additional annotations
				STR_product = STR_product.replace("%2",",")
				for LIST_additional_annots in STR_product.split("^^"):
					items = LIST_additional_annots.split("::")
					DICT_of_annots[items[0]] = items[1]

				prokka_product = DICT_of_annots.get('product')

		## Fill optional columns
		if 'locus_tag' in DICT_of_annots.keys(): locus_tag = DICT_of_annots.get('locus_tag')
		if 'gene' in DICT_of_annots.keys(): prokka_gene = DICT_of_annots.get('gene')
		if 'eC_number' in DICT_of_annots.keys(): prokka_EC_number = DICT_of_annots.get('eC_number')
		if 'eggNOG' in DICT_of_annots.keys(): swissprot_eggNOG = DICT_of_annots.get('eggNOG')
		if 'KO' in DICT_of_annots.keys(): swissprot_KO = DICT_of_annots.get('KO')
		if 'Pfam' in DICT_of_annots.keys(): swissprot_Pfam = DICT_of_annots.get('Pfam')
		if 'CAZy' in DICT_of_annots.keys(): swissprot_CAZy = DICT_of_annots.get('CAZy')
		if 'TIGRFAMs' in DICT_of_annots.keys(): swissprot_TIGRFAMs = DICT_of_annots.get('TIGRFAMs')

		print(contig_id,locus_tag,ftype,start,end,strand,prokka_gene,prokka_EC_number,prokka_product,swissprot_gene,swissprot_EC_number,swissprot_product,swissprot_eggNOG,swissprot_KO,swissprot_Pfam,swissprot_CAZy,swissprot_TIGRFAMs,sep=TAB,file=output)

        
##### II. Prokka table -> DF ######## ___________________________________________________________________

DF_prokka=pd.read_table(STR_sample_id + "_prokka_gff_tabulated.tsv", sep=TAB)
if 'CDS' not in DF_prokka['ftype'].unique():
    sys.error("No proteins were predicted from this SAG, so it will not be used in the classifier")
    
DF_prokka=DF_prokka.rename(columns={'locus_tag':'protein_id'})
DF_prokka['contig_id'] = DF_prokka['contig_id'].str.replace('-','').str.replace('_','_C')

##### III.  Add DNA sequences to DF ###### ___________________________________________________
LIST_contigs = []
for record in SeqIO.parse(PATH_contigs, "fasta"): LIST_contigs.append([record.id, str(record.seq)])
DF_contigs = pd.DataFrame(LIST_contigs, columns=["contig_id", "contig_sequence"])

DF_prokka_contigs = DF_prokka.merge(DF_contigs, on='contig_id', how='left')

##### IV. Add protein sequences to DF ###### _________________________________________________
LIST_proteins = []
for record in SeqIO.parse(PATH_proteins, "fasta"): LIST_proteins.append([record.id, str(record.seq)])
DF_proteins = pd.DataFrame(LIST_proteins, columns=["protein_id", "protein_sequence"])
DF_all = DF_prokka_contigs.merge(DF_proteins, on='protein_id', how='left')

DF_all.to_csv(PATH_out, index=False)
