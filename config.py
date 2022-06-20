"""
This file would hold all of the constant values needed that would be used in the project
"""
import os

working_dir = os.path.dirname(os.path.abspath(__file__))

PSEUDOMONAS = 'Pseudomonas_aeruginosa'
STRAIN_P = 'PAO1'

RAW_DATA_P = os.path.join(working_dir, "PAO1_raw_data.csv")

CDS_FILE_P = os.path.join(working_dir, "Data", "GCF_000006765.1_ASM676v1_cds_from_genomic.fna.gz")

# TF file, containing all the TF relevant for the trail in a table:
# TF | operon | target | regulation type | technology | article pub id | strain
NETWORK_FILE = os.path.join(working_dir, "Data", "NewNetwork_250821.txt")

# The directory of all the data (has all the trails directories)
HIGH_LOST_PATH_P = os.path.join(working_dir, "Data", "{0}".format(PSEUDOMONAS))

EXTRA_FILES = os.path.join(working_dir, "Data", "additional_annotations_converted.txt")

# This file contains all the paths that are relevant for this specific virus
RESULT_PROPER_FILE = os.path.join(working_dir, "Data", "results_proper.txt")

UP_REG_COLOR = 'thistle'

DOWN_REG_COLOR = 'hotpink'

TABLE_S1 = os.path.join(working_dir, "Data", "Table_S1.txt")

# This file contains all the propagation and heterogeneity results
SUPPLEMENTARY_TABLE_1 = os.path.join(working_dir, "{0}_S1.csv".format(STRAIN_P))

# KEGG raw data
KEGG_PATH = os.path.join(working_dir, "Data", "pae.list")
CONVERT_KEGG_PATH = os.path.join(working_dir, "Data", "map_title.tab")

# GO raw data
GO_PATH = os.path.join(working_dir, "Data", "Pseudomonas_aeruginosa_go_annotations.txt")
CONVERT_GO_PATH = os.path.join(working_dir, "Data", "go.obo")

