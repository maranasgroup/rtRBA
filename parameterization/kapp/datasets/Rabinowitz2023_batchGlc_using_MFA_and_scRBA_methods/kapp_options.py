# Model-specific settings for kapp calculations, as well as steps that must be run each time
import pandas as pd
import numpy as np
import requests,json,os,shutil,sys

model_root_path = '../../../../'
path_model = model_root_path+'GAMS/model/'
pycore_path = model_root_path+'pycore/'
for path in [model_root_path, pycore_path]:
    if path not in sys.path:
        sys.path.append(path)
# retrieve name of the directory containing the options file
dir_name = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

from pycore.utils import metabolites_dict_from_reaction_equation_RBA, extract_details_from_rxnid
from pycore.simulate import * 
from pycore.gsm_custom_functions import *
from pycore.cobrapy_custom_extras import *
vmax = 1e3 # max flux in either direction, before applying nscale
nscale=1e5 # will be used as scale when running GAMS
nscale2=1e5 # 2nd scale, if different one needed; otherwise, set to same as nscale

# Load path
path_gen = model_root_path+'build_model/'
path_gams = model_root_path+'GAMS/'

sm_path = path_gen + 'input/iRhtoBatch-Rabinowitz.json'

prot_path = path_gen + 'model/PROTEIN_stoich_curation.tsv'
unknown_prot_path = path_gen + 'input/PROTEIN_stoich_curation_unknown.xlsx'
model_xlsx_path = path_model + 'RBA_stoichiometry.tsv'
ribonuc_path = path_gen + 'input/RIBOSOME_nucleus.xlsx'
ribomito_path = path_gen + 'input/RIBOSOME_mitochondria.xlsx'
metab_rxns_path = path_model + 'RBA_rxns_rxnmetabolicnetwork.txt'
gsm_rxn_ids_path = path_model + 'SM_rxn_ids.txt'
sij_path = path_model + 'RBA_sij.txt'
aa_mapping_path = path_gen + 'input/PROTEIN_amino_acid_map.txt'
aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

flux_path = '../raw_data_files/Rabinowitz-flux.xlsx'
fluxes_to_ignore = ['Ht_c_e','EX_h_e','ALAt_c_m'] # Excluded 'Ht_c_e','EX_h_e' since we're unsure if this species acidifies its environment, or proton export is due to other reasons (e.g., missing ATP maintenance systems, Ht_c_e being used to lump together proton transport from other processes, etc.), and 'ALAt_c_m' since its direction in the MFA flux appears to be wrong according to yeast9, yeast8.3.4, iRhtoC, and scRBA.
# col_LB = 'mfaLB'
# col_UB = 'mfaUB'
col_LB = 'pfba'
col_UB = 'pfba'

# must match the growth rate in other files
mu = 0.38
# mu = 0.3865768 
# protein fraction (disable by uncommenting "ptot = 1" unless composition varies w/ growth rate)
ptot = 0.468
# ptot = 0.4555 

# protein data
df_raw = pd.read_excel('../raw_data_files/Rabinowitz-BatchGlc-abridged.xlsx',
                         sheet_name='for-kapps-v4', usecols=[0,1,2,3])
df_raw.index = df_raw['TRUE best match'].to_list() # name of protein
# If encountering issues with averaging multiple replicates, could try using only one column in case the average results in an infeasible AA composition
cols_data = ['B_frac_final'] # where protein abundance data is stored
uniprot_col = '' # set to '' if no column with accession names provided
proteomics_units = 'g/g protein' # supported options: 'g/g protein' (originally built with this in mind), 'g/gDW'
use_ribo_data = False # False to ignore ribosome subunit abundance data (still used for MPFS constraints)

# if True, proteins with measured translation rates of 0 can be made at a flux equal to prosynSlackUB, thus discouraging their use where possible. If you want such proteins to be made, you thus may have to specify them in the phenotype.txt file as having fixed values for prosynSlackUB. 
# if False, proteins with measured translation rates of 0 can't be made. This may cause issues if no isozymes exist for such proteins but their rxns are essential, and such proteins may still be present but unmeasured due to measurement errors.
allow_trans_when_measurement_is_0 = True 

nonmodeled_proteome_allocation = 0
# nonmodeled_proteome_allocation = 0.503973558
dummy_protein = {'id':'PROSYN-PROTDUMMY','AA abundances':dict()}

path_data = './protein_data.tsv'
if os.path.exists(path_data):
    df_data = read_spreadsheet(path_data)
    df_data.index = df_data['id'].to_list()
    df_data = df_data[df_data['conc (g/gDW)'] > 0]
    # Excluding ribosome protein subunit (conflicting if fit to both enzymatic and ribosomal protein data)
    # make perfect copy of df_data
    df_data_full = df_data.copy()
    if not use_ribo_data:
        df_data = df_data[(df_data.type == 'truedata_enz') | (df_data.type == 'gapfill_subunit')]

biom_id = 'BIOSYN-BIODILAERO'
#biom_id = 'BIOSYN-BIODILAERO-NOGAM'

# set to True to avoid printing .lst files and other large log files
hide_output = True
#doesn't work for GAMS 42 (on Roar cluster): output_redirect_str = ' writeOutput=0 ll=0' if hide_output else ''
#output_redirect_str = ' o=/dev/null ll=0 lo=0' if hide_output else ''
output_redirect_str = ' o=/dev/null' if hide_output else ''

spont_rxn_suffixes = ['SPONT', 'UNKNOWN']
rxns_to_ignore_for_kapps = [] # ignored when printing rxns during C1
