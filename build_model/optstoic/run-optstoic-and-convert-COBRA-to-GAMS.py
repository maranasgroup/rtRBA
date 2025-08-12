import pandas as pd
import cobra
from collections import OrderedDict
from copy import deepcopy
import numpy as np
from scipy.stats import linregress
import os,sys
model_root_path = '../../'
sys.path.append(model_root_path+'pycore/')
from gsm_custom_functions import *

curr_dir = os.getcwd()
# model = cobra.io.load_json_model(curr_dir + '/../../optional-extras/rt (kim et al. 2020)/Supplementary File S2/Rt_IFO0880.json')
model = cobra.io.load_json_model('/Users/ejm6426/Downloads/iML1515.json')
target_product = 'lipopb_c' # make blank if none meant to be added
# check for any demand reactions for the target product
if target_product in [met.id for met in model.metabolites]:
	demand_rxn = None
	for rxn in model.reactions:
		if target_product in rxn.reactants and rxn.upper_bound > 0:
			demand_rxn = rxn
			break
		elif target_product in rxn.products and rxn.lower_bound < 0:
			demand_rxn = rxn
			break
	if demand_rxn is None:
		# create a demand reaction for the target product
		demand_rxn = cobra.Reaction('DM_' + target_product)
		demand_rxn.name = 'Demand reaction for ' + target_product
		demand_rxn.lower_bound = 0
		demand_rxn.upper_bound = 1000
		model.add_reactions([demand_rxn])
		demand_rxn.add_metabolites({target_product: -1})

convert_cobra_to_gams(model,curr_dir + '/data/',add_slashes=False,mets_filename='metabolites_modified.txt',rxns_filename='reactions_modified.txt',sij_filename='Sij_modified.txt',rxntype_filename='rxntype_modified.txt')