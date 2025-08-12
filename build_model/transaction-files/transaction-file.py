import pandas as pd
import cobra
from collections import OrderedDict
from copy import deepcopy
import os,sys
model_root_path = '../../'
sys.path.append(model_root_path)
from pycore.gsm_custom_functions import *

# Applies changes to models based on transaction files, to allow easier documentation of updates
# Metabolic model (COBRApy json)
modelpath = '../input/iRhtoC.json'
model = cobra.io.load_json_model(modelpath)

# if the input model is known, only need the transaction file path
transaction_files_and_model_versions = [
	{'file':'./transaction-file-iRhto.csv'},
	{'file':'./transaction-file-ClimToNlim.csv',
	'id':'iRhtoN',
	'name':'iRhtoN',
	'model_output':'../input/iRhtoN.json',
	'biomass_changes':{
		'reaction':'1.505405 13BDglucan_c + 0.39616 16BDglucan_en + 0.022241 G00006_c + 0.121529 alatrna_c + 0.000101 amet_c + 0.048014 argtrna_c + 0.057717 asntrna_c + 0.057717 asptrna_c + 183.672546 atp_c + 0.000101 btn_c + 0.001297 ca2_c + 0.035592 chtn_c + 0.000101 coa_c + 0.025956 ctp_c + 0.000115 cu2_c + 0.001741 cystrna_c + 0.015209 datp_c + 0.023738 dctp_c + 0.023738 dgtp_c + 0.000202 docosa_c + 0.015209 dttp_c + 0.511075 epist_c + 0.000101 fad_c + 0.00067 fe2_c + 0.096277 glntrna_c + 0.096277 glutrna_c + 0.110582 glytrna_c + 0.025571 gpianchorSC_r + 0.022743 gtp_c + 183.6717 h2o_c + 0.002439 hdca_c + 0.000101 hemeA_c + 0.024007 histrna_c + 0.073265 iletrna_c + 0.014554 ipc_g + 0.608732 k_c + 0.099636 leutrna_c + 0.003081 linoea_c + 0.001219 linolen_c + 0.081724 lystrna_c + 0.01418 mettrna_c + 0.064142 mg2_c + 9.5e-05 mn2_c + 0.000101 nad_c + 0.000101 nadp_c + 0.002072 ocdca_c + 0.006333 ocdcea_c + 0.071858 oglycanSC_g + 0.004254 pail_c + 0.007985 pc_c + 0.001929 pe_c + 0.046771 phetrna_c + 0.052492 protrna_c + 0.001034 ps_c + 0.000101 ribflv_c + 0.0663 sertrna_c + 0.031486 so4_c + 0.000101 spmd_c + 0.063231 tag_c + 0.000101 thf_c + 0.000101 thmpp_c + 0.069285 thrtrna_c + 0.008085 trptrna_c + 0.000807 ttcosa_c + 0.02438 tyrtrna_c + 0.013605 utp_c + 0.091177 valtrna_c + 0.001589 zn2_c --> 183.6717 adp_c + 183.6717 h_c + 183.565156 pi_c + 0.157009 ppi_c + 0.121529 trnaala_c + 0.048014 trnaarg_c + 0.057717 trnaasn_c + 0.057717 trnaasp_c + 0.001741 trnacys_c + 0.096277 trnagln_c + 0.096277 trnaglu_c + 0.110582 trnagly_c + 0.024007 trnahis_c + 0.073265 trnaile_c + 0.099636 trnaleu_c + 0.081724 trnalys_c + 0.01418 trnamet_c + 0.046771 trnaphe_c + 0.052492 trnapro_c + 0.0663 trnaser_c + 0.069285 trnathr_c + 0.008085 trnatrp_c + 0.02438 trnatyr_c + 0.091177 trnaval_c',
		'name':'Biomass in chemostat (continuous) under nitrogen (NH4) limitation'
		}
	},
	{'file':'./transaction-file-ClimToBatch.csv',
	'id':'iRhtoBatch-Rabinowitz',
	'name':'iRhtoBatch-Rabinowitz',
	'model_output':'../input/iRhtoBatch-Rabinowitz.json',
	'biomass_changes':{
		'name':'Biomass in batch cultivation, modified to fit kapp calculation results for RBA model',
		'reaction':'0.667059797182664 13BDglucan_c + 0.175541846742192 16BDglucan_en + 0.432822757129132 alatrna_c + 0.242142743774166 argtrna_c + 0.136009600298269 asntrna_c + 0.229695228907398 asptrna_c + 251.299438 atp_c + 0.00125314645456112 ca2_c + 0.0157705460502466 chtn_c + 0.116254182209833 ctp_c + 0.000111087632830457 cu2_c + 0.0394116891383016 cystrna_c + 0.00585548759366856 datp_c + 0.009139393932604 dctp_c + 0.009139393932604 dgtp_c + 1.35753472356956e-05 docoscoa_c + 0.00585548759366856 dttp_c + 0.0417718733030105 ergst_c + 0.0103623323466586 ergstest_c + 0.000647036738591434 fe2_c + 0.158742150812667 glntrna_c + 0.275386918489312 glutrna_c + 0.228812264496074 glycogen_c + 0.331296909013928 glytrna_c + 0.101863461493711 gtp_c + 247.144701 h2o_c + 0.0854918252602606 histrna_c + 0.000974452919565412 hxccoa_c + 0.197690946911468 iletrna_c + 0.588228504895661 k_c + 0.370077561546643 leutrna_c + 0.000257507322625944 linocoa_c + 6.34327105074009e-05 linolncoa_c + 0.26940006845879 lystrna_c + 0.477628098524987 mannan_c + 0.0849858749385989 mettrna_c + 0.0619820268547972 mg2_c + 9.15985744391488e-05 mn2_c + 0.00292335875869624 nadph_c + 0.000974452919565416 o2_c + 0.000149965231871062 odecoa_c + 0.00948922253072798 pail_c + 0.0280476783838513 pc_c + 0.00752277653904498 pe_c + 0.155477321197495 phetrna_c + 0.00207155457628959 pmtcoa_c + 0.229924447191349 protrna_c + 0.00647329074467303 ps_c + 0.293336514123922 sertrna_c + 0.0304253435075909 so4_c + 8.16736330810605e-05 stcoa_c + 0.0102259089379194 tag_c + 0.237686677320607 thrtrna_c + 0.134188988194594 tre_c + 0.042649375902274 trptrna_c + 0.00100966845632259 ttccoa_c + 0.116922264264054 tyrtrna_c + 0.0609393622307783 utp_c + 0.313653417351731 valtrna_c + 0.00153573780123509 zn2_c --> 251.232837 adp_c + 0.00194890583913083 co2_c + 0.00462183019749875 coa_c + 0.00194890583913082 dag_c + 251.234556 h_c + 0.00292335875869624 nadp_c + 251.129881 pi_c + 0.375646728227392 ppi_c + 0.432822757129132 trnaala_c + 0.242142743774166 trnaarg_c + 0.136009600298269 trnaasn_c + 0.229695228907398 trnaasp_c + 0.0394116891383016 trnacys_c + 0.158742150812667 trnagln_c + 0.275386918489312 trnaglu_c + 0.331296909013928 trnagly_c + 0.0854918252602606 trnahis_c + 0.197690946911468 trnaile_c + 0.370077561546643 trnaleu_c + 0.26940006845879 trnalys_c + 0.0849858749385989 trnamet_c + 0.155477321197495 trnaphe_c + 0.229924447191349 trnapro_c + 0.293336514123922 trnaser_c + 0.237686677320607 trnathr_c + 0.042649375902274 trnatrp_c + 0.116922264264054 trnatyr_c + 0.313653417351731 trnaval_c'
		}
	},
]
for item in transaction_files_and_model_versions:
	with model as m:
		transaction_file = item['file']
		model_output = item['model_output'] if 'model_output' in item else modelpath
		biomass_changes = item['biomass_changes'] if 'biomass_changes' in item else {}
		if 'id' in item:
			m.id = item['id']
		if 'name' in item:
			m.name = item['name']
		try:
			df_cmds = pd.read_csv(transaction_file)
			df_cmds.index = df_cmds.id.to_list()
			model_copy = deepcopy(m)
			execute_command(m, model_copy, df_cmds,verbose=True)
		except Exception as e:
			print(e)
		for attr in biomass_changes:
			setattr(m.reactions.BIOMASS, attr, biomass_changes[attr])
		cobra.io.save_json_model(m, model_output)