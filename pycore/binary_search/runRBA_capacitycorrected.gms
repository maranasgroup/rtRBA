************************* Run RBA model ********************
*       Authors: (v1) Hoang Dinh, (v2) Eric Mooney
************************************************************

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"
$setGlobal nscale 1e5

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1

* Ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0.45

options
	LP = cplex /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 600 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
prosyn(j)
$include "%prosyn_path%"
prowaste(j)
$include "%prowaste_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_biomass(j)
$include "%biomass_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
kapp(j)
$include "%kapp_path%"
;

Variables
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e3 * %nscale%;
* bounds from GSM model
$include %model_root_path%GAMS/model/GSM_rxn_bounds.txt

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0; v.up('PROWASTE-TOTALPROTEIN') = inf; v.up('PROWASTE-PROTCYT') = inf;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;
* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e3 * %nscale%;

$include "%phenotype_path%"
* Set all organism-specific and/or condition-specific parameters in phenotype.txt (e.g., NGAM and biomass composition)

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, ModelProtAlloCorrection
$include %enz_cap_declares_path%
;

Obj..			z =e= v('BIOSYN-PROTMODELED');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
NonModelProtAllo..	v('BIOSYN-PROTMODELED') + v('BIOSYN-PROTDUMMY2') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
ModelProtAlloCorrection..	v('BIOSYN-PROTMODELED') =l= %corrected_protein_capacity_percentage% * (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
$include %enz_cap_eqns_path%

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, ModelProtAlloCorrection
$include %enz_cap_declares_path%
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff /runRBA.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /runRBA.flux.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j) / %nscale%):0:15/;
	);
);
putclose ff2;

