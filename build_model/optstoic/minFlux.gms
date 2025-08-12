$Ontext
minFlux

Authors: Anupam Chowdhury, Chiam Yu Ng

DESCRIPTION
The objective of this code (minFlux.gms) is to identify a pathway with
minimal total flux that conform to the given stoichiometry.

In this sample code, we use the reaction equation identified using optStoic for
acetate production from glucose:
                            glucose -> 3 acetate + 3 H+

OUTPUT FILE
minFlux_output.txt - a list of reactions constituting the pathway and the flux through
                     each reaction

For more information, please refer to the paper below:
Anupam Chowdhury and Costas D. Maranas
"Designing overall stoichiometric conversions and intervening metabolic reactions"
Scientific Reports, 2015
$Offtext

$INLINECOM /*  */
$set fpath data/
$set fsuffix _modified

options
    limrow = 5000
    optCR = 1E-6
    optCA = 1E-6
    iterlim = 100000
    decimals = 8
    reslim = 1200
    mip = cplex;

***************************************************
Sets

i   Metabolites
/
$include "%fpath%metabolites%fsuffix%.txt"
/

j   Reactions
/
$include "%fpath%reactions%fsuffix%.txt"
/

rm_rxn(j)   "reactions with potential stoichiometric imbalance"
$include "%fpath%inc_stoic.txt"

k           "iteration of minFlux/minRxn" /1*10/

blocked(j)  blocked reactions
$include "%fpath%blocked_reactions.txt"

atp_irr(j)  "atp-related reactions with defined forward directionality"
$include "%fpath%atp_irreversible.txt"

backward_irr_rxn(j)
$include "%fpath%backward_irreversible_rxn.txt"

;


***************************************************
Parameters

S(i,j) Stoichiometry matrix
/
$include "%fpath%Sij%fsuffix%.txt"
/

rxntype(j) Reaction Type
/
$include "%fpath%rxntype%fsuffix%.txt"
/

store(k,j)  matrix to store pathway identified

LB(j)       lower bound of flux for reaction j

UB(j)       upper bound of flux for reaction j

active(k)   indicate if the integer cut constraint should be active at iteration k
;


*Initialize all result to zero
store(k,j) = 0;
active(k) = 0;

****************************************************
Variables

zr       "objective function value for minRxn"
zf       "objective function value for minFlux"
;

Integer variables

v(j)    reaction flux

x(j)    absolute value of reaction flux
;

Binary Variables

y(j)    indicator of whether a reaction participate in the pathway
;

***************************************************
Equations

objRxn          objective function for minRxn
objFlux         objective function for minFlux
stoic           stoichiometric balance for metabolite i
con1            lower bound constraints for reaction flux
con2            upper bound constraints for reaction flux
con3            define x as absolute value of reaction flux
con4            define x as absolute value of reaction flux
con5            integer cut constraint
con6            "setting upper bound for summation of number of reaction in pathway (optional)"
;

objRxn..        zr =e= sum(j$(rxntype(j) ne 4), y(j));
objFlux..       zf =e= sum(j$(rxntype(j) ne 4), x(j));
stoic(i)..      sum(j, S(i,j) * v(j) ) =e= 0;
con1(j)..       v(j) =g= LB(j)*y(j);
con2(j)..       v(j) =l= UB(j)*y(j);
con3(j)..       x(j) =g= v(j);
con4(j)..       x(j) =g= -v(j);
con5(k)$(active(k) = 1)..   sum(j$(store(k,j) = 1), 1 - y(j)) =g= 1;
con6..          sum(j, y(j)) =l= 12;

**************************************************
Scalar
vmax            maximum flux value  /1000/;

**************************************************
*Updating the bounds for all the parameters and variables

*Changing reaction type for reactions in set atp_irr
*rxntype(j)$atp_irr(j) = 0;
*rxntype(j)$backward_irr_rxn(j) = 2;

*Setting and lower bound (LB) and upper bound (UB) of all reaction
*based on their reaction type (0 = forward irreversible;
*1 = reversible; 2 = backward irreversible; 4 or 5 = exchange)
LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = vmax;
*LB(j)$(rxntype(j) = 2) = -vmax;
*UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = vmax;
*LB(j)$(rxntype(j) = 1) = -vmax;
LB(j)$(rxntype(j) = 1) = 0;
UB(j)$(rxntype(j) = 1) = vmax;
*LB(j)$(rxntype(j) = 4) = 0;
*UB(j)$(rxntype(j) = 4) = 0;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = vmax;
LB(j)$(rxntype(j) = 5) = 0;
UB(j)$(rxntype(j) = 5) = vmax;

*The uptake and export flux of metabolites are set to the
*stoichiometric coefficients in the reaction equation.
*1 glucose -> 3 acetate + 3 H+
*LB('RXN-EX_glc__D_e_REV-SPONT')=1;
*UB('RXN-EX_glc__D_e_REV-SPONT')=1;

*LB('EX_accoa_m')=1;
*UB('EX_accoa_m')=1;
LB('RXN-DM_lipopb_m_FWD-SPONT')=1;
UB('RXN-DM_lipopb_m_FWD-SPONT')=1;

*LB('DM_accoa_c')=1;
*UB('DM_accoa_c')=1000;

*LB('RXN-EX_ppa_e_REV-SPONT')=1;
*UB('RXN-EX_ppa_e_REV-SPONT')=1;
*LB('DM_succ_m')=1;
*UB('DM_succ_m')=1;

*Turn off reactions that are undesirable
*LB(j)$rm_rxn(j) = 0;
*UB(j)$rm_rxn(j) = 0;
*v.fx(j)$rm_rxn(j) = 0;
*y.fx(j)$rm_rxn(j) = 0;

v.up(j) = UB(j);
v.lo(j) = LB(j);

x.lo(j) = 0;
x.up(j) = vmax;

*remove ad_hoc reactions
**************************************************
Model
minrxn
/
objRxn
Stoic
con1
con2
con5
*constraint 6 may speed up the time
con6
/

minflux
/
objFlux
Stoic
con1
con2
con3
con4
con5
/
;


minrxn.optfile = 1;
minrxn.holdfixed = 1;

minflux.optfile = 1;
minflux.holdfixed = 1;

Scalar
continue    continue while loop /1/
n           current iteration   /0/
count       reaction count      /0/
nstop       terminate loop after nstop-th iteration    /10/
;

********************************************************************************
*                                   minFlux                                    *
********************************************************************************
File file2 /minFlux_output.txt/;
Put file2;

Put "*** network of minimum total flux assuring 100% conversion"//;

While( continue = 1,
        Solve minflux Using mip Minimizing zf;
        n = n + 1;
        count = 0;

        /*Write the solution to output file*/
        put "iteration no: ", @30, n/;
        put "sum of fluxes: ", @30, zf.l/;
        put "modelstat: ", @30, minflux.modelstat/;

        loop(j$(v.l(j) ne 0),
            store(k, j)$(ord(k) = n) = 1;
            put "'"j.tl:0"'", system.tab, system.tab, system.tab, system.tab, v.l(j):0:8/;
            count = count + 1;
        );
        put "no. of rxns: ", @30, count/;
        put "***************"//;

        /*Turn on the integer cut constraint at the next iteration*/
        active(k)$(ord(k) = n) = 1;

        /*Stop the program if model status is not optimal or not an integer solution*/
        if(n gt nstop or (minflux.modelstat ne 1 and minflux.modelstat ne 8),
                continue = 0;
        );
);
putclose file2;

