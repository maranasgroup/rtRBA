# *rt*RBA: *R. toruloides* Resource Balance Analysis model
This repository provides resources and model files for the genome-scale resource balance analysis (RBA) model *rt*RBA.
Eric J. Mooney and Costas D. Maranas<br>
Additional info is present in the publication for this model [^1].<br>
Program and packages requirements for the provided software implementation: Python 3 (plus packages: cobra, pandas, numpy, matplotlib, jupyter, scikit-learn, cplex (compiled from installed CPLEX linear programming solver), and all automatically-installed associated dependencies to those already listed) and GAMS with ("built-in") Soplex solver. (Tested on Python 3.6 and GAMS 39.1.0. Python package installation could be done through "pip" and could take up to 1 hour)
gsm_custom_functions.py contains many useful functions, and may be easier to use if you move it into your $PYTHONPATH directory.”<br><br>
To build this model, I modified and expanded upon code from another model - the paper below has more info about it:<br>
Evaluating proteome allocation of Saccharomyces cerevisiae phenotypes with resource balance analysis<br>
Hoang V. Dinh and Costas D. Maranas<br>
Metabolic Engineering, 2023, 77:242-255; doi: https://doi.org/10.1016/j.ymben.2023.04.009<br>
Formulation and software usage guide are available at suppMat/scRBA_suppMethods_2022-09-12.docx<br>
# Installation
To install directly via the environment.yml file, Anaconda must already be installed. To create a self-contained environment for RBA, activate the Anaconda prompt and type:
```conda env create -f <path to environment.yml file>```
# Building an RBA model
See "RBA User Guide.docx" under suppMat.
# Funding
This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Biological and Environmental Research Program under Award Number DE-SC0018420). Any opinions, findings, and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy. Computations for this research were performed on the Pennsylvania State University’s Roar Collab supercomputer. The authors of this work recognize the Penn State Institute for Computational and Data Sciences (RRID:SCR_025154) for providing access to computational research infrastructure within the Roar Core Facility (RRID: SCR_026424).
# References
[^1]: Eric J. Mooney, Patrick F. Suthers, Wheaton L. Schroeder, Hoang V. Dinh, Xi Li, Yihui Shen, Tianxia Xiao, Catherine M. Call, Heide Baron, Arjuna M. Subramanian, Daniel R. Weilandt, Felix C. Keber, Martin Wühr, Joshua D. Rabinowitz, Costas D. Maranas. "Metabolic flux and resource balance in the oleaginous yeast Rhodotorula toruloides." Submitted.