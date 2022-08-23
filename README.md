# HalocarbonBanks
Code and data to creating figures and tables in "Bayesian assessment of chlorofluorocarbon (CFC), hydrochlorofluorocarbon (HCFC) and halon banks suggest large reservoirs still present in old equipment."  Lickley, Reimann, Daniel, Flemming, and Solomon, (submitted to ACP March 2022)

Banks and emissions model based on methods developed
Lickley et al. (2020, 2021) Nature Communications. 

Gases for modeling: 

HCFC-22
HCFC-141b
HCFC-142b
Halon 1211
Halon 1301
CFC-11
CFC-12
CFC-113
CFC-114
CFC-115

Differences between this analysis and Lickley et al. (2021): 
 - lifetimes are fixed in present analysis (not inferred) 
 - gases are all modeled independently
 - feedstock values are included in the analysis for HCFC-22, CFC-113 and HCFC-142b


Code for generating figures in main text. 

Contained in Code Directory

AGAGE_published.mat contains global mean mole fractions from the Advanced Global Atmospheric Gas Experiment (AGAGE; https://agage.mit.edu, last
accessed: March 29, 2022). (Prinn et al., 2000, 2018)

Develop_priors.m: matlab script to that generates prior distributions for production, RF and DE for each of chemical by equipment type, following AFEAS breakdown. 

Master.m: matlab script that runs the BPE model for all chemicals and generates posterior estimates

PostProcess.m: matlab script that makes the figures after the prior and posterior samples have been generated for each scenario and gas. 

SPARC_lifetimes.mat: file containing SPARC atmospheric lifetimes for each gas. 


boundedline.m: a function used to create the uncertainty shading in figures. 

CFC11, CFC12, CFC113, CFC114, CFC115, HCFC22, HCFC141b, HCFC142b, halons folders: contains InputFolders.  Each Input folder for each gas contains AFEAS and UNEP production data. 

%%%%%%%%%%%%%%%%%%  INSTRUCTION FOR RUNNING THE BPE MODEL %%%%%%%%%%%%%%%%%%%%%%

Step 1:  Run Develop_priors.m script.  (takes ~20 minutes to run on a normal desktop computer)

Step 2:  Run Master.m  (takes ~20 minutes to run on a normal desktop computer)

Step 3: Run PostProcess.m (takes ~ 5 minutes to run)
	
References: 
Prinn, R. G., Weiss, R. F., Fraser, P. J., Simmonds, P. G., Cunnold, D. M., Alyea, F. N., O’Doherty, S., Salameh, P., Miller, B. R., Huang, J., Wang, R. H. J., Hartley, D. E., Harth, C., Steele, L. P., Sturrock, G., Midgley, P. M., and McCulloch, A.: A history of chemically and radiatively important gases in air deduced from ALE/GAGE/AGAGE, J. Geophys. Res.-Atmos., 105, 17751–17792, https://doi.org/10.1029/2000JD900141, 2000.

Prinn, R. G., Weiss, R. F., Arduini, J., Arnold, T., DeWitt, H. L., Fraser, P. J., Ganesan, A. L., Gasore, J., Harth, C. M., Hermansen, O., Kim, J., Krummel, P. B., Li, S., Loh, Z. M., Lunder, C. R., Maione, M., Manning, A. J., Miller, B. R., Mitrevski, B., Mühle, J., O’Doherty, S., Park, S., Reimann, S., Rigby, M.,
Saito, T., Salameh, P. K., Schmidt, R., Simmonds, P. G., Steele, L. P., Vollmer, M. K., Wang, R. H., Yao, B., Yokouchi, Y., Young, D., and Zhou, L.: History of chemically and radiatively important atmospheric gases from the Advanced Global Atmospheric Gases Experiment (AGAGE), Earth Syst. Sci. Data, 10, 9851018, https://doi.org/10.5194/essd-10-985-2018, 2018
