>> Order of scripts: If you don't run the scripts exactly in the order as described below things won't work.

>> Software: You need R, R Studio, Jupyter Lab, and Python.

>> Set up: Download all files and put them in one folder (e.g., "OTpathway_evolution-main"). 
           CREATE AN R PROJECT IN THAT FOLDER. 
           Operate from within that R project, otherwise the relative paths will NOT work.

## start of workflow ##
01BLASTp_evaluation.R - Open this script first and run the whole script.
  Output: scaled_max_scores.xlsx & homology_boolean.xlsx (basis for creating sup_mat_01, sheet 1), OTpthwy_res02_vertebrate_threshold.xlsx (--> sup_mat_03)

02dnds_analyses.R - Run the whole script.
  Output: sup_mat_12b.pdf, dndsOT_tree_mq_phylo6.pdf (figure 3), sup_mat_12a.pdf, sup_mat_05.pdf

03FUMA_results_eval.R - Run the whole script.
  Output: FUMA_tissueEnrich02.pdf (figure 4 B-D)

#### IMPORTANT INFO ####
The next part of the script is using the toolbox 'abagen' (https://github.com/rmarkello/abagen), which is written in PYTHON.
Read https://abagen.readthedocs.io/en/stable/ carefully for optimal use.
Do not attempt to run in R, it will not work. Run in python shell in the terminal or in Jupyter Lab/Notebook.

!!! Previous installation and set up strictly necessary!!! 
To start jupyter lab, open command line/bash and enter "jupyter lab" or "jupyter-lab", then it will open automatically in your browser. 
########################

04abagen_analyses_pub.ipynb - Run the whole script.
  Output: AHBAdk_220513_9861/10021/12876/14380/15496/15697.csv, DS_values_abagen.csv

05Abagen_DK_atlas_expression.R - Run the whole script.
  Output: abagen_analysis01_ScC_combined02.pdf (figure 5), abagen_analysis02_cSc_03.pdf (figure 6)

06differential_stability.R - Run the whole script.
  Output: sup_mat_09.xlsx

## end of workflow ## 








