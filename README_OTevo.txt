>> Order of scripts: If you don't run the scripts exactly in the order as described below things won't work.

>> Software: You need R, R Studio, Jupyter Lab, and Python.

>> Set up: Download all files and put them in one folder (e.g., "OTpathway_evolution-main"). Create an R Project in that folder. Operate from within 
           that R project, otherwise the relative paths will NOT work.

## start of workflow ##
01BLASTp_evaluation.R - Open this script first and run the whole script. 

02dnds_analyses.R - Run the whole script.

03FUMA_results_eval.R - Run the whole script.

#### IMPORTANT INFO ####
The next part of the script is using the toolbox 'abagen' (https://github.com/rmarkello/abagen), which is written in PYTHON.
Read https://abagen.readthedocs.io/en/stable/ carefully for optimal use.
Do not attempt to run in R, it will not work. Run in python shell in the terminal or in Jupyter Lab/Notebook.

!!! Previous installation and set up strictly necessary!!! 
To start jupyter lab, open command line/bash and enter "jupyter lab" or "jupyter-lab", then it will open automatically in your browser. 
########################

04abagen_analyses_pub.ipynb - Run the whole script.

05Abagen_DK_atlas_expression.R - Run the whole script.

06differential_stability.R - Run the whole script.

## end of workflow ## 








