>> Order of scripts. If you don't run them exactly in the order as described below things won't work.

>> Software. You need R, R Studio, Jupyter Lab, and Python.

## workflow ##
01OTevo_mainscript.R - run until checkpoint "BLASTp/microsynteny analyses"

02BLASTp_evaluation.R - run the whole script

01OTevo_mainscript.R - pick up where you left AFTER "BLASTp/microsynteny results visualisation prep" and 
                       run "dN/dS values from Dumas et al., 2021 --- GenEvo tool" until checkpoint "FUMA analyses".

03FUMA_analyses.R - run the whole script

01OTevo_mainscript.R - pick up where you left after "FUMA analyses" and continue reading 
                       "Abagen toolbox expression data pre-processing". Then run --->

04abagen_analyses_updated.ipynb - run the whole script

05Abagen_DK_atlas_expression_analyses02.R - run the whole script

01OTevo_mainscript.R - pick up where you left AFTER "Visualization of AHBA cortical data with ggseg" and continue 
                       with running the code from the checkpoint "Differential stability analysis"

## end of workflow ## 








