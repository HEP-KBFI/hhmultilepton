law run cf.PlotVariables1D \
    --version testnanov15_2024__ver4 \
    --producers default \
    --variables nmu \
    --datasets qcd_mu_pt15to20_pythia \
    --view-cmd imgcat \
    --configs 24_v15_central \
    $1 
    #--limit-dataset-files 1 \ # work-in-progress
    #--configs 23preBPix_v12_central \ 

# Useful tips:
# for testing purposes, limit the number of files to 1 using : --limit-dataset-files 1

