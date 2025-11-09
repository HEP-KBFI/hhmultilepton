law run cf.PlotVariables1D \
    --version prod1 \
    --producers default \
    --variables nmu \
    --categories ceormu \
    --datasets all_backgrounds \
    --view-cmd imgcat \
    --configs 24_v15_central \
    --workflow slurm \
    --parallel-jobs 200 \
    $1 
    #--version testnanov15_2024__ver4 \
    #--limit-dataset-files 1 \ # work-in-progress, limit the number of files to 1

# options: 

#   --configs: 
#        22preEE_v14_private, 22postEE_v14_private, 23preBPix_v14_private, 23postBPix_v14_private
#        22preEE_v12_central, 22postEE_v12_central, 23preBPix_v12_central, 23postBPix_v12_central, 24_v15_central

#   --processes:
#       all_data, all_signals, all_backgrounds,       
#       resonant, nonresonant, nonresonant_ggf, nonresonant_vbf
#       ggf_4v, ggf_4t, ggf_2t2v, vbf_4v, vbf_4t, vbf_2t2v
#       4v, 4t, 2t2v

#   --datasets:
#       all_data, all_backgrounds, all_signals
#       ttbar, single_top, dy, wjets, qcd, zz, single_higgs, vvv, others



