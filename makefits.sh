# DrawUncertainty = False, DrawRatioPanel=False
#python fitSThists.py QCD.root QCD_Unshaded.root FitNormRanges_2016.txt QCD_unshaded.json -b
# DrawUncertainty = True, DrawRatioPanel= True  Lower norm ranges to 2 TeV
#python fitSThists.py QCD.root QCD_Exc3ratio.root FitNormRanges_2016.txt useMET -b
# DrawUncertainty = True, DrawRatioPanel= False 
#python fitSThists.py QCD.root QCD_shaded.root FitNormRanges_2016.txt useMET -b


#Full Data plots
# DrawUncertainty = False, special range
#python fitSThists.py dataPreApprove.root fullData_preApprove_Unshaded.root  FitNormRanges_2016_fullRange.txt data_unshaded.json -b
#python fitSThists.py dataPreApprove.root fullData_preApprove.root  FitNormRanges_2016.txt useMET -b
# DrawUncertainty = True, DrawRatioPanel= True
#python fitSThists.py dataPreApprove.root fullData_preApprove_pulls.root  FitNormRanges_2016.txt data_pulls.json useMET -b
#python fitSThists.py dataPreApprove.root fullData_preApprove_pulls.root  FitNormRanges_2016.txt useMET -b
#python fitSThists.py dataPreApprove.root fullData_preApprove_signals.root  FitNormRanges_2016.txt useMET -b
#python fitSThists.py dataPreApprove.root fullData_allFunctions.root  FitNormRanges_2016.txt useMET -b
#python fitSThists.py dataPreApprove.root fullDataShaded.root  FitNormRanges_2016.txt useMET -b



# Update plots for approval
python fitSThists.py dataPreApprove.root fullData_preApprove.root  FitNormRanges_2016.txt data_pulls.json useMET -b
#python fitSThists.py dataPreApprove.root fullData_preApprove_BHsignals.root  FitNormRanges_2016_fullRange.txt data_BHsignals.json useMET -b
#python fitSThists.py dataPreApprove.root fullData_preApprove_SphaleronSignals.root  FitNormRanges_2016.txt data_SphaleronSignal.json useMET -b
