#BHMacros II
##1) Create one big ntuple! (Faster than TChains)
Use hadd to combine all the data ntuples. For example:
```
eosmount ~/eos
hadd allMyData.root ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015CreMiniAODv1_28Nov15/* ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015DreMiniAODv1_27Nov15/* ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015DpmptRecov4_27Nov15/*
```
##2) Make one big MET filtering list!
Visit the MET twiki:
https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Event_Lists
Download the JetHT MET filtering lists, put them in a convenient
directory, untar them, and then use cat to combine them, e.g.:
```
tar -xzvf JetHT_Nov14.tar.gz
cat JetHT_ecalscn1043093.txt JetHT_csc2015.txt > allMyMETfilteredEvents.txt
```

##3) Run the flatTuplizer in compiled mode!
The BHflatTuplizer takes three arguments:
* the name of your input BH ntuple (e.g. the one you hadded in step 1)
* the name you want for your output files
* the name of your MET filtering list

Run it with the -q option (to close ROOT once the process is finished),
the -l option (to get rid of that silly splash), and the + option at the 
end of the name of the macro to run in compiled mode (this makes it go faster.)

e.g.:
```
root -l -q 'BHflatTuplizer.cc+("allMyData.root","myBHflatTuple.root","allMyMETfilteredEvents.txt")'
```
##4) Define your fitting and normalization regions!
You will need to create a text file to define your fit ranges and it
should look something like this:
```
FitRange,  exc,  2,  1300, 3000
FitRange,  exc,  3,  1300, 3000
NormRange, inc,  2,  2000, 2300
NormRange, exc,  2,  2000, 2300
NormRange, inc,  3,  2000, 2300
```
Don't worry about whitespace, but make sure your commas are the same as
above.
##5) Fit the ST distributions and make tables for the combine tool!
fitSThists.py takes four arguments:
* the name of your input BH flatTuple (e.g. the one you made in step 3)
* the name you want for your output files
* the name of the textfile of your fit and normalization regions
* either "useMET" or "useMHT" to select ST distributions that
  incorporate either MET or MHT.

Run it with the -b option to suppress pyroot's frenzied flashing
plots.

e.g.:
```
python fitSThists.py myBHflatTuple.root myOutputFile.root myFitAndNormRanges.txt useMET -b
```
Your output histograms and textfiles for the combine tool will appear in
the output directory. i.e. 
```
mkdir ./output
```
