#BHMacros II
##1) Create one big ntuple! (Faster than TChains)
Use hadd to combine all the data ntuples. For example:
```
eosmount ~/eos
hadd allMyData.root ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015CreMiniAODv1_28Nov15/* ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015DreMiniAODv1_27Nov15/* ~/eos/cms/store/group/phys_exotica/BH_RunII/Data/BH_Ntuples_Run2015DpmptRecov4_27Nov15/*
```
Tested on lxplus.
##2) Make one big MET filtering list!
Visit the MET twiki:
https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Event_Lists
Download the JetHT MET filtering lists, put them in your home directory,
untar them, and then use cat to combine them, e.g.:
```
tar -xzvf JetHT_Nov14.tar.gz
cat JetHT_ecalscn1043093.txt JetHT_csc2015.txt > allMyMETfilteredEvents.tx
```

##3) Run the flatTuplizer in compiled mode!
The macro takes three arguments:
* the name of your input BH ntuple (e.g. the one you hadded in step 1)
* the name you want for your output file
* the name of your MET filtering list

Run it with the -q option (to close ROOT once the process is finished),
the -l option (to get rid of that silly splash), and + option at the end of the name of the macro so it
compiles (this makes it go faster.)

e.g.:
```
root -l -q 'BHflatTuplizer.cc+("allMyData.root","myBHflatTuple.root","allMyMETfilteredEvents.txt")'
```
