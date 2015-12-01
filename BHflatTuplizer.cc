#include <stdio.h>
#include <iostream>
#include "Riostream.h"
#include "TBranch.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <TROOT.h>
#include <TMath.h>

// Macro to make a BH flat ntuple. John Hakala 12/1/2015

void BHflatTuplizer(std::string inFilename, std::string outFilename, std::string metListFilename);
float dR(float eta1, float phi1, float eta2, float phi2);
std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName);

void BHflatTuplizer(std::string inFilename, std::string outFilename, std::string metListFilename) {
  std::map<unsigned, std::set<unsigned> > list = readEventList(metListFilename.c_str());
  bool debugFlag     = false ;
  int  eventsToDump  = 10    ;  // if debugFlag is true, then stop once the number of dumped events reaches eventsToDump
  bool dumpBigEvents = true  ;
  bool dumpIsoInfo   = true  ;
  int  nDumpedEvents = 0     ;

  // define output textfile
  ofstream outTextFile;
  std::string outTextFilename  = outFilename+"_log.tx";
  outTextFile.open(outTextFilename.c_str());

  // define output histograms
  TH2F METvsMHT                = TH2F("METvsMHT"                ,  "METvsMHT"                ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2            = TH2F("METvsMHTinc2"            ,  "METvsMHTinc2"            ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasMuon     = TH2F("METvsMHTinc2hasMuon"     ,  "METvsMHTinc2hasMuon"     ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasPhoton   = TH2F("METvsMHTinc2hasPhoton"   ,  "METvsMHTinc2hasPhoton"   ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2hasElectron = TH2F("METvsMHTinc2hasElectron" ,  "METvsMHTinc2hasElectron" ,  1000,  0.,  20000.,  1000,  0.,  20000.);
  TH2F METvsMHTinc2onlyJets    = TH2F("METvsMHTinc2onlyJets"    ,  "METvsMHTinc2onlyJets"    ,  1000,  0.,  20000.,  1000,  0.,  20000.);

  TH1F MuonJetIso1             = TH1F("MuonJetIso1"             ,  "MuonJetIso1"             ,  300,   0.,  3);
  TH1F MuonJetIso2             = TH1F("MuonJetIso2"             ,  "MuonJetIso2"             ,  300,   0.,  3);
  TH1F MuonJetIso3             = TH1F("MuonJetIso3"             ,  "MuonJetIso3"             ,  300,   0.,  3);
  TH1F MuonJetIso4             = TH1F("MuonJetIso4"             ,  "MuonJetIso4"             ,  300,   0.,  3);
  TH1F MuonJetoverlapdR1       = TH1F("MuonJetIso4overlapdR1"   ,  "MuonJetIso4overlapdR1"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR2       = TH1F("MuonJetIso4overlapdR2"   ,  "MuonJetIso4overlapdR2"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR3       = TH1F("MuonJetIso4overlapdR3"   ,  "MuonJetIso4overlapdR3"   ,  300,   0., .3);
  TH1F MuonJetoverlapdR4       = TH1F("MuonJetIso4overlapdR4"   ,  "MuonJetIso4overlapdR4"   ,  300,   0., .3);
  TH1F ElectronJetIso1         = TH1F("ElectronJetIso1"         ,  "ElectronJetIso1"         ,  300,   0.,  3);
  TH1F ElectronJetIso2         = TH1F("ElectronJetIso2"         ,  "ElectronJetIso2"         ,  300,   0.,  3);
  TH1F ElectronJetIso3         = TH1F("ElectronJetIso3"         ,  "ElectronJetIso3"         ,  300,   0.,  3);
  TH1F ElectronJetIso4         = TH1F("ElectronJetIso4"         ,  "ElectronJetIso4"         ,  300,   0.,  3);
  TH1F ElectronJetoverlapdR1   = TH1F("ElectronJetoverlapdR1"   ,  "ElectronJetoverlapdR1"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR2   = TH1F("ElectronJetoverlapdR2"   ,  "ElectronJetoverlapdR2"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR3   = TH1F("ElectronJetoverlapdR3"   ,  "ElectronJetoverlapdR3"   ,  300,   0., .3);
  TH1F ElectronJetoverlapdR4   = TH1F("ElectronJetoverlapdR4"   ,  "ElectronJetoverlapdR4"   ,  300,   0., .3);
  TH1F PhotonJetIso1           = TH1F("PhotonJetIso1"           ,  "PhotonJetIso1"           ,  300,   0.,  3);
  TH1F PhotonJetIso2           = TH1F("PhotonJetIso2"           ,  "PhotonJetIso2"           ,  300,   0.,  3);
  TH1F PhotonJetIso3           = TH1F("PhotonJetIso3"           ,  "PhotonJetIso3"           ,  300,   0.,  3);
  TH1F PhotonJetIso4           = TH1F("PhotonJetIso4"           ,  "PhotonJetIso4"           ,  300,   0.,  3);
  TH1F PhotonJetoverlapdR1     = TH1F("PhotonJetoverlapdR1"     ,  "PhotonJetoverlapdR1"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR2     = TH1F("PhotonJetoverlapdR2"     ,  "PhotonJetoverlapdR2"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR3     = TH1F("PhotonJetoverlapdR3"     ,  "PhotonJetoverlapdR3"     ,  300,   0., .3);
  TH1F PhotonJetoverlapdR4     = TH1F("PhotonJetoverlapdR4"     ,  "PhotonJetoverlapdR4"     ,  300,   0., .3);

  // loop to create ST histograms for inclusive and exclusive multiplicities from 2 up to multMax
  TH1F stHist = TH1F("stHist", "ST", 100, 500, 10500);
  int mult=2;
  int multMax = 12;
  TH1F *stIncHist[multMax-2];
  TH1F *stExcHist[multMax-2];
  TH1F stHistMHT = TH1F("stHistMHT", "ST using MHT", 100, 500, 10500);
  TH1F *stIncHistMHT[multMax-2];
  TH1F *stExcHistMHT[multMax-2];
  char *histTitle = new char[20];
  // These use pat::slimmedMETs
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHist", mult);
    stIncHist[iHist] = new TH1F(histTitle, "Inclusive ST", 100, 500, 10500);
    sprintf(histTitle, "stExc%02dHist", mult);
    stExcHist[iHist] = new TH1F(histTitle, "Exclusive ST", 100, 500, 10500);
    ++mult;
  }
  mult=2;
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    sprintf(histTitle, "stInc%02dHistMHT", mult);
    stIncHistMHT[iHist] = new TH1F(histTitle, "Inclusive ST using MHT", 100, 500, 10500);
    sprintf(histTitle, "stExc%02dHistMHT", mult);
    stExcHistMHT[iHist] = new TH1F(histTitle, "Exclusive ST using MHT", 100, 500, 10500);
    ++mult;
  }

  // variables calculated in the loop
  float OurMet          = 0.            ;
  float Px              = 0.            ;
  float Py              = 0.            ;
  float ST              = 0.            ;
  float STMHTnoMET      = 0.            ;
  int multiplicity      = 0             ;
  bool passIso          = true          ;
  char *messageBuffer   = new char[400] ;
  bool eventHasMuon     = false         ;
  bool eventHasPhoton   = false         ;
  bool eventHasElectron = false         ;
  float JetMuonEt       = 0.            ;
  float JetElectronEt   = 0.            ;
  float JetPhotonEt     = 0.            ;

  // variables accessed from the tree
  Bool_t     firedHLT_PFHT800_v2       ;
  Bool_t     passed_CSCTightHaloFilter ;
  Bool_t     passed_goodVertices       ;
  Bool_t     passed_eeBadScFilter      ;
  int        runno                     ;
  long long  evtno                     ;
  int        lumiblock                 ;
  float      JetEt [25]                ;
  float      JetPx [25]                ;
  float      JetPy [25]                ;
  float      JetEta[25]                ;
  float      JetPhi[25]                ;
  float      JetNeutHadFrac[25]        ;
  float      JetNeutEMFrac[25]         ;
  float      JetChgHadFrac[25]         ;
  float      JetMuFrac[25]             ;
  float      JetChgEMFrac[25]          ;
  int        JetNConstituents[25]      ;
  int        JetNNeutConstituents[25]  ;
  int        JetNChgConstituents[25]   ;
  float      EleEt[25]                 ;
  float      ElePx[25]                 ;
  float      ElePy[25]                 ;
  float      EleEta[25]                ;
  float      ElePhi[25]                ;
  float      PhEt[25]                  ;
  float      PhPx[25]                  ;
  float      PhPy[25]                  ;
  float      PhEta[25]                 ;
  float      PhPhi[25]                 ;
  float      MuEt[25]                  ;
  float      MuPx[25]                  ;
  float      MuPy[25]                  ;
  float      MuEta[25]                 ;
  float      MuPhi[25]                 ;
  float      MuPFdBiso[25]             ;
  float      Met                       ;

  // tree branches
  TBranch  *b_firedHLT_PFHT800_v2       ;
  TBranch  *b_passed_CSCTightHaloFilter ;
  TBranch  *b_passed_goodVertices       ;
  TBranch  *b_passed_eeBadScFilter      ;
  TBranch  *b_JetEt                     ;
  TBranch  *b_JetPx                     ;
  TBranch  *b_JetPy                     ;
  TBranch  *b_JetEta                    ;
  TBranch  *b_JetPhi                    ;
  TBranch  *b_JetNeutHadFrac            ;
  TBranch  *b_JetNeutEMFrac             ;
  TBranch  *b_JetChgHadFrac             ;
  TBranch  *b_JetMuFrac                 ;
  TBranch  *b_JetChgEMFrac              ;
  TBranch  *b_JetNConstituents          ;
  TBranch  *b_JetNNeutConstituents      ;
  TBranch  *b_JetNChgConstituents       ;
  TBranch  *b_EleEt                     ;
  TBranch  *b_ElePx                     ;
  TBranch  *b_ElePy                     ;
  TBranch  *b_EleEta                    ;
  TBranch  *b_ElePhi                    ;
  TBranch  *b_PhEt                      ;
  TBranch  *b_PhPx                      ;
  TBranch  *b_PhPy                      ;
  TBranch  *b_PhEta                     ;
  TBranch  *b_PhPhi                     ;
  TBranch  *b_MuEt                      ;
  TBranch  *b_MuPx                      ;
  TBranch  *b_MuPy                      ;
  TBranch  *b_MuEta                     ;
  TBranch  *b_MuPhi                     ;
  TBranch  *b_MuPFdBiso                 ;
  TBranch  *b_Met                       ;
  TBranch  *b_runno                     ;
  TBranch  *b_evtno                     ;
  TBranch  *b_lumiblock                 ;

  //create a chain by looping over the input filename
  TChain chain("bhana/t");
  //ifstream infile;
  //infile.open(inFilename.c_str());
  //std::string buffer;
  //const char *eosURL = "root://eoscms.cern.ch/";
  //chain.SetMakeClass(1);
  //while (std::getline(infile, buffer)) {
  //  std::string ntupleURL = eosURL + buffer;
  chain.Add(inFilename.c_str());
  //}

  cout << "Opened chain: " << chain.GetName() << endl;

  // set all branch addresses
  chain.SetBranchAddress( "firedHLT_PFHT800_v2"       ,  &firedHLT_PFHT800_v2       ,  &b_firedHLT_PFHT800_v2       );
  chain.SetBranchAddress( "passed_CSCTightHaloFilter" ,  &passed_CSCTightHaloFilter ,  &b_passed_CSCTightHaloFilter );
  chain.SetBranchAddress( "passed_goodVertices"       ,  &passed_goodVertices       ,  &b_passed_goodVertices       );
  chain.SetBranchAddress( "passed_eeBadScFilter"      ,  &passed_eeBadScFilter      ,  &b_passed_eeBadScFilter      );
  chain.SetBranchAddress( "runno"                     ,  &runno                     ,  &b_runno                     );
  chain.SetBranchAddress( "lumiblock"                 ,  &lumiblock                 ,  &b_lumiblock                 );
  chain.SetBranchAddress( "evtno"                     ,  &evtno                     ,  &b_evtno                     );
  chain.SetBranchAddress( "JetEt"                     ,  JetEt                      ,  &b_JetEt                     );
  chain.SetBranchAddress( "JetPx"                     ,  JetPx                      ,  &b_JetPx                     );
  chain.SetBranchAddress( "JetPy"                     ,  JetPy                      ,  &b_JetPy                     );
  chain.SetBranchAddress( "JetEta"                    ,  JetEta                     ,  &b_JetEta                    );
  chain.SetBranchAddress( "JetPhi"                    ,  JetPhi                     ,  &b_JetPhi                    );
  chain.SetBranchAddress( "JetNeutHadFrac"            ,  JetNeutHadFrac             ,  &b_JetNeutHadFrac            );
  chain.SetBranchAddress( "JetNeutEMFrac"             ,  JetNeutEMFrac              ,  &b_JetNeutEMFrac             );
  chain.SetBranchAddress( "JetChgHadFrac"             ,  JetChgHadFrac              ,  &b_JetChgHadFrac             );
  chain.SetBranchAddress( "JetMuFrac"                 ,  JetMuFrac                  ,  &b_JetMuFrac                 );
  chain.SetBranchAddress( "JetChgEMFrac"              ,  JetChgEMFrac               ,  &b_JetChgEMFrac              );
  chain.SetBranchAddress( "JetNConstituents"          ,  JetNConstituents           ,  &b_JetNConstituents          );
  chain.SetBranchAddress( "JetNNeutConstituents"      ,  JetNNeutConstituents       ,  &b_JetNNeutConstituents      );
  chain.SetBranchAddress( "JetNChgConstituents"       ,  JetNChgConstituents        ,  &b_JetNChgConstituents       );
  chain.SetBranchAddress( "EleEt"                     ,  EleEt                      ,  &b_EleEt                     );
  chain.SetBranchAddress( "ElePx"                     ,  ElePx                      ,  &b_ElePx                     );
  chain.SetBranchAddress( "ElePy"                     ,  ElePy                      ,  &b_ElePy                     );
  chain.SetBranchAddress( "EleEta"                    ,  EleEta                     ,  &b_EleEta                    );
  chain.SetBranchAddress( "ElePhi"                    ,  ElePhi                     ,  &b_ElePhi                    );
  chain.SetBranchAddress( "PhEt"                      ,  PhEt                       ,  &b_PhEt                      );
  chain.SetBranchAddress( "PhPx"                      ,  PhPx                       ,  &b_PhPx                      );
  chain.SetBranchAddress( "PhPy"                      ,  PhPy                       ,  &b_PhPy                      );
  chain.SetBranchAddress( "PhEta"                     ,  PhEta                      ,  &b_PhEta                     );
  chain.SetBranchAddress( "PhPhi"                     ,  PhPhi                      ,  &b_PhPhi                     );
  chain.SetBranchAddress( "MuEt"                      ,  MuEt                       ,  &b_MuEt                      );
  chain.SetBranchAddress( "MuPx"                      ,  MuPx                       ,  &b_MuPx                      );
  chain.SetBranchAddress( "MuPy"                      ,  MuPy                       ,  &b_MuPy                      );
  chain.SetBranchAddress( "MuEta"                     ,  MuEta                      ,  &b_MuEta                     );
  chain.SetBranchAddress( "MuPhi"                     ,  MuPhi                      ,  &b_MuPhi                     );
  chain.SetBranchAddress( "MuPFdBiso"                 , MuPFdBiso                   ,  &b_MuPFdBiso                 );
  chain.SetBranchAddress( "Met"                       ,  &Met                       ,  &b_Met                       );

  const int nEvents = chain.GetEntries();
  cout << "Number of events in chain is: " << nEvents << endl;

  // loop over all events
  char percent[5];
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (iEvent%50000==0) {
      sprintf(percent, "%i", (100*iEvent)/nEvents);
      cout << percent << "\% done: Scanned " << iEvent << " events." << endl;
    }

    // reset variables
    OurMet           = 0.    ;
    Px               = 0.    ;
    Py               = 0.    ;
    ST               = 0.    ;
    multiplicity     = 0     ;
    eventHasMuon     = false ;
    eventHasPhoton   = false ;
    eventHasElectron = false ;

    chain.GetEntry(iEvent);
    // apply trigger and filter requirements
    if (    !firedHLT_PFHT800_v2 || !passed_CSCTightHaloFilter
         || !passed_goodVertices || !passed_eeBadScFilter      ) continue;

    // use Yutaro's method for applying the event filter
    auto rItr(list.find(runno));
    if (rItr != list.end()) {
      if (rItr->second.find(evtno) != rItr->second.end()){
        if (dumpBigEvents) {
          sprintf(messageBuffer, "Event in MET list skipped: run number %d lumi section %d event number %lld\n", runno, lumiblock, evtno);
          outTextFile << messageBuffer;
        }
        continue;
      }
    }
    // apply isolation requirement and calculate ST and MHT.
    //Jets
    for (int iJet = 0; iJet < 25; ++iJet) {
      passIso=true;
      JetMuonEt     =0;
      JetElectronEt =0;
      JetPhotonEt   =0;
      if (JetEt[iJet]>50.) {
        for (int iMuon = 0; iMuon < 25; ++iMuon ) {
          if (MuEt[iMuon]>50 && MuPFdBiso[iMuon]<0.15) {
            eventHasMuon = true;
            if (JetEt[iJet] && dR(JetEta[iJet],JetPhi[iJet], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
              JetMuonEt+=MuEt[iMuon];
              if (MuEt[iMuon]<150) {
                MuonJetIso1.Fill(MuEt[iMuon]/JetEt[iJet]);
                MuonJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
              }
              if (150<=MuEt[iMuon] && MuEt[iMuon]<250) {
                MuonJetIso2.Fill(MuEt[iMuon]/JetEt[iJet]);
                MuonJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
              }
              if (250<=MuEt[iMuon] && MuEt[iMuon]<400) {
                MuonJetIso3.Fill(MuEt[iMuon]/JetEt[iJet]);
                MuonJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
              }
              if (400<=MuEt[iMuon]) {
                MuonJetIso4.Fill(MuEt[iMuon]/JetEt[iJet]);
                MuonJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],MuEta[iMuon],MuPhi[iMuon]));
              }
              if (JetMuonEt>0.8*JetEt[iJet]) {
                passIso = false;
                if (dumpIsoInfo) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iJet, iMuon, runno, lumiblock, evtno);
                  outTextFile << messageBuffer;
                }
                break;
              }
            }
          }
        }
        for (int iElectron = 0; iElectron < 25; ++iElectron ) {
          if (EleEt[iElectron]>50) {
            eventHasElectron = true;
            if (dR(JetEta[iJet],JetPhi[iJet], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
              JetElectronEt+=EleEt[iElectron];
              if (EleEt[iElectron]<150) {
                ElectronJetIso1.Fill(EleEt[iElectron]/JetEt[iJet]);
                ElectronJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
              }
              if (150<=EleEt[iElectron] && EleEt[iElectron]<250) {
                ElectronJetIso2.Fill(EleEt[iElectron]/JetEt[iJet]);
                ElectronJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
              }
              if (250<=EleEt[iElectron] && EleEt[iElectron]<400) {
                ElectronJetIso3.Fill(EleEt[iElectron]/JetEt[iJet]);
                ElectronJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
              }
              if (400<=EleEt[iElectron]) {
                ElectronJetIso4.Fill(EleEt[iElectron]/JetEt[iJet]);
                ElectronJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],EleEta[iElectron],ElePhi[iElectron]));
              }
              if (JetElectronEt > 0.7*JetEt[iJet] ) {
                passIso = false;
                if (dumpIsoInfo) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %lld\n", iJet, iElectron, runno, lumiblock, evtno);
                  outTextFile << messageBuffer;
                }
                break;
              }
            }
          }
        }
        for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          if (PhEt[iPhoton]>50) {
            eventHasPhoton = true;
            if (dR(JetEta[iJet],JetPhi[iJet], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
              JetPhotonEt+=PhEt[iPhoton];
              if (PhEt[iPhoton]<150) {
                PhotonJetIso1.Fill(PhEt[iPhoton]/JetEt[iJet]);
                PhotonJetoverlapdR1.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
              }
              if (150<=PhEt[iPhoton] && PhEt[iPhoton]<250) {
                PhotonJetIso2.Fill(PhEt[iPhoton]/JetEt[iJet]);
                PhotonJetoverlapdR2.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
              }
              if (250<=PhEt[iPhoton] && PhEt[iPhoton]<400) {
                PhotonJetIso3.Fill(PhEt[iPhoton]/JetEt[iJet]);
                PhotonJetoverlapdR3.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
              }
              if (400<=PhEt[iPhoton]) {
                PhotonJetIso4.Fill(PhEt[iPhoton]/JetEt[iJet]);
                PhotonJetoverlapdR4.Fill(dR(JetEta[iJet],JetPhi[iJet],PhEta[iPhoton],PhPhi[iPhoton]));
              }
              if (JetPhotonEt>0.5*JetEt[iJet] ) {
                passIso = false;
                if (dumpIsoInfo) {
                  sprintf(messageBuffer, "Jet number %d failed isolation with Photon number %d  in run number %d lumi section %d event number %lld\n", iJet, iPhoton, runno, lumiblock, evtno);
                  outTextFile << messageBuffer;
                }
                break;
              }
            }
          }
        }
        if (!passIso) continue;

        if (debugFlag) outTextFile << "    JetEt for jet number " << iJet << " is: " << JetEt[iJet] << endl;
        ST += JetEt[iJet];
        multiplicity+=1;
        if (debugFlag && dumpIsoInfo) {
          sprintf(messageBuffer, "Jet number %d passed isolation in run number %d lumi section %d event number %lld.\n       It had Px=%f and Py=%f\n", iJet, runno, lumiblock, evtno, JetPx[iJet], JetPy[iJet]);
          outTextFile << messageBuffer;
        }
        Px += JetPx[iJet];
        Py += JetPy[iJet];
        if (debugFlag && dumpIsoInfo) {
          sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
          outTextFile << messageBuffer;
        }
      }
      else break;
    }

    //Electrons
    if (eventHasElectron) {
      for (int iElectron = 0; iElectron < 25; ++iElectron) {
        passIso=true;
        if (EleEt[iElectron]>50.) {
          for (int iJet = 0; iJet < 25; ++iJet ) {
            if (JetEt[iJet]>50 && dR(EleEta[iElectron],ElePhi[iElectron], JetEta[iJet], JetPhi[iJet]) < 0.3) {
              if (EleEt[iElectron]<0.7*JetEt[iJet]) {
                passIso = false;
                if (dumpIsoInfo) {
                  sprintf(messageBuffer, "Electron number %d failed isolation with Jet number %d  in run number %d lumi section %d event number %lld\n", iElectron, iJet, runno, lumiblock, evtno);
                  outTextFile << messageBuffer;
                }
                break;
              }
            }
          }
          if (!passIso) continue;

          // Throw away photon if there is an electron/photon overlap
          //for (int iPhoton = 0; iPhoton < 25; ++iPhoton ) {
          //  if (PhEt[iPhoton]>50 && dR(EleEta[iElectron],ElePhi[iElectron], PhEta[iPhoton], PhPhi[iPhoton]) < 0.3) {
          //    if (dumpIsoInfo) {
          //      sprintf(messageBuffer, "Electron number %d failed isolation with Photon number %d  in run number %d lumi section %d event number %lld\n", iElectron, iPhoton, runno, lumiblock, evtno);
          //      outTextFile << messageBuffer;
          //    }
          //    passIso = false;
          //    break;
          //  }
          //}
          //if (!passIso) continue;

          // Throw away electron if there's an electron/muon overlap.
          for (int iMuon = 0; iMuon < 25; ++iMuon ) {
            if (MuEt[iMuon]>50 && MuPFdBiso[iMuon]<0.15 && dR(EleEta[iElectron],ElePhi[iElectron], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
              passIso = false;
              if (dumpIsoInfo) {
                sprintf(messageBuffer, "Electron number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iElectron, iMuon, runno, lumiblock, evtno);
                outTextFile << messageBuffer;
              }
              break;
            }
          }
          if (!passIso) continue;

          if (debugFlag) cout << "    EleEt for electron number " << iElectron << " is: " << EleEt[iElectron] << endl;
          ST += EleEt[iElectron];
          multiplicity+=1;
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "Ele number %d passed isolation in run number %d lumi section %d event number %lld.      \n It had Px=%f and Py=%f\n", iElectron, runno, lumiblock, evtno, ElePx[iElectron], ElePy[iElectron]);
            outTextFile << messageBuffer;
          }
          Px += ElePx[iElectron];
          Py += ElePy[iElectron];
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
            outTextFile << messageBuffer;
          }
        }
        else break;
      }
    }

    //Photons
    if (eventHasPhoton) {
      for (int iPhoton = 0; iPhoton < 25; ++iPhoton) {
        passIso=true;
        if (PhEt[iPhoton]>50.) {
          for (int iJet = 0; iJet < 25; ++iJet ) {
            if (JetEt[iJet]>50 && dR(PhEta[iPhoton],PhPhi[iPhoton], JetEta[iJet], JetPhi[iJet]) < 0.3) {
              if (PhEt[iPhoton]<0.5*JetEt[iJet]) {
                passIso = false;
                if (dumpIsoInfo) {
                  sprintf(messageBuffer, "Photon number %d failed isolation with Jet number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iJet, runno, lumiblock, evtno);
                  outTextFile << messageBuffer;
                }
                break;
              }
            }
          }
          if (!passIso) continue;

          // Throw out photon if there's a photon/muon overlap
          for (int iMuon = 0; iMuon < 25; ++iMuon ) {
            if (MuEt[iMuon]>50 && MuPFdBiso[iMuon]<0.15 && dR(PhEta[iPhoton], PhPhi[iPhoton], MuEta[iMuon], MuPhi[iMuon]) < 0.3) {
              if (dumpIsoInfo) {
                sprintf(messageBuffer, "Photon number %d failed isolation with Muon number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iMuon, runno, lumiblock, evtno);
                outTextFile << messageBuffer;
              }
              passIso = false;
              break;
            }
          }
          if (!passIso) continue;

          // Throw out photon if there's a photon/electron overlap
          for (int iElectron = 0; iElectron < 25; ++iElectron ) {
            if (EleEt[iElectron]>50 && dR(PhEta[iPhoton], PhPhi[iPhoton], EleEta[iElectron], ElePhi[iElectron]) < 0.3) {
              if (dumpIsoInfo) {
                sprintf(messageBuffer, "Photon number %d failed isolation with Electron number %d  in run number %d lumi section %d event number %lld\n", iPhoton, iElectron, runno, lumiblock, evtno);
                outTextFile << messageBuffer;
              }
              passIso = false;
              break;
            }
          }
          if (!passIso) continue;

          if (debugFlag) cout << "    PhEt for photon number " << iPhoton << " is: " << PhEt[iPhoton] << endl;
          ST += PhEt[iPhoton];
          multiplicity+=1;
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "Photon number %d passed isolation in run number %d lumi section %d event number %lld.\n      It had Px=%f and Py=%f\n", iPhoton, runno, lumiblock, evtno, PhPx[iPhoton], PhPy[iPhoton]);
            outTextFile << messageBuffer;
          }
          Px += PhPx[iPhoton];
          Py += PhPy[iPhoton];
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
            outTextFile << messageBuffer;
          }
        }
        else break;
      }
    }

    //Muons
    if (eventHasMuon) {
      for (int iMuon = 0; iMuon < 25; ++iMuon) {
        passIso=true;
        if (MuEt[iMuon]>50. && MuPFdBiso[iMuon]<0.15) {
          if (debugFlag) cout << "    MuEt for muon number " << iMuon << " is: " << MuEt[iMuon] << endl;
          ST += MuEt[iMuon];
          multiplicity+=1;
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "Muon number %d passed isolation in run number %d lumi section %d event number %lld.\n       It had Px=%f and Py=%f\n", iMuon, runno, lumiblock, evtno, MuPx[iMuon], MuPy[iMuon]);
            outTextFile << messageBuffer;
          }
          Px += MuPx[iMuon];
          Py += MuPy[iMuon];
          if (debugFlag && dumpIsoInfo) {
            sprintf(messageBuffer, "   Cumulative: Px=%f and Py=%f\n", Px, Py);
            outTextFile << messageBuffer;
          }
        }
        else break;
      }
    }


    //debug info and big ST printing
    if (debugFlag) cout << "    Met from PAT collection is: " << Met << endl;
    STMHTnoMET = ST + OurMet;
    ST += Met;
    OurMet = std::sqrt(Px*Px + Py*Py);
    if (debugFlag) cout << "    Met calculated according to my recipe is: " << OurMet << endl;
    stHist.Fill(ST);
    stHistMHT.Fill(STMHTnoMET);
    for (int iHist = 0; iHist<multMax-2; ++iHist) {
      if (multiplicity == iHist+2) stExcHist[iHist]->Fill(ST);
      if (multiplicity >= iHist+2) stIncHist[iHist]->Fill(ST);
    }
    for (int iHist = 0; iHist<multMax-2; ++iHist) {
      if (multiplicity == iHist+2) stExcHistMHT[iHist]->Fill(STMHTnoMET);
      if (multiplicity >= iHist+2) stIncHistMHT[iHist]->Fill(STMHTnoMET);
    }
    METvsMHT.Fill(OurMet,Met);
    if (multiplicity>=2){
      METvsMHTinc2.Fill(OurMet,Met);
      if (eventHasMuon)                                           METvsMHTinc2hasMuon.Fill(OurMet, Met);
      if (eventHasPhoton)                                         METvsMHTinc2hasPhoton.Fill(OurMet, Met);
      if (eventHasElectron)                                       METvsMHTinc2hasElectron.Fill(OurMet, Met);
      if (!eventHasMuon && !eventHasPhoton && !eventHasElectron)  METvsMHTinc2onlyJets.Fill(OurMet, Met);
    }
    if (dumpIsoInfo && fabs(OurMet-Met)>300) {
      sprintf(messageBuffer, "MET-MHT is %f in run number %d lumi section %d event number %lld. ST is %f and multiplicity is %d\n", Met-OurMet, runno, lumiblock, evtno, ST, multiplicity);
      outTextFile << messageBuffer;
      if (debugFlag) cout << messageBuffer;
    }


    // dump info on events with
    if ((ST>3000 || STMHTnoMET>4000 || Met > 500 || OurMet > 500 || fabs(OurMet-Met)>300) && dumpBigEvents) {
      if (debugFlag) cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " ST is:" << ST << endl;
      if (debugFlag) cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " ST with MHT is:" << STMHTnoMET << endl;
    }
    if (( STMHTnoMET>4000 || Met > 500 || OurMet > 500 || fabs(OurMet-Met)>300) && dumpBigEvents) {
      if (debugFlag) cout << "In run number " << runno << " lumi section " << lumiblock << " event number " << evtno << " sT is:" << ST << endl;
      if (dumpIsoInfo) {
        sprintf(messageBuffer, "In run number %d lumi section %d event number %lld ST is %f and multiplicity is %d\n", runno, lumiblock, evtno, ST, multiplicity);
        outTextFile << messageBuffer;
      }
      if (debugFlag) cout << messageBuffer;

      // dump all object info
      for (int j=0; j<25; ++j) {
        if(debugFlag && dumpIsoInfo && JetEt[j]>0.000) {
          sprintf(messageBuffer, "    Jet %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, JetEt[j], JetPx[j], JetPy[j], JetEta[j], JetPhi[j]);
          outTextFile << messageBuffer;
        }
        if (debugFlag) cout  << messageBuffer;
      }
      for (int j=0; j<25; ++j) {
        if(debugFlag && dumpIsoInfo && EleEt[j]>0.000) {
          sprintf(messageBuffer, "    Ele %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, EleEt[j], ElePx[j], ElePy[j], EleEta[j], ElePhi[j]);
          outTextFile << messageBuffer;
        }
        if (debugFlag) cout  << messageBuffer;
      }
      for (int j=0; j<25; ++j) {
        if(debugFlag && dumpIsoInfo && PhEt[j]>0.000) {
          sprintf(messageBuffer, "    Ph %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, PhEt[j], PhPx[j], PhPy[j], PhEta[j], PhPhi[j]);
          outTextFile << messageBuffer;
        }
        if (debugFlag) cout  << messageBuffer;
      }
      for (int j=0; j<25; ++j) {
        if(debugFlag && dumpIsoInfo && MuEt[j]>0.000) {
          sprintf(messageBuffer, "    Mu %d has Et=%f, Px=%f, Py=%f, Eta=%f, Phi=%f\n", j, MuEt[j], MuPx[j], MuPy[j], MuEta[j], MuPhi[j]);
          outTextFile << messageBuffer;
        }
        if (debugFlag) cout  << messageBuffer;
      }
      if (debugFlag && dumpIsoInfo) {
        sprintf(messageBuffer, "    our Px is=%f\n", Px);
        outTextFile << messageBuffer;
        sprintf(messageBuffer, "    our Py is=%f\n", Py);
        outTextFile << messageBuffer;
        sprintf(messageBuffer, "    our MHT is=%f\n", OurMet);
        outTextFile << messageBuffer;
        sprintf(messageBuffer, "    MET is=%f\n", Met);
        outTextFile << messageBuffer;
        if (debugFlag) cout  << messageBuffer;
        sprintf(messageBuffer, "\n\n\n\n");
        outTextFile << messageBuffer;
      }
    }
    nDumpedEvents+=1;
    if (debugFlag && nDumpedEvents==eventsToDump) break;
  }
  // write output textfile
  outTextFile.close();
  // write output root file
  TFile* outRootFile = new TFile(outFilename.c_str(), "RECREATE");
  outRootFile->cd();
  outRootFile->mkdir("ST");
  outRootFile->mkdir("MET-MHT");
  outRootFile->mkdir("Isolation");

  outRootFile->cd("ST");
  stHist.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHist[iHist]->Write();
    stIncHist[iHist]->Write();
  }
  stHistMHT.Write();
  for (int iHist = 0; iHist<multMax-2; ++iHist) {
    stExcHistMHT[iHist]->Write();
    stIncHistMHT[iHist]->Write();
  }

  outRootFile->cd("MET-MHT");
  METvsMHT.Write();
  METvsMHTinc2.Write();
  METvsMHTinc2hasMuon.Write();
  METvsMHTinc2hasPhoton.Write();
  METvsMHTinc2hasElectron.Write();
  METvsMHTinc2onlyJets.Write();

  outRootFile->cd("Isolation");
  MuonJetIso1.Write();
  MuonJetIso2.Write();
  MuonJetIso3.Write();
  MuonJetIso4.Write();
  MuonJetoverlapdR1.Write();
  MuonJetoverlapdR2.Write();
  MuonJetoverlapdR3.Write();
  MuonJetoverlapdR4.Write();
  ElectronJetIso1.Write();
  ElectronJetIso2.Write();
  ElectronJetIso3.Write();
  ElectronJetIso4.Write();
  ElectronJetoverlapdR1.Write();
  ElectronJetoverlapdR2.Write();
  ElectronJetoverlapdR3.Write();
  ElectronJetoverlapdR4.Write();
  PhotonJetIso1.Write();
  PhotonJetIso2.Write();
  PhotonJetIso3.Write();
  PhotonJetIso4.Write();
  PhotonJetoverlapdR1.Write();
  PhotonJetoverlapdR2.Write();
  PhotonJetoverlapdR3.Write();
  PhotonJetoverlapdR4.Write();

  outRootFile->Close();
}

//Event list object for MET filtering


// function to calculate dR between two objects
float dR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + std::pow(TMath::ATan2(TMath::Sin( phi1 - phi2), TMath::Cos(phi1-phi2)),2) );
}


std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName) {
  std::map<unsigned, std::set<unsigned> > list;
  ifstream listFile(_fileName);
  if (!listFile.is_open())
    throw std::runtime_error(_fileName);

  unsigned iL(0);
  std::string line;
  while (true) {
    std::getline(listFile, line);
    if (!listFile.good())
      break;

    if (line.find(":") == std::string::npos || line.find(":") == line.rfind(":"))
      continue;

    unsigned run(std::atoi(line.substr(0, line.find(":")).c_str()));
    unsigned event(std::atoi(line.substr(line.rfind(":") + 1).c_str()));

    list[run].insert(event);

    ++iL;
  }

  std::cout << "Loaded " << iL << " events" << std::endl;

  return list;
}

