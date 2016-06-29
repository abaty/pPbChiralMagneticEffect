#include "TMath.h"
#include "TComplex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "Settings.h"
#include <iostream>
#include <string>
#include <vector>

void makeQs(){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  Settings s;

  int nTrk;
  int nVtx;
  int nTrkTimesnVtx;
  bool highPurity[50000];
  float trkPt[50000];
  float trkPtError[50000];
  float trkEta[50000];
  float trkPhi[50000];
  float trkDxy1[50000];
  float trkDxyError1[50000];
  float trkDz1[50000];
  float trkDzError1[50000];
  float trkDzOverDzError[500000];
  float trkDxyOverDxyError[500000];
  float trkChi2[50000];
  float zVtx[20];
  unsigned char trkNHit[50000];
  unsigned char trkNlayer[50000];
  unsigned char trkNdof[50000];
  unsigned char trkAlgo[50000];
  unsigned char trkOriginalAlgo[50000];
  int trkCharge[50000];

  int mult100, mult130, mult160, mult190, mult220;

  int pVtx;
  int pBeamScrape;

  int HFn;
  float HFEt[20000];
  float HFeta[20000];
  float HFphi[20000];
 
  TFile * input = new TFile("HiForest.root","read");
  TTree * trkCh = (TTree*)input->Get("ppTrack/trackTree");
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("nVtx",&nVtx);
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkNHit",&trkNHit);
  trkCh->SetBranchAddress("trkPtError",&trkPtError);
  trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
  trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
  trkCh->SetBranchAddress("trkDz1",&trkDz1);
  trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
  trkCh->SetBranchAddress("trkChi2",&trkChi2);
  trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
  trkCh->SetBranchAddress("trkNdof",&trkNdof);
  trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
  trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
  trkCh->SetBranchAddress("zVtx",&zVtx);
  trkCh->SetBranchAddress("trkCharge",&trkCharge);
  
  TTree * hltCh = (TTree*)input->Get("hltanalysis/HltTree");
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity100_v2",&mult100);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity130_v2",&mult130);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity160_v2",&mult160);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity190_v2",&mult190);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity220_v2",&mult220);

  TTree * evtCh = (TTree*)input->Get("skimanalysis/HltTree");
  evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
  evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
 
  TTree * towerCh = (TTree*)input->Get("rechitanalyzer/tower"); 
  towerCh->SetBranchAddress("n",&HFn);
  towerCh->SetBranchAddress("et",&HFEt);
  towerCh->SetBranchAddress("eta",&HFeta);
  towerCh->SetBranchAddress("phi",&HFphi);

  for(int i = 0; i<evtCh->GetEntries(); i++){
    evtCh->GetEntry(i);
    if(!(pVtx && pBeamScrape)) continue;

    hltCh->GetEntry(i);
    if(!(mult100 || mult130 || mult160 || mult190)) continue;

    trkCh->GetEntry(i);
    if(nTrk<s.nTrkMin || nTrk>s.nTrkMax) continue;
    towerCh->GetEntry(i);
  }
  
}
