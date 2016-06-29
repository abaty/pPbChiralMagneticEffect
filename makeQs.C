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

    //trk loop for Qs
    TComplex Q2trk_Both = TComplex(0,0);
    TComplex Q1trk_Plus = TComplex(0,0);
    TComplex Q1trk_Minus = TComplex(0,0);
    double totalW_Both = 0;
    double totalW_Plus = 0;
    double totalW_Minus = 0;
    for(int j = 0; j<nTrk; j++){
      if(TMath::Abs(zVtx[0]>15)) continue;
      if(TMath::Abs(trkEta[j])>s.trkEtaCut) continue;
      if(trkPt[j]<s.ptMin || trkPt[j]>s.ptMax) continue;
      if(!highPurity[j]) continue;
      if(trkPtError[j]/trkPt[j]>0.1) continue;
      
      double weight = 1;
      //Q2 for v2
      TComplex q2 = TComplex(weight,2*trkPhi[j],true);
      Q2trk_Both += q2;
      totalW_Both += weight;
      //Q1
      if(trkCharge[j]>0){
        TComplex q1 = TComplex(weight,1*trkPhi[j],true);
        Q1trk_Plus += q1;
        totalW_Plus += weight;
      }else{
        TComplex q1 = TComplex(weight,1*trkPhi[j],true);
        Q1trk_Minus += q1;
        totalW_Minus += weight;
      }
    }   
    Q2trk_Both  = Q2trk_Both/totalW_Both;
    Q1trk_Plus  = Q1trk_Plus/totalW_Plus;
    Q1trk_Minus = Q1trk_Minus/totalW_Minus;  
 
    //HF Tower loop for Qs
    towerCh->GetEntry(i);
    TComplex Q2hf_Plus = TComplex(0,0);
    TComplex Q2hf_Minus = TComplex(0,0);
    double totalEt_Plus = 0;
    double totalEt_Minus = 0;
    for(int j = 0; j<HFn; j++){
      if(TMath::Abs(HFeta[j])>s.HFetaMax || TMath::Abs(HFeta[j])<s.HFetaMin) continue;
      if(HFeta[j]>0){
        TComplex q2 = TComplex(HFEt[j],2*HFphi[j],true);
        Q2hf_Plus += q2;
        totalEt_Plus += HFEt[j];
      }else{
        TComplex q2 = TComplex(HFEt[j],2*HFphi[j],true);
        Q2hf_Minus += q2;
        totalEt_Minus += HFEt[j];
      }
    }
    Q2hf_Plus  = Q2hf_Plus/totalEt_Plus;
    Q2hf_Minus = Q2hf_Minus/totalEt_Minus;
  
    std::cout << Q2trk_Both.Re() << " " << Q1trk_Plus.Re()  << " " << Q1trk_Minus.Re() << " " << Q2hf_Plus.Re() << " " << Q2hf_Minus.Re() << " " << std::endl;
  }  
}
