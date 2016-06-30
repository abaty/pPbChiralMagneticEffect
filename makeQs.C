#include "TMath.h"
#include "TComplex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
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
 
  TFile * input = new TFile("/mnt/hadoop/cms/store/user/abaty/unmergedForests/pPb_HighMultiplicityForests_AfterBeamReversal_NoJets/store/user/abaty/unmergedForests/pPb_HighMultiplicityForests_AfterBeamReversal_NoJets/PAHighPt/crab_20160629_200757/160629_180823/0000/HiForest_11.root","read");
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

  TFile * output = TFile::Open("output.root","recreate");
  std::string QSkimVars;
  QSkimVars=   "QaQbQ2c_ppPlus:QaQbQ2c_ppMinus:QaQbQ2c_pmPlus:QaQbQ2c_pmMinus:QhfpQhfm:QhfmQhfp:QhfpQtrk:QhfmQtrk";
  TNtuple * QSkim[s.trkEtaGaps];
  for(int g = 0; g<s.trkEtaGaps; g++) QSkim[g]  = new TNtuple(Form("QSkim_%d",g),"",QSkimVars.data()); 
 
  for(int i = 0; i<evtCh->GetEntries(); i++){
    evtCh->GetEntry(i);
    if(i%1000==0) std::cout << i << "/" << evtCh->GetEntries() << std::endl;
    if(!(pVtx && pBeamScrape)) continue;

    hltCh->GetEntry(i);
    if(!(mult100 || mult130 || mult160 || mult190)) continue;

    trkCh->GetEntry(i);
    if(nTrk<s.nTrkMin || nTrk>s.nTrkMax) continue;

    //trk loop for Qs
    TComplex Q2trk_Both = TComplex(0,0);
    double totalW_Both = 0;
    TComplex Q1trk_pp[s.trkEtaGaps]; 
    TComplex Q1trk_pm[s.trkEtaGaps]; 
    TComplex Q1trk_mm[s.trkEtaGaps]; 
    double totalW_pp[s.trkEtaGaps] = {0};
    double totalW_pm[s.trkEtaGaps] = {0};
    double totalW_mm[s.trkEtaGaps] = {0};
    for(int g = 0; g<s.trkEtaGaps; g++){
      Q1trk_pp[g] = TComplex(0,0);
      Q1trk_pm[g] = TComplex(0,0);
      Q1trk_mm[g] = TComplex(0,0); 
    }
    //selection skim
    double weight[20000] = {0};
    for(int j = 0; j<nTrk; j++){
      if(TMath::Abs(zVtx[0]>15)) continue;
      if(TMath::Abs(trkEta[j])>s.trkEtaCut) continue;
      if(trkPt[j]<s.ptMin || trkPt[j]>s.ptMax) continue;
      if(!highPurity[j]) continue;
      if(trkPtError[j]/trkPt[j]>0.1) continue;
      weight[j]=1;
    }
    for(int j = 0; j<nTrk; j++){
      if(weight[j]==0) continue;
      
      //Q2 for v2
      TComplex q2 = TComplex(weight[j],2*trkPhi[j],true);
      Q2trk_Both += q2;
      totalW_Both += weight[j];

      for(int jj = j+1; jj<nTrk; jj++){
        if(weight[jj]==0) continue;
        for(int g = 0; g<s.trkEtaGaps; g++){
          if(TMath::Abs(trkEta[j]-trkEta[jj])<s.etaGaps[g]) continue;
          //if(TMath::Abs(trkEta[j]-trkEta[jj])>=s.etaGaps[g+1]) continue;
 
          //Q1
          //if(trkCharge[j]>0 && trkCharge[jj]>0){
          if(trkCharge[j]==trkCharge[jj]){
            TComplex q1 = TComplex(weight[j]*weight[jj],trkPhi[j]+trkPhi[jj],true);
            Q1trk_pp[g] += q1;
            totalW_pp[g] += weight[j]*weight[jj];
          }else if(trkCharge[j]== -trkCharge[jj]){
            TComplex q1 = TComplex(weight[j]*weight[jj],trkPhi[j]+trkPhi[jj],true);
            Q1trk_pm[g] += q1;
            totalW_pm[g] += weight[j]*weight[jj];
          /*}else{
            TComplex q1 = TComplex(weight[j]*weight[jj],trkPhi[j]+trkPhi[jj],true);
            Q1trk_mm += q1;
            totalW_mm += weight[j]*weight[jj];*/
          }
        }
      }
    }   
    Q2trk_Both  = Q2trk_Both/totalW_Both;
    for(int g = 0; g<s.trkEtaGaps; g++){
      Q1trk_pp[g]  = Q1trk_pp[g]/totalW_pp[g];
      Q1trk_pm[g]  = Q1trk_pm[g]/totalW_pm[g];  
      //Q1trk_mm  = Q1trk_mm/totalW_mm;  
    }
 
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
    
    //std::cout << Q2trk_Both.Re() << " " << Q1trk_pp[0].Re()  << " " << Q1trk_pm[0].Re() << " " << Q2hf_Plus.Re() << " " << Q2hf_Minus.Re() << " " << std::endl;

    // for numerator
    for(int g = 0; g<s.trkEtaGaps; g++){
      TComplex QaQbQ2c_ppPlus  = Q1trk_pp[g]*TComplex::Conjugate(Q2hf_Plus);
      TComplex QaQbQ2c_ppMinus = Q1trk_pp[g]*TComplex::Conjugate(Q2hf_Minus);
      TComplex QaQbQ2c_pmPlus  = Q1trk_pm[g]*TComplex::Conjugate(Q2hf_Plus);
      TComplex QaQbQ2c_pmMinus = Q1trk_pm[g]*TComplex::Conjugate(Q2hf_Minus);

      //for v2
      TComplex QhfpQhfm = Q2hf_Plus*TComplex::Conjugate(Q2hf_Minus);
      TComplex QhfmQhfp = Q2hf_Minus*TComplex::Conjugate(Q2hf_Plus);
      TComplex QhfpQtrk = Q2hf_Plus*TComplex::Conjugate(Q2trk_Both);
      TComplex QhfmQtrk = Q2hf_Minus*TComplex::Conjugate(Q2trk_Both);

      float skimEntry[] = {(float)QaQbQ2c_ppPlus.Re(),(float)QaQbQ2c_ppMinus.Re(),(float)QaQbQ2c_pmPlus.Re(),(float)QaQbQ2c_pmMinus.Re(),(float)QhfpQhfm.Re(),(float)QhfmQhfp.Re(),(float)QhfpQtrk.Re(),(float)QhfmQtrk.Re()};
      QSkim[g]->Fill(skimEntry);
     }
  }
  for(int g = 0; g<s.trkEtaGaps; g++) QSkim[g]->Write();
  output->Close();
  input->Close();  
}
