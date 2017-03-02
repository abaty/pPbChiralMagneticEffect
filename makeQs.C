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
#include "TLorentzVector.h"

float roundEta(float eta){
  return (int)(10*eta)/10.0+((eta>=0)?0.05:-0.05);
}

void makeQs(std::vector<std::string> inputFiles, int job){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  
  Settings s;

  TNtuple * QSkim[s.trkEtaGaps];
  std::string QSkimVars;
  QSkimVars=   "QaQbQ2c_ppPlus:QaQbQ2c_ppMinus:QaQbQ2c_pmPlus:QaQbQ2c_pmMinus:QhfpQhfm:QhfmQhfp:QhfpQtrk:QhfmQtrk:wQaQbQ2c_ppPlus:wQaQbQ2c_ppMinus:wQaQbQ2c_pmPlus:wQaQbQ2c_pmMinus:wQhfpQhfm:wQhfmQhfp:wQhfpQtrk:wQhfmQtrk:posFrac:negFrac:nTrkOffline";

  TFile * f = TFile::Open("EPOS_eff.root","Read");
  TH2D * trkCorr = (TH2D*)f->Get("recoHist");
  trkCorr->SetDirectory(0);
  f->Close();

  int nTrk;
  int nVtx;
  int nRun;
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
  int trkNPixlayer[50000];
  unsigned char trkNdof[50000];
  float trkAlgo[50000];
  unsigned char trkOriginalAlgo[50000];
  int trkCharge[50000];

  int mult100, mult130, mult160, mult190, mult220;
  int mult100_v2, mult130_v2, mult160_v2, mult190_v2, mult220_v2;
  int mult100_v1, mult130_v1, mult160_v1, mult190_v1, mult220_v1;

  int pVtx;
  int pBeamScrape;

  int HFn;
  float HFEt[20000];
  float HFeta[20000];
  float HFphi[20000];
  int pVtxGplus = 0;
  int phfCoinc = 0;

  TFile * input, * output;
  std::cout << inputFiles.size() << std::endl;
  for(int f = 0; f<inputFiles.size(); f++){ 
  input = TFile::Open(inputFiles.at(f).c_str(),"read");

  TTree * trkCh = (TTree*)input->Get("ppTrack/trackTree");
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("nVtx",&nVtx);
  trkCh->SetBranchAddress("nRun",&nRun);
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
  trkCh->SetBranchAddress("trkNPixlayer",&trkNPixlayer);
  trkCh->SetBranchAddress("trkNdof",&trkNdof);
  trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
  trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
  trkCh->SetBranchAddress("zVtx",&zVtx);
  trkCh->SetBranchAddress("trkCharge",&trkCharge);
  
  TTree * hltCh = (TTree*)input->Get("hltanalysis/HltTree");
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity100_v1",&mult100_v1);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity130_v1",&mult130_v1);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity160_v1",&mult160_v1);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity190_v1",&mult190_v1);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity220_v1",&mult220_v1);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity100_v2",&mult100_v2);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity130_v2",&mult130_v2);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity160_v2",&mult160_v2);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity190_v2",&mult190_v2);
  hltCh->SetBranchAddress("HLT_PAPixelTracks_Multiplicity220_v2",&mult220_v2);

  TTree * evtCh = (TTree*)input->Get("skimanalysis/HltTree");
  evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
  evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
  evtCh->SetBranchAddress("pVertexFilterCutGplus",&pVtxGplus);
  evtCh->SetBranchAddress("phfCoincFilter",&phfCoinc);
 
  TTree * towerCh = (TTree*)input->Get("rechitanalyzer/tower"); 
  towerCh->SetBranchAddress("n",&HFn);
  towerCh->SetBranchAddress("et",&HFEt);
  towerCh->SetBranchAddress("eta",&HFeta);
  towerCh->SetBranchAddress("phi",&HFphi);

  output = TFile::Open(Form("pPbCME_output__%d.root",f),"recreate");
  TH1D * cutEvt = new TH1D("cutEvt","cutEvt",20,0,20);
  TH1D * chPhip = new TH1D("chPhip","chPhip",360,-TMath::Pi(),TMath::Pi());
  TH1D * chPhin = new TH1D("chPhin","chPhin",360,-TMath::Pi(),TMath::Pi());

  TH1D * kshortSame = new TH1D("kshortSame","kshortSame",50,0.4,0.6);
  TH1D * kshortOppo = new TH1D("kshortOppo","kshortOppo",50,0.4,0.6);
  for(int g = 0; g<s.trkEtaGaps; g++){ QSkim[g]  = new TNtuple(Form("QSkim%d",g),"",QSkimVars.data());}

  /*
  std::string QSkimVars;
  QSkimVars=   "QaQbQ2c_ppPlus:QaQbQ2c_ppMinus:QaQbQ2c_pmPlus:QaQbQ2c_pmMinus:QhfpQhfm:QhfmQhfp:QhfpQtrk:QhfmQtrk:wQaQbQ2c_ppPlus:wQaQbQ2c_ppMinus:wQaQbQ2c_pmPlus:wQaQbQ2c_pmMinus:wQhfpQhfm:wQhfmQhfp:wQhfpQtrk:wQhfmQtrk:posFrac:negFrac:nTrkOffline";
  for(int g = 0; g<s.trkEtaGaps; g++){ QSkim[g]  = new TNtuple(Form("QSkim_%d",g),"",QSkimVars.data());}
  */ 

  for(int i = 0; i<evtCh->GetEntries(); i++){
    evtCh->GetEntry(i);
    cutEvt->Fill(1);
    if(i%1000==0) std::cout << i << "/" << evtCh->GetEntries() << std::endl;
    if(!(pVtx && pBeamScrape && pVtxGplus && phfCoinc)) continue;
    cutEvt->Fill(2);

    hltCh->GetEntry(i);
    mult100 = (mult100_v1==1) || (mult100_v2==1);
    mult130 = (mult130_v1==1) || (mult130_v2==1);
    mult160 = (mult160_v1==1) || (mult160_v2==1);
    mult190 = (mult190_v1==1) || (mult190_v2==1);
    mult220 = (mult220_v1==1) || (mult220_v2==1);
    if(!(mult100 || mult130 || mult160)) continue;
    cutEvt->Fill(3);

    trkCh->GetEntry(i);
    if(TMath::Abs(zVtx[0])>15) continue;
    cutEvt->Fill(4);
  
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
    double trkOffline = 0;
    double weight[20000] = {0};
    float totP = 0, totN = 0, tot = 0;
    for(int j = 0; j<nTrk; j++){
      if(TMath::Abs(trkEta[j])>s.trkEtaCut) continue;
      if(trkPt[j]<s.ptMin) continue;
      if(!highPurity[j]) continue;
      if(trkPtError[j]/trkPt[j]>0.1) continue;
      if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3 ) continue;
      if(trkPt[j]>0.4) trkOffline++;

      if(trkNPixlayer[j]<1) continue;
      if(trkPt[j]>s.ptMax) continue;

      if(trkCharge[j]>0){
        totP++;
        chPhip->Fill(trkPhi[j]);
      }
      if(trkCharge[j]<0){
        totN++;
        chPhin->Fill(trkPhi[j]);
      }
      tot++;
      weight[j]=1.;
      if(s.doPlusTrkCorr==true && trkCharge[j]>0) weight[j]=1./trkCorr->GetBinContent(trkCorr->GetXaxis()->FindBin(trkEta[j]),trkCorr->GetYaxis()->FindBin(trkPt[j]));
      if(s.doMinusTrkCorr==true && trkCharge[j]<0) weight[j]=1./trkCorr->GetBinContent(trkCorr->GetXaxis()->FindBin(trkEta[j]),trkCorr->GetYaxis()->FindBin(trkPt[j]));
    }

    //std::cout << nTrk << " " << trkOffline << std::endl; 
    if(trkOffline<s.nTrkMin || trkOffline>=s.nTrkMax) continue;
    cutEvt->Fill(5);
    
    for(int j = 0; j<nTrk; j++){
      if(weight[j]==0) continue;
      
      //Q2 for v2
      TComplex q2 = TComplex(weight[j],2*trkPhi[j],true);
      Q2trk_Both += q2; 
      totalW_Both += weight[j];
        
      TLorentzVector *pi1 = new TLorentzVector(0,0,0,0);
      TLorentzVector *pi2 = new TLorentzVector(0,0,0,0);


      for(int jj = j+1; jj<nTrk; jj++){
        if(weight[jj]==0) continue;

        if(trkCharge[j]==trkCharge[jj]){
          pi1->SetPtEtaPhiM(trkPt[j],trkEta[j],trkPhi[j],0.13957);
          pi2->SetPtEtaPhiM(trkPt[jj],trkEta[jj],trkPhi[jj],0.13957);
          *pi1 = *pi1 + *pi2;
          kshortSame->Fill(pi1->M());
        }else{
          pi1->SetPtEtaPhiM(trkPt[j],trkEta[j],trkPhi[j],0.13957);
          pi2->SetPtEtaPhiM(trkPt[jj],trkEta[jj],trkPhi[jj],0.13957);
          *pi1 = *pi1 + *pi2;
          kshortOppo->Fill(pi1->M()); 
        }

        for(int g = 0; g<s.trkEtaGaps; g++){
          if(!s.doDiscreteEta && TMath::Abs(trkEta[j]-trkEta[jj])<s.etaGaps[g]) continue;
          if(!s.doDiscreteEta && TMath::Abs(trkEta[j]-trkEta[jj])>=s.etaGaps[g+1]) continue;
          if(s.doDiscreteEta && TMath::Abs(roundEta(trkEta[j])-roundEta(trkEta[jj]))<s.etaGaps[g]) continue;
          if(s.doDiscreteEta && TMath::Abs(roundEta(trkEta[j])-roundEta(trkEta[jj]))>=s.etaGaps[g+1]) continue;
 
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
      delete pi1;
      delete pi2;
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
    int rev = 1;
    if(nRun<=211312) rev = -1;//swapp eta for early runs (because beam was reversed)

    for(int j = 0; j<HFn; j++){
      if(TMath::Abs(HFeta[j])>s.HFetaMax || TMath::Abs(HFeta[j])<s.HFetaMin) continue;
      if(rev*HFeta[j]>0){
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
    cutEvt->Fill(6);
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

      if(TMath::IsNaN(QaQbQ2c_ppPlus.Re()) || TMath::IsNaN(QaQbQ2c_ppMinus.Re()) || TMath::IsNaN(QaQbQ2c_pmPlus.Re()) || TMath::IsNaN(QaQbQ2c_pmMinus.Re()) || TMath::IsNaN(QhfpQhfm.Re()) || TMath::IsNaN(QhfmQhfp.Re()) || TMath::IsNaN(QhfpQtrk.Re()) || TMath::IsNaN(QhfmQtrk.Re())){std::cout << "ERROR!!! GOT NAN" << std::endl; continue;}  
       

      float skimEntry[] = {(float)QaQbQ2c_ppPlus.Re(),(float)QaQbQ2c_ppMinus.Re(),(float)QaQbQ2c_pmPlus.Re(),(float)QaQbQ2c_pmMinus.Re(),(float)QhfpQhfm.Re(),(float)QhfmQhfp.Re(),(float)QhfpQtrk.Re(),(float)QhfmQtrk.Re(),(float)(totalW_pp[g]*totalEt_Plus),(float)(totalW_pp[g]*totalEt_Minus),(float)(totalW_pm[g]*totalEt_Plus),(float)(totalW_pm[g]*totalEt_Minus),(float)(totalEt_Plus*totalEt_Minus),(float)(totalEt_Minus*totalEt_Plus),(float)(totalEt_Plus*totalW_Both),(float)(totalEt_Minus*totalW_Both),totP/tot,totN/tot,(float)trkOffline};
      QSkim[g]->Fill(skimEntry);
     }
    //std::cout << QSkim[0]->GetEntries() << std::endl;
  }//close evt loop
    for(int g = 0; g<s.trkEtaGaps; g++){  QSkim[g]->Write();}
    cutEvt->Write();
    chPhip->Write();
    chPhin->Write();
    kshortSame->Write();
    kshortOppo->Write();
    output->Close();
//    delete cutEvt;
//    for(int g = 0; g<s.trkEtaGaps; g++){delete  QSkim[g];}
    input->Close();  
  }//close file loop
}

//***************************************************************************************************************

int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: countTracks <fileList>  <job>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  makeQs(listOfFiles,job);
  return 0; 
}

