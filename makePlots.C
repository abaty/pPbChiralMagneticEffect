#include "TMath.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "Settings.h"
#include <iostream>
#include <string>
#include <vector>

void makePlots(){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  Settings s;

  TH1D * samePlus = new TH1D("samePlus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * sameMinus = new TH1D("sameMinus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoPlus = new TH1D("oppoPlus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoMinus = new TH1D("oppoMinus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);

  TFile * output = TFile::Open("output.root","read");
  TNtuple * QSkim[s.trkEtaGaps];
  for(int g = 0; g<s.trkEtaGaps; g++) QSkim[g] = (TNtuple*) output->Get(Form("QSkim_%d",g));

  for(int g = 0; g<s.trkEtaGaps; g++){
    float QaQbQ2c_ppPlus, avg_QaQbQ2c_ppPlus=0; 
    float QaQbQ2c_ppMinus, avg_QaQbQ2c_ppMinus=0; 
    float QaQbQ2c_pmPlus, avg_QaQbQ2c_pmPlus=0; 
    float QaQbQ2c_pmMinus, avg_QaQbQ2c_pmMinus=0;
    float QhfpQhfm, avg_QhfpQhfm=0;
    float QhfmQhfp, avg_QhfmQhfp=0;
    float QhfpQtrk, avg_QhfpQtrk=0;
    float QhfmQtrk, avg_QhfmQtrk=0;
    QSkim[g]->SetBranchAddress("QaQbQ2c_ppPlus",&QaQbQ2c_ppPlus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_ppMinus",&QaQbQ2c_ppMinus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_pmPlus",&QaQbQ2c_pmPlus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_pmMinus",&QaQbQ2c_pmMinus);
    
    QSkim[g]->SetBranchAddress("QhfpQhfm",&QhfpQhfm);
    QSkim[g]->SetBranchAddress("QhfmQhfp",&QhfmQhfp);
    QSkim[g]->SetBranchAddress("QhfpQtrk",&QhfpQtrk);
    QSkim[g]->SetBranchAddress("QhfmQtrk",&QhfmQtrk);

    float n = 0;
    for(int i = 0; i<QSkim[g]->GetEntries(); i++){
      QSkim[g]->GetEntry(i);
      avg_QaQbQ2c_ppPlus +=QaQbQ2c_ppPlus; 
      avg_QaQbQ2c_ppMinus+=QaQbQ2c_ppMinus; 
      avg_QaQbQ2c_pmPlus +=QaQbQ2c_pmPlus; 
      avg_QaQbQ2c_pmMinus+=QaQbQ2c_pmMinus;

      avg_QhfpQhfm+=QhfpQhfm;
      avg_QhfmQhfp+=QhfmQhfp;
      avg_QhfpQtrk+=QhfpQtrk;
      avg_QhfmQtrk+=QhfmQtrk; 
      n++;
    }
    avg_QaQbQ2c_ppPlus = avg_QaQbQ2c_ppPlus/n;
    avg_QaQbQ2c_ppMinus = avg_QaQbQ2c_ppMinus/n;
    avg_QaQbQ2c_pmPlus = avg_QaQbQ2c_pmPlus/n;
    avg_QaQbQ2c_pmMinus = avg_QaQbQ2c_pmMinus/n;

    avg_QhfpQhfm = avg_QhfpQhfm/n;
    avg_QhfmQhfp = avg_QhfmQhfp/n;
    avg_QhfpQtrk = avg_QhfpQtrk/n;
    avg_QhfmQtrk = avg_QhfmQtrk/n;

    float v2p = TMath::Power(avg_QhfpQhfm*avg_QhfpQtrk/avg_QhfmQtrk,0.5);
    float v2m = TMath::Power(avg_QhfmQhfp*avg_QhfmQtrk/avg_QhfpQtrk,0.5);
   
    std::cout << v2p << " " << v2m << std::endl;
 
    samePlus->SetBinContent(g+1,avg_QaQbQ2c_ppPlus/v2p);
    sameMinus->SetBinContent(g+1,avg_QaQbQ2c_ppPlus/v2m);
    oppoPlus->SetBinContent(g+1,avg_QaQbQ2c_pmPlus/v2p);
    oppoMinus->SetBinContent(g+1,avg_QaQbQ2c_pmPlus/v2m);
  }
  samePlus->Print("All");
  sameMinus->Print("All");
  oppoPlus->Print("All");
  oppoMinus->Print("All");
  TCanvas * c1 = new TCanvas("C1","C1");
  samePlus->Draw("p");
  sameMinus->SetMarkerStyle(24);
  sameMinus->Draw("same p");
  oppoPlus->SetMarkerColor(kRed);
  oppoPlus->Draw("same p");
  oppoMinus->SetMarkerColor(kRed);
  oppoMinus->SetMarkerStyle(24);
  oppoMinus->Draw("same p");
  c1->SaveAs("CME.png");
  c1->SaveAs("CME.pdf");
}
