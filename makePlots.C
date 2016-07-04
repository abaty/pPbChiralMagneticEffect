#include "TMath.h"
#include "TComplex.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TProfile.h"
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

  TProfile * avg[s.trkEtaGaps];

  TH1D * samePlus = new TH1D("samePlus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * sameMinus = new TH1D("sameMinus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoPlus = new TH1D("oppoPlus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoMinus = new TH1D("oppoMinus",";#Delta #eta",s.trkEtaGaps,s.etaGaps);

  TFile * output = TFile::Open("/export/d00/scratch/abaty/pPbCME_OutputTotal.root","read");
  TNtuple * QSkim[s.trkEtaGaps];
  for(int g = 0; g<s.trkEtaGaps; g++) QSkim[g] = (TNtuple*) output->Get(Form("QSkim_%d",g));

  for(int g = 0; g<s.trkEtaGaps; g++){
    avg[g] = new TProfile(Form("avg_%d",g),";",10,0,10);
    float QaQbQ2c_ppPlus, wQaQbQ2c_ppPlus; 
    float QaQbQ2c_ppMinus, wQaQbQ2c_ppMinus; 
    float QaQbQ2c_pmPlus, wQaQbQ2c_pmPlus; 
    float QaQbQ2c_pmMinus, wQaQbQ2c_pmMinus;
    float QhfpQhfm, wQhfpQhfm;
    float QhfmQhfp, wQhfmQhfp;
    float QhfpQtrk, wQhfpQtrk;
    float QhfmQtrk, wQhfmQtrk;
    QSkim[g]->SetBranchAddress("QaQbQ2c_ppPlus",&QaQbQ2c_ppPlus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_ppMinus",&QaQbQ2c_ppMinus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_pmPlus",&QaQbQ2c_pmPlus);
    QSkim[g]->SetBranchAddress("QaQbQ2c_pmMinus",&QaQbQ2c_pmMinus);
    QSkim[g]->SetBranchAddress("wQaQbQ2c_ppPlus",&wQaQbQ2c_ppPlus);
    QSkim[g]->SetBranchAddress("wQaQbQ2c_ppMinus",&wQaQbQ2c_ppMinus);
    QSkim[g]->SetBranchAddress("wQaQbQ2c_pmPlus",&wQaQbQ2c_pmPlus);
    QSkim[g]->SetBranchAddress("wQaQbQ2c_pmMinus",&wQaQbQ2c_pmMinus);
    
    QSkim[g]->SetBranchAddress("QhfpQhfm",&QhfpQhfm);
    QSkim[g]->SetBranchAddress("QhfmQhfp",&QhfmQhfp);
    QSkim[g]->SetBranchAddress("QhfpQtrk",&QhfpQtrk);
    QSkim[g]->SetBranchAddress("QhfmQtrk",&QhfmQtrk);
    QSkim[g]->SetBranchAddress("wQhfpQhfm",&wQhfpQhfm);
    QSkim[g]->SetBranchAddress("wQhfmQhfp",&wQhfmQhfp);
    QSkim[g]->SetBranchAddress("wQhfpQtrk",&wQhfpQtrk);
    QSkim[g]->SetBranchAddress("wQhfmQtrk",&wQhfmQtrk);

    for(int i = 0; i<QSkim[g]->GetEntries(); i++){
    //for(int i = 0; i<50000; i++){
      QSkim[g]->GetEntry(i);
      avg[g]->Fill(1,QaQbQ2c_ppPlus,wQaQbQ2c_ppPlus);
      avg[g]->Fill(2,QaQbQ2c_ppMinus,wQaQbQ2c_ppMinus);
      avg[g]->Fill(3,QaQbQ2c_pmPlus,wQaQbQ2c_pmPlus);
      avg[g]->Fill(4,QaQbQ2c_pmMinus,wQaQbQ2c_pmMinus);

      avg[g]->Fill(5,QhfpQhfm,wQhfpQhfm);
      avg[g]->Fill(6,QhfmQhfp,wQhfmQhfp);
      avg[g]->Fill(7,QhfpQtrk,wQhfpQtrk);
      avg[g]->Fill(8,QhfmQtrk,wQhfmQtrk);
    }
    
//have to add 1 using GenBinCOntent
    std::cout << avg[g]->GetBinContent(5+1) << " " << avg[g]->GetBinContent(6+1)<< " " << avg[g]->GetBinContent(7+1) << " " << avg[g]->GetBinContent(8+1) << std::endl;
    float v2p = TMath::Power(avg[g]->GetBinContent(5+1)*avg[g]->GetBinContent(7+1)/avg[g]->GetBinContent(8+1),0.5);
    float v2m = TMath::Power(avg[g]->GetBinContent(6+1)*avg[g]->GetBinContent(8+1)/avg[g]->GetBinContent(7+1),0.5);
  
    avg[g]->Print("All"); 
    std::cout << v2p << " " << v2m << std::endl;
 
    samePlus->SetBinContent(g+1,avg[g]->GetBinContent(1+1)/v2p);
    sameMinus->SetBinContent(g+1,avg[g]->GetBinContent(2+1)/v2m);
    oppoPlus->SetBinContent(g+1,avg[g]->GetBinContent(3+1)/v2p);
    oppoMinus->SetBinContent(g+1,avg[g]->GetBinContent(4+1)/v2m);
    samePlus->SetBinError(g+1,avg[g]->GetBinError(1+1)/v2p);
    sameMinus->SetBinError(g+1,avg[g]->GetBinError(2+1)/v2m);
    oppoPlus->SetBinError(g+1,avg[g]->GetBinError(3+1)/v2p);
    oppoMinus->SetBinError(g+1,avg[g]->GetBinError(4+1)/v2m);
  }
  samePlus->Print("All");
  sameMinus->Print("All");
  oppoPlus->Print("All");
  oppoMinus->Print("All");
  TCanvas * c1 = new TCanvas("C1","C1");
  sameMinus->SetMarkerColor(kRed);
  sameMinus->Draw("p");
  samePlus->SetMarkerStyle(24);
  samePlus->SetMarkerColor(kRed);
  samePlus->Draw("same p");
  oppoMinus->SetMarkerColor(kBlue);
  oppoMinus->Draw("same p");
  oppoPlus->SetMarkerStyle(24);
  oppoPlus->SetMarkerColor(kBlue);
  oppoPlus->Draw("same p");
  c1->SaveAs("CME.png");
  c1->SaveAs("CME.pdf");
}
