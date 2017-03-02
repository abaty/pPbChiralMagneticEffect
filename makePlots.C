#include "TLegend.h"
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

void makePlots(std::string readFile="/export/d00/scratch/abaty/pPbCME_OutputTotal.root"){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  Settings s;

  TProfile * avg[s.trkEtaGaps];

  TH1D * samePlus = new TH1D("samePlus",";#Delta #eta;#LT cos( #phi_{#alpha}+#phi_{#beta}-2 #phi_{c} )#GT/v_{2,c}",s.trkEtaGaps,s.etaGaps);
  TH1D * sameMinus = new TH1D("sameMinus",";#Delta #eta;#LT cos( #phi_{#alpha}+#phi_{#beta}-2 #phi_{c} )#GT/v_{2,c}",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoPlus = new TH1D("oppoPlus",";#Delta #eta;#LT cos( #phi_{#alpha}+#phi_{#beta}-2 #phi_{c} )#GT/v_{2,c}",s.trkEtaGaps,s.etaGaps);
  TH1D * oppoMinus = new TH1D("oppoMinus",";#Delta #eta;#LT cos( #phi_{#alpha}+#phi_{#beta}-2 #phi_{c} )#GT/v_{2,c}",s.trkEtaGaps,s.etaGaps);

  TFile * output = TFile::Open(readFile.c_str(),"read");
  TNtuple * QSkim[s.trkEtaGaps];
  for(int g = 0; g<s.trkEtaGaps; g++) QSkim[g] = (TNtuple*) output->Get(Form("QSkim%d",g));

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

  TFile * out = TFile::Open("PlotsOutput.root","recreate");
  samePlus->Write();
  sameMinus->Write();
  oppoPlus->Write();
  oppoMinus->Write();

  TFile * kong = TFile::Open("Kong_CME_pPb.root","read");
  TH1D * k_sameMinus = (TH1D*) kong->Get("temp1");
  TH1D * k_samePlus =  (TH1D*) kong->Get("temp3");
  TH1D * k_oppoMinus = (TH1D*) kong->Get("temp2");
  TH1D * k_oppoPlus =  (TH1D*) kong->Get("temp4");

  float labelSize = 0.04;

  TCanvas * c1 = new TCanvas("C1","C1");
  sameMinus->SetMarkerColor(kRed);
  sameMinus->SetLineColor(kRed);
  sameMinus->SetLineWidth(1);
  sameMinus->GetYaxis()->SetLabelSize(labelSize);
  sameMinus->GetYaxis()->SetRangeUser(-0.001,0.0015);
  sameMinus->Draw("p");
  samePlus->SetMarkerStyle(24);
  samePlus->SetMarkerColor(kRed);
  samePlus->SetLineColor(kRed);
  samePlus->SetLineWidth(1);
  samePlus->Draw("same p");
  oppoMinus->SetMarkerColor(kBlue);
  oppoMinus->SetLineColor(kBlue);
  oppoMinus->SetLineWidth(1);
  oppoMinus->Draw("same p");
  oppoPlus->SetMarkerStyle(24);
  oppoPlus->SetMarkerColor(kBlue);
  oppoPlus->SetLineColor(kBlue);
  oppoPlus->SetLineWidth(1);
  oppoPlus->Draw("same p");

  TLegend * leg1 = new TLegend(0.2,0.65,0.7,0.92);
  leg1->AddEntry((TObject*)0,"185<N_{trk}^{offline}<220 pPb Crosscheck","");  
  leg1->AddEntry(oppoPlus,"OS p-going","p"); 
  leg1->AddEntry(oppoMinus,"OS Pb-going","p");  
  leg1->AddEntry(samePlus,"SS p-going","p");  
  leg1->AddEntry(sameMinus,"SS Pb-going","p");  
  leg1->Draw("same"); 

  c1->SaveAs("CME.png");
  c1->SaveAs("CME.pdf");

  sameMinus->SetMarkerColor(kBlack);
  sameMinus->SetLineColor(kBlack);
  sameMinus->SetMarkerStyle(24);
  sameMinus->GetYaxis()->SetLabelSize(labelSize);
  sameMinus->Draw("p");
  k_sameMinus->Draw("p same");
  oppoMinus->SetMarkerStyle(25);
  oppoMinus->SetMarkerColor(kBlack);
  oppoMinus->SetLineColor(kBlack);
  oppoMinus->Draw("p same");
  k_oppoMinus->Draw("p same");
  sameMinus->Draw("p same");
  oppoMinus->Draw("p same");
  
  TLegend * leg2 = new TLegend(0.2,0.65,0.7,0.92);
  leg2->AddEntry((TObject*)0,"185<N_{trk}^{offline}<220 Pb-going","");  
  leg2->AddEntry(k_oppoMinus,"OS Paper","p"); 
  leg2->AddEntry(oppoMinus,"OS Crosscheck","p");  
  leg2->AddEntry(k_sameMinus,"SS Paper","p");  
  leg2->AddEntry(sameMinus,"SS Crosscheck","p");  
  leg2->Draw("same"); 
  c1->SaveAs("CME_MinusComparison.png");
  c1->SaveAs("CME_MinusComparison.pdf");
  
  samePlus->SetMarkerColor(kBlack);
  samePlus->SetLineColor(kBlack);
  samePlus->GetYaxis()->SetLabelSize(labelSize);
  samePlus->GetYaxis()->SetRangeUser(-0.001,0.0015);
  samePlus->Draw("p");
  //k_samePlus->SetMarkerStyle(24);
  k_samePlus->Draw("same p");
  oppoPlus->SetMarkerStyle(25);
  oppoPlus->SetMarkerColor(kBlack);
  oppoPlus->SetLineColor(kBlack);
  oppoPlus->Draw("same p");
  k_oppoPlus->SetMarkerStyle(21);
  k_oppoPlus->Draw("same p");
  samePlus->Draw("same p");
  oppoPlus->Draw("same p");
  
  TLegend * leg3 = new TLegend(0.2,0.65,0.7,0.92);
  leg3->AddEntry((TObject*)0,"185<N_{trk}^{offline}<220 p-going","");  
  leg3->AddEntry(k_oppoPlus,"OS Paper","p"); 
  leg3->AddEntry(oppoPlus,"OS Crosscheck","p");  
  leg3->AddEntry(k_samePlus,"SS Paper","p");  
  leg3->AddEntry(samePlus,"SS Crosscheck","p");  
  leg3->Draw("same"); 
 
  c1->SaveAs("CME_PlusComparison.png");
  c1->SaveAs("CME_PlusComparison.pdf");
}
