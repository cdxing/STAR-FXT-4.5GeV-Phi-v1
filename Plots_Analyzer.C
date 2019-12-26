#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"

void Plots_Analyzer(const char* inFile = "./result/97D2AFB119464E2BB781CF04F1524B05.picoDst.result.combined.root")
{
  TFile* F = TFile::Open("./result/97D2AFB119464E2BB781CF04F1524B05.picoDst.result.combined.root");
  if (!F || F->IsZombie()){
    std::cout<<"Cannot open the file"<<inFile<<std::endl;
    return;
  }

  TH2D* H2 = (TH2D*)F->Get("h2_m2_QA_pq");

  TH1D* H1 = H2->ProjectionY("h1_m2",0,-1);

  double d_PRO_m2  = 0.880354;
  double d_K_m2    = 0.24371698;
  double d_PI_m2   = 0.019479835;

  TF1* F1 = new TF1("FitF1","[0] * TMath::Gaus(x, 0.880354, [1])", d_PRO_m2-0.13,  d_PRO_m2+0.15);
  F1->SetParameter(0,37095);
  F1->SetParameter(1,0.8896);
  F1->SetParameter(2,0.04846);

  TF1* F2 = new TF1("FitF2","[0] * TMath::Gaus(x, 0.24371698, [1])", d_K_m2-0.027,  d_K_m2+0.03);
  F2->SetParameter(0,8000);
  F2->SetParameter(1,0.24);
  F2->SetParameter(2,0.0176);

  TF1* F3 = new TF1("FitF3","[0] * TMath::Gaus(x, 0.019479835, [1])", 0.019-0.004,  0.019+0.007);
  F3->SetParameter(0,600000);
  F3->SetParameter(1,0.019);
  F3->SetParameter(2,0.004);

  // TF1* F4 = new TF1("FitF4","[0] +[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x^7+[8]*x^8+[9]*x^9",0,1.2);


  H1->Fit(F1->GetName(), "ER0+", "kRed", d_PRO_m2-0.1, d_PRO_m2+0.12);
  std::cout<< "Chisquare = "<<F1->GetChisquare()<<std::endl;
  std::cout<< "NDF = "<<F1->GetNDF()<<std::endl;

  H1->Fit(F2->GetName(), "ER0+", "kGreen", d_K_m2-0.015, d_K_m2+0.015);
  std::cout<< "Chisquare = "<<F2->GetChisquare()<<std::endl;
  std::cout<< "NDF = "<<F2->GetNDF()<<std::endl;

  H1->Fit(F3->GetName(), "ER0+", "kBlue", 0.019-0.001, 0.019+0.003);
  std::cout<< "Chisquare = "<<F3->GetChisquare()<<std::endl;
  std::cout<< "NDF = "<<F3->GetNDF()<<std::endl;

  // H1->Fit(F4->GetName(), "ER0+", "kBlue", -0.2, 1.6);



  TCanvas* c1 = new TCanvas("c1", "Analyzer", 1600, 800);
  // c1->Divide(2,1);
  // c1->cd(1);
  // H2->DrawCopy("colz");
  // c1->cd(2);
  c1->cd();
  H1->DrawCopy();
  F1->SetLineColor(2);
  F2->SetLineColor(2);
  F3->SetLineColor(2);
  // F4->SetLineColor(4);
  F1->SetLineStyle(2);
  F2->SetLineStyle(2);
  F3->SetLineStyle(2);
  // F4->SetLineStyle(2);
  F1->Draw("same");
  F2->Draw("same");
  F3->Draw("same");
  // F4->Draw("same");
  // TLegend* leg1 = new TLegend(0.1,0.7,0.48,0.9);
  // leg1->AddEntry(F2,"K_m2","lp");
  // leg1->Draw("same");

  return;
}
