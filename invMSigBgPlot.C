#include <TStyle.h>

void invMSigBgPlot()
{
  TFile *f1 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_primary.root","READ");
  TFile *f2 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_out_single_tot.root","READ");
  TFile *f3 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_systematic_error_analysis_invM/phi_meson_normal_event_invariant_mass_primary.root","READ");

  TH1D *h_1 = (TH1D*) f2->Get("h_sub");
  TH1D *h_2 = (TH1D*) f1->Get("h_TPC2_TOF3_invM");
  TH1D *h_3 = (TH1D*) f1->Get("h_mx_TPC2_TOF3_invM");
  TH1D *h_4 = (TH1D*) f3->Get("h_mx_TPC2_TOF3_invM");

  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double d_norm  = 0.00296682;//primary//0.0159051;//rebin//0.0159;//5mx_test//0.57346;
  double d_mean  = 1.02087;//primary//1.02094;//rebin//1.02094;//5mx_test//1.02017;//secondfit//1.0202;//prefit//
  double d_sigma = 0.00444335;//primary//0.00452375;//rebin//0.00452732;//5mx_test//0.00463761;//secondfit//0.00483768;//prefit//

  TF1 * tf1 = h_3->GetFunction("tf1_pol");
  TF1 * tf2 = new TF1("tf2","[0] + [1]*x + [2]*x**2",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  double  d_pol_scale_p0 = (tf1->GetParameter(0)) /** 0.5*/ * d_norm;
  double  d_pol_scale_p1 = (tf1->GetParameter(1)) /** 0.5*/ * d_norm;
  double  d_pol_scale_p2 = (tf1->GetParameter(2)) /** 0.5*/ * d_norm;
  // h_3->Fit(tf2, "E","R",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  tf2 -> SetParameter(0,d_pol_scale_p0);
  tf2 -> SetParameter(1,d_pol_scale_p1);
  tf2 -> SetParameter(2,d_pol_scale_p2);
  tf2 -> SetLineColor(kBlue);

  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
  c1->DrawFrame(0.9, 0., 1.1, 500);
  c1->cd();
  h_2->Draw();
  h_2->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
  h_2->SetTitle("");
  h_3->Draw("same");
  h_1->Draw("same");
  tf2->Draw("same");


}
