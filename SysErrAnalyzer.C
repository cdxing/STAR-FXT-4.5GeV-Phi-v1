#include <iostream>
#include "Riostream.h"
//#include <unistd.h>
#include <cstdlib>
//#include <pthread.h>

// #include <cstdlib>
// #include <pthread.h>

// //#include <iomanip>
#include <string>
#include <vector>
// //#include <algorithm>
// //#include <cstdlib>
// //#include "gsl/gsl_rng.h"
#include <math.h>
#include <set>
#include <map>

//---
//need this stuff to compile in linux:
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>
#include "TSystemDirectory.h"
//---

#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// #include "TLorentzVector.h"
// #include "TVector3.h"
#include "TMath.h"
// #include "TGraphErrors.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGaxis.h"
#include "TSpectrum.h"
// #include "TLorentzVector.h"
// #include "TVector3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "TProfile.h"

#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TRandom3.h"
#include "TDatime.h"

#include "Math/MinimizerOptions.h"

using namespace std;

void SysErrAnalyzer()
{
  // TFile * tf_input0  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_pol0.root","READ");
  TFile * tf_input1  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_primary.root","READ");
  TFile * tf_input2  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_eff_corr_EP_.root","READ");

  // TFile * tf_input3  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_pol3.root","READ");

  // TFile * tf_input0  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_0_4_.root","READ");
  // TFile * tf_input1  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_1_4_.root","READ");
  // TFile * tf_input2  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_2_4_.root","READ");
  //
  // TFile * tf_input3  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_3_4_.root","READ");
  // TFile * tf_input4  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_4_4_.root","READ");
  // TFile * tf_input5  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_5_4_.root","READ");
  //
  // TFile * tf_input6  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_6_4_.root","READ");
  // TFile * tf_input7  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_7_4_.root","READ");
  // TFile * tf_input8  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_8_4_.root","READ");
  //
  // TFile * tf_input9  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_9_4_.root","READ");
  // TFile * tf_input10  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_sysErr_merged_1_10_4_.root","READ");

  // TH1D * h_input0   = (TH1D*) tf_input0 -> Get("h_phi_v1_vs_y");
  TH1D * h_input1   = (TH1D*) tf_input1 -> Get("h_phi_v1_vs_y");
  TH1D * h_input2   = (TH1D*) tf_input2 -> Get("h_phi_v1_vs_y");

  // TH1D * h_input3   = (TH1D*) tf_input3 -> Get("h_phi_v1_vs_y");
  // TH1D * h_input4   = (TH1D*) tf_input4 -> Get("h_phi_v1_vs_y");
  // TH1D * h_input5   = (TH1D*) tf_input5 -> Get("h_phi_v1_vs_y");
  //
  // TH1D * h_input6   = (TH1D*) tf_input6 -> Get("h_phi_v1_vs_y");
  // TH1D * h_input7   = (TH1D*) tf_input7 -> Get("h_phi_v1_vs_y");
  // TH1D * h_input8   = (TH1D*) tf_input8 -> Get("h_phi_v1_vs_y");
  //
  // TH1D * h_input9   = (TH1D*) tf_input9 -> Get("h_phi_v1_vs_y");
  // TH1D * h_input10   = (TH1D*) tf_input10 -> Get("h_phi_v1_vs_y");

  // TH1D * h_flip_input0   = (TH1D*) tf_input0 -> Get("h_phi_v1_vs_y_flip");
  TH1D * h_flip_input1   = (TH1D*) tf_input1 -> Get("h_phi_v1_vs_y_flip");
  TH1D * h_flip_input2   = (TH1D*) tf_input2 -> Get("h_phi_v1_vs_y_flip");

  // TH1D * h_flip_input3   = (TH1D*) tf_input3 -> Get("h_phi_v1_vs_y_flip");
  // TH1D * h_flip_input4   = (TH1D*) tf_input4 -> Get("h_phi_v1_vs_y_flip");
  // TH1D * h_flip_input5   = (TH1D*) tf_input5 -> Get("h_phi_v1_vs_y_flip");
  //
  // TH1D * h_flip_input6   = (TH1D*) tf_input6 -> Get("h_phi_v1_vs_y_flip");
  // TH1D * h_flip_input7   = (TH1D*) tf_input7 -> Get("h_phi_v1_vs_y_flip");
  // TH1D * h_flip_input8   = (TH1D*) tf_input8 -> Get("h_phi_v1_vs_y_flip");
  //
  // TH1D * h_flip_input9   = (TH1D*) tf_input9 -> Get("h_phi_v1_vs_y_flip");
  // TH1D * h_flip_input10   = (TH1D*) tf_input10 -> Get("h_phi_v1_vs_y_flip");

  // TF1  *f0 = new TF1("f0","[0] + [1]*x",-1.53,1.53);
  TF1  *f1 = new TF1("f1","[0] + [1]*x",-1.53,1.53);
  TF1  *f2 = new TF1("f2","[0] + [1]*x",-1.53,1.53);
  // TF1  *f3 = new TF1("f3","[0] + [1]*x",-1.53,1.53);
  // TF1  *f4 = new TF1("f4","[0] + [1]*x",-1.53,1.53);
  // TF1  *f5 = new TF1("f5","[0] + [1]*x",-1.53,1.53);
  // TF1  *f6 = new TF1("f6","[0] + [1]*x",-1.53,1.53);
  // TF1  *f7 = new TF1("f7","[0] + [1]*x",-1.53,1.53);
  // TF1  *f8 = new TF1("f8","[0] + [1]*x",-1.53,1.53);
  // TF1  *f9 = new TF1("f9","[0] + [1]*x",-1.53,1.53);
  // TF1  *f10 = new TF1("f10","[0] + [1]*x",-1.53,1.53);

  // TF1  *f_gaus1 = new TF1("f_gaus1","gaus",-0.3,0.3);
  // TF1  *f_gaus2 = new TF1("f_gaus2","gaus",-0.3,0.3);
  // TF1  *f_gaus3 = new TF1("f_gaus3","gaus",-0.3,0.3);
  // TF1  *f_gaus4 = new TF1("f_gaus4","gaus",-0.3,0.3);

  TProfile  *TP_profile1 = new TProfile("TP_profile1","Systematic Error of v1(y) in rapidity bin [-1.5, -1.0]",1,0,1);
  TProfile  *TP_profile2 = new TProfile("TP_profile2","Systematic Error of v1(y) in rapidity bin [-1.0, -0.5]",1,0,1);
  TProfile  *TP_profile3 = new TProfile("TP_profile3","Systematic Error of v1(y) in rapidity bin [-0.5, 0]",1,0,1);
  TProfile  *TP_profile4 = new TProfile("TP_profile4","Systematic Error of dv1/dy",1,0,1);


  // h_flip_input0 -> Fit(f0,"E","R",-1.53,1.53);
  h_flip_input1 -> Fit(f1,"E","R",-1.53,1.53);
  h_flip_input2 -> Fit(f2,"E","R",-1.53,1.53);
  // h_flip_input3 -> Fit(f3,"E","R",-1.53,1.53);
  // h_flip_input4 -> Fit(f4,"E","R",-1.53,1.53);
  // h_flip_input5 -> Fit(f5,"E","R",-1.53,1.53);
  // h_flip_input6 -> Fit(f6,"E","R",-1.53,1.53);
  // h_flip_input7 -> Fit(f7,"E","R",-1.53,1.53);
  // h_flip_input8 -> Fit(f8,"E","R",-1.53,1.53);
  // h_flip_input9 -> Fit(f9,"E","R",-1.53,1.53);
  // h_flip_input10 -> Fit(f10,"E","R",-1.53,1.53);

  TH1D * h_bin4   = new TH1D("h_bin4","Systematic Error in rapidity bin (-1.5, -1.0 )",10000000,-0.3,0.3);
  h_bin4->GetXaxis()->SetTitle("v1");

  TH1D * h_bin5   = new TH1D("h_bin5","Systematic Error in rapidity bin (-1.0, -0.5 )",10000000,-0.3,0.3);
  h_bin5->GetXaxis()->SetTitle("v1");

  TH1D * h_bin6   = new TH1D("h_bin6","Systematic Error in rapidity bin (-0.5, 0 )",10000000,-0.3,0.3);
  h_bin6->GetXaxis()->SetTitle("v1");

  TH1D * h_slope   = new TH1D("h_slope","Systematic Error of slope dv_{1}/dy",10000000,-0.3,0.3);
  h_slope->GetXaxis()->SetTitle("dv1/dy");


  // h_bin4->Fill((h_input0->GetBinContent(4)));
  h_bin4->Fill((h_input1->GetBinContent(4)));
  h_bin4->Fill((h_input2->GetBinContent(4)));

  // h_bin4->Fill((h_input3->GetBinContent(4)));
  // h_bin4->Fill((h_input4->GetBinContent(4)));
  // h_bin4->Fill((h_input5->GetBinContent(4)));
  //
  // h_bin4->Fill((h_input6->GetBinContent(4)));
  // h_bin4->Fill((h_input7->GetBinContent(4)));
  // h_bin4->Fill((h_input8->GetBinContent(4)));
  //
  // h_bin4->Fill((h_input9->GetBinContent(4)));
  // h_bin4->Fill((h_input10->GetBinContent(4)));

// TProfile1
  // TP_profile1->Fill(0.5,(h_input0->GetBinContent(4)));
  TP_profile1->Fill(0.5,(h_input1->GetBinContent(4)));
  TP_profile1->Fill(0.5,(h_input2->GetBinContent(4)));

  // TP_profile1->Fill(0.5,(h_input3->GetBinContent(4)));
  // TP_profile1->Fill(0.5,(h_input4->GetBinContent(4)));
  // TP_profile1->Fill(0.5,(h_input5->GetBinContent(4)));
  //
  // TP_profile1->Fill(0.5,(h_input6->GetBinContent(4)));
  // TP_profile1->Fill(0.5,(h_input7->GetBinContent(4)));
  // TP_profile1->Fill(0.5,(h_input8->GetBinContent(4)));
  //
  // TP_profile1->Fill(0.5,(h_input9->GetBinContent(4)));
  // TP_profile1->Fill(0.5,(h_input10->GetBinContent(4)));

  // h_bin5->Fill((h_input0->GetBinContent(5)));
  h_bin5->Fill((h_input1->GetBinContent(5)));
  h_bin5->Fill((h_input2->GetBinContent(5)));

  // h_bin5->Fill((h_input3->GetBinContent(5)));
  // h_bin5->Fill((h_input4->GetBinContent(5)));
  // h_bin5->Fill((h_input5->GetBinContent(5)));
  //
  // h_bin5->Fill((h_input6->GetBinContent(5)));
  // h_bin5->Fill((h_input7->GetBinContent(5)));
  // h_bin5->Fill((h_input8->GetBinContent(5)));
  //
  // h_bin5->Fill((h_input9->GetBinContent(5)));
  // h_bin5->Fill((h_input10->GetBinContent(5)));

//TProfile2
  // TP_profile2->Fill(0.5,(h_input0->GetBinContent(5)));
  TP_profile2->Fill(0.5,(h_input1->GetBinContent(5)));
  TP_profile2->Fill(0.5,(h_input2->GetBinContent(5)));

  // TP_profile2->Fill(0.5,(h_input3->GetBinContent(5)));
  // TP_profile2->Fill(0.5,(h_input4->GetBinContent(5)));
  // TP_profile2->Fill(0.5,(h_input5->GetBinContent(5)));
  //
  // TP_profile2->Fill(0.5,(h_input6->GetBinContent(5)));
  // TP_profile2->Fill(0.5,(h_input7->GetBinContent(5)));
  // TP_profile2->Fill(0.5,(h_input8->GetBinContent(5)));
  //
  // TP_profile2->Fill(0.5,(h_input9->GetBinContent(5)));
  // TP_profile2->Fill(0.5,(h_input10->GetBinContent(5)));

  // h_bin6->Fill((h_input0->GetBinContent(6)));
  h_bin6->Fill((h_input1->GetBinContent(6)));
  h_bin6->Fill((h_input2->GetBinContent(6)));

  // h_bin6->Fill((h_input3->GetBinContent(6)));
  // h_bin6->Fill((h_input4->GetBinContent(6)));
  // h_bin6->Fill((h_input5->GetBinContent(6)));
  //
  // h_bin6->Fill((h_input6->GetBinContent(6)));
  // h_bin6->Fill((h_input7->GetBinContent(6)));
  // h_bin6->Fill((h_input8->GetBinContent(6)));
  //
  // h_bin6->Fill((h_input9->GetBinContent(6)));
  // h_bin6->Fill((h_input10->GetBinContent(6)));

//TProfile3
  // TP_profile3->Fill(0.5,(h_input0->GetBinContent(6)));
  TP_profile3->Fill(0.5,(h_input1->GetBinContent(6)));
  TP_profile3->Fill(0.5,(h_input2->GetBinContent(6)));

  // TP_profile3->Fill(0.5,(h_input3->GetBinContent(6)));
  // TP_profile3->Fill(0.5,(h_input4->GetBinContent(6)));
  // TP_profile3->Fill(0.5,(h_input5->GetBinContent(6)));
  //
  // TP_profile3->Fill(0.5,(h_input6->GetBinContent(6)));
  // TP_profile3->Fill(0.5,(h_input7->GetBinContent(6)));
  // TP_profile3->Fill(0.5,(h_input8->GetBinContent(6)));
  //
  // TP_profile3->Fill(0.5,(h_input9->GetBinContent(6)));
  // TP_profile3->Fill(0.5,(h_input10->GetBinContent(6)));

  // h_slope->Fill((f0->GetParameter(1)));
  h_slope->Fill((f1->GetParameter(1)));
  h_slope->Fill((f2->GetParameter(1)));
  // h_slope->Fill((f3->GetParameter(1)));
  // h_slope->Fill((f4->GetParameter(1)));
  // h_slope->Fill((f5->GetParameter(1)));
  // h_slope->Fill((f6->GetParameter(1)));
  // h_slope->Fill((f7->GetParameter(1)));
  // h_slope->Fill((f8->GetParameter(1)));
  // h_slope->Fill((f9->GetParameter(1)));
  // h_slope->Fill((f10->GetParameter(1)));

// TProfile4
  // TP_profile4->Fill(0.5,(f0->GetParameter(1)));
  TP_profile4->Fill(0.5,(f1->GetParameter(1)));
  TP_profile4->Fill(0.5,(f2->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f3->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f4->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f5->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f6->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f7->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f8->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f9->GetParameter(1)));
  // TP_profile4->Fill(0.5,(f10->GetParameter(1)));

  // h_bin4->Fit("f_gaus1");
  // h_bin5->Fit("f_gaus2");
  // h_bin6->Fit("f_gaus3");
  // h_slope->Fit("f_gaus4");

// TPaveText on the final result
  TPaveText * ptext1 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  Double_t d_mean_bin4 = TP_profile1        -> GetBinContent(1);
  Double_t d_err_on_mean_bin4 = TP_profile1 -> GetBinError(1);
  ptext1 -> AddText(Form("Mean: %.4f",d_mean_bin4));
  ptext1 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_bin4));

  TPaveText * ptext2 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  Double_t d_mean_bin5 = TP_profile2        -> GetBinContent(1);
  Double_t d_err_on_mean_bin5 = TP_profile2 -> GetBinError(1);
  ptext2 -> AddText(Form("Mean: %.4f",d_mean_bin5));
  ptext2 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_bin5));

  TPaveText * ptext3 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  Double_t d_mean_bin6 = TP_profile3        -> GetBinContent(1);
  Double_t d_err_on_mean_bin6 = TP_profile3 -> GetBinError(1);
  ptext3 -> AddText(Form("Mean: %.4f",d_mean_bin6));
  ptext3 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_bin6));

  TPaveText * ptext4 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  Double_t d_mean_slope = TP_profile4        -> GetBinContent(1);
  Double_t d_err_on_mean_slope = TP_profile4 -> GetBinError(1);
  ptext4 -> AddText(Form("Mean: %.4f",d_mean_slope));
  ptext4 -> AddText(Form("Error on mean: %.4f",d_err_on_mean_slope));

  TCanvas * TC_inv1 = new TCanvas("TC_inv1","TC_inv1",1280,720);
  TC_inv1->cd();
  h_bin4->Draw();
  ptext1->Draw("same");

  TCanvas * TC_inv2 = new TCanvas("TC_inv2","TC_inv2",1280,720);
  TC_inv2->cd();
  h_bin5->Draw();
  ptext2->Draw("same");

  TCanvas * TC_inv3 = new TCanvas("TC_inv3","TC_inv3",1280,720);
  TC_inv3->cd();
  h_bin6->Draw();
  ptext3->Draw("same");

  TCanvas * TC_inv4 = new TCanvas("TC_inv4","TC_inv4",1280,720);
  TC_inv4->cd();
  h_slope->Draw();
  ptext4->Draw("same");


  // TFile * outFile = new TFile("sys_err_21cuts_merged_1.root","RECREATE" );
  TFile * outFile = new TFile("sys_err_21cuts_eff_corr_EP.root","RECREATE" );

  outFile->cd();
  h_bin4->Write();
  h_bin5->Write();
  h_bin6->Write();
  h_slope->Write();

  TP_profile1->Write();
  TP_profile2->Write();
  TP_profile3->Write();
  TP_profile4->Write();

  TC_inv1->Write();
  TC_inv2->Write();
  TC_inv3->Write();
  TC_inv4->Write();

  // h_input0 ->Write();
  h_input1 ->Write();
  h_input2 ->Write();

  // h_input3 ->Write();
  // h_input4 ->Write();
  // h_input5 ->Write();
  //
  // h_input6 ->Write();
  // h_input7 ->Write();
  // h_input8 ->Write();
  //
  // h_input9 ->Write();
  // h_input10 ->Write();

  // h_flip_input0 ->Write();
  h_flip_input1 ->Write();
  h_flip_input2 ->Write();
  // h_flip_input3 ->Write();
  // h_flip_input4 ->Write();
  // h_flip_input5 ->Write();
  // h_flip_input6 ->Write();
  // h_flip_input7 ->Write();
  // h_flip_input8 ->Write();
  // h_flip_input9 ->Write();
  // h_flip_input10 ->Write();
}
