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

#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TRandom3.h"
#include "TDatime.h"

#include "Math/MinimizerOptions.h"

using namespace std;

void phiFlowOverlap()
{
  TFile * tf_const  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_flow_fitting_out_const.root","READ");
  TFile * tf_pol1   = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_flow_fitting_out_pol1.root","READ");
  TFile * tf_pol2   = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_flow_fitting_out_pol2.root","READ");
  TFile * tf_pol3   = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_flow_fitting_out_pol3.root","READ");

  TH1D  * h_phi_v1_vs_y_flip_const = (TH1D*) tf_const -> Get("h_phi_v1_vs_y_flip");
  TH1D  * h_phi_v1_vs_y_flip_pol1  = (TH1D*) tf_pol1  -> Get("h_phi_v1_vs_y_flip");
  TH1D  * h_phi_v1_vs_y_flip_pol2  = (TH1D*) tf_pol2  -> Get("h_phi_v1_vs_y_flip");
  TH1D  * h_phi_v1_vs_y_flip_pol3  = (TH1D*) tf_pol3  -> Get("h_phi_v1_vs_y_flip");

  TCanvas * TC_v1_y_flip_overlap = new TCanvas("TC_v1_y_flip_overlap","TC_v1_y_flip_overlap",1280,720);
  TC_v1_y_flip_overlap->cd();
  h_phi_v1_vs_y_flip_const->SetMarkerColor(kBlack);
  h_phi_v1_vs_y_flip_pol1 ->SetMarkerColor(kBlue);
  h_phi_v1_vs_y_flip_pol2 ->SetMarkerColor(kGreen);
  h_phi_v1_vs_y_flip_pol3 ->SetMarkerColor(kRed);

  h_phi_v1_vs_y_flip_const->Draw();
  h_phi_v1_vs_y_flip_pol1 ->Draw("same");
  h_phi_v1_vs_y_flip_pol2 ->Draw("same");
  h_phi_v1_vs_y_flip_pol3 ->Draw("same");

  TF1 * tf1_pol_const = new TF1("tf1_pol_const","[0] + [1]*x",-1.0,1.0);
  TF1 * tf1_pol_pol1  = new TF1("tf1_pol_pol1" ,"[0] + [1]*x",-1.0,1.0);
  TF1 * tf1_pol_pol2  = new TF1("tf1_pol_pol2" ,"[0] + [1]*x",-1.0,1.0);
  TF1 * tf1_pol_pol3  = new TF1("tf1_pol_pol3" ,"[0] + [1]*x",-1.0,1.0);

  tf1_pol_const->SetLineColor(kBlack);
  tf1_pol_pol1->SetLineColor(kBlue);
  tf1_pol_pol2->SetLineColor(kGreen);
  tf1_pol_pol3->SetLineColor(kRed);

  h_phi_v1_vs_y_flip_const->Fit(tf1_pol_const,"E","R",-1.0,1.0);
  h_phi_v1_vs_y_flip_pol1 ->Fit(tf1_pol_pol1,"E","R",-1.0,1.0);
  h_phi_v1_vs_y_flip_pol2 ->Fit(tf1_pol_pol2,"E","R",-1.0,1.0);
  h_phi_v1_vs_y_flip_pol3 ->Fit(tf1_pol_pol3,"E","R",-1.0,1.0);
  // gStyle->SetLegendTextSize(0.);

  TLegend * TL_legend = new TLegend(0.1,0.7,0.48,0.9);
  // gStyle->SetLegendTextSize(0);
  TL_legend->SetTextSize(0.03);
  TL_legend->AddEntry(h_phi_v1_vs_y_flip_const,"v_{1}^{Sig}(y): v_{1}^{Bg} as a constant","lep");
  TL_legend->AddEntry(tf1_pol_const,"Linear fit: v_{1}^{Bg} as a constant","l");

  TL_legend->AddEntry(h_phi_v1_vs_y_flip_pol1,"v_{1}^{Sig}(y): v_{1}^{Bg} as a Pol1","lep");
  TL_legend->AddEntry(tf1_pol_pol1,"Linear fit: v_{1}^{Bg} as a Pol1","l");

  TL_legend->AddEntry(h_phi_v1_vs_y_flip_pol2,"v_{1}^{Sig}(y): v_{1}^{Bg} as a Pol2","lep");
  TL_legend->AddEntry(tf1_pol_pol2,"Linear fit: v_{1}^{Bg} as a Pol2","l");

  TL_legend->AddEntry(h_phi_v1_vs_y_flip_pol3,"v_{1}^{Sig}(y): v_{1}^{Bg} as a Pol3","lep");
  TL_legend->AddEntry(tf1_pol_pol3,"Linear fit: v_{1}^{Bg} as a Pol3","l");
  TL_legend->Draw("same");

}
