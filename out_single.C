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

void out_single()
{
  TFile * tf_evt_in  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/kaonTTree/merged_F6A9F929788E732058DC337FA808DEE6.picoDst.result.root","READ");

  // mixed event invM
  TFile * tf_evt_mx_in  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_systematic_error_analysis_invM/phi_meson_mixed_events_invariant_mass_primary.root","READ");

  TH1D * h_evt   = (TH1D*) tf_evt_in -> Get("h_evt");
  double d_N_evt    = h_evt -> Integral();

  // normal event invM
  // TFile * tf_total = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_systematic_error_analysis_invM/","READ");
  TFile * tf_total = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_test/test_all.readKTree.result.root","READ");

  TH1D * h_TPC2_TOF3_invM = (TH1D*) tf_total -> Get("h_TPC2_TOF3_invM");
  h_TPC2_TOF3_invM->Rebin();
  // TFile * tf_mx = new TFile("h_mx_TPC2_TOF3_invM_out_tot.root","READ");
  // TH1D  * h_mx_TPC2_TOF3_invM = (TH1D*) tf_mx -> Get("h_mx_TPC2_TOF3_invM");
  TH1D  * h_mx_TPC2_TOF3_invM = (TH1D*) tf_evt_mx_in -> Get("h_mx_TPC2_TOF3_invM");
  h_mx_TPC2_TOF3_invM->Rebin();



  // TFile * tf_total_y_bin2 = new TFile("h_TPC2_TOF3_invM_y_bin2_out_tot.root","READ");
  TH1D * h_TPC2_TOF3_invM_y_bin2 = (TH1D*) tf_total -> Get("h_TPC2_TOF3_invM_y_bin2");
  h_TPC2_TOF3_invM_y_bin2->Rebin();

  // TFile * tf_mx_y_bin2 = new TFile("h_mx_TPC2_TOF3_invM_y_bin2_out_tot.root","READ");
  // TH1D  * h_mx_TPC2_TOF3_invM_y_bin2 = (TH1D*) tf_mx_y_bin2 -> Get("h_mx_TPC2_TOF3_invM_y_bin2");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin2 = (TH1D*) tf_evt_mx_in -> Get("h_mx_TPC2_TOF3_invM_y_bin2");
  h_mx_TPC2_TOF3_invM_y_bin2->Rebin();


  // TFile * tf_total_y_bin3 = new TFile("h_TPC2_TOF3_invM_y_bin3_out_tot.root","READ");
  TH1D * h_TPC2_TOF3_invM_y_bin3 = (TH1D*) tf_total -> Get("h_TPC2_TOF3_invM_y_bin3");
  h_TPC2_TOF3_invM_y_bin3->Rebin();

  // TFile * tf_mx_y_bin3 = new TFile("h_mx_TPC2_TOF3_invM_y_bin3_out_tot.root","READ");
  // TH1D  * h_mx_TPC2_TOF3_invM_y_bin3 = (TH1D*) tf_mx_y_bin3 -> Get("h_mx_TPC2_TOF3_invM_y_bin3");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin3 = (TH1D*) tf_evt_mx_in -> Get("h_mx_TPC2_TOF3_invM_y_bin3");
  h_mx_TPC2_TOF3_invM_y_bin3->Rebin();


  // TFile * tf_total_y_bin4 = new TFile("h_TPC2_TOF3_invM_y_bin4_out_tot.root","READ");
  TH1D * h_TPC2_TOF3_invM_y_bin4 = (TH1D*) tf_total -> Get("h_TPC2_TOF3_invM_y_bin4");
  h_TPC2_TOF3_invM_y_bin4->Rebin();

  // TFile * tf_mx_y_bin4 = new TFile("h_mx_TPC2_TOF3_invM_y_bin4_out_tot.root","READ");
  // TH1D  * h_mx_TPC2_TOF3_invM_y_bin4 = (TH1D*) tf_mx_y_bin4 -> Get("h_mx_TPC2_TOF3_invM_y_bin4");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin4 = (TH1D*) tf_evt_mx_in -> Get("h_mx_TPC2_TOF3_invM_y_bin4");
  h_mx_TPC2_TOF3_invM_y_bin4->Rebin();


  //--- norm
  double a_d_int_range[4] =
    {
      0.99,
      1.008,
      1.04,//1.06,
      1.09
    };//double a_d_int_range[4] =

  int a_iBin_range[4];
  for(int ijk = 0; ijk < 4; ijk++) a_iBin_range[ijk] =  h_mx_TPC2_TOF3_invM -> FindFixBin(a_d_int_range[ijk]);
  gStyle->SetOptDate(0);
  TCanvas * TC_inv = new TCanvas("TC_inv","TC_inv",1280,720);
  h_TPC2_TOF3_invM   -> SetXTitle("M_{K+K-}");
  h_TPC2_TOF3_invM   -> SetYTitle("Counts");
  TC_inv->cd();
  // h_TPC2_TOF3_invM    ->Rebin();
  // h_mx_TPC2_TOF3_invM ->Rebin();
  h_TPC2_TOF3_invM    -> Draw();
  h_mx_TPC2_TOF3_invM -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM -> SetFillStyle(3001);

  TCanvas * TC_inv_y_bin2 = new TCanvas("TC_inv_y_bin2","TC_inv_y_bin2",1280,720);
  h_TPC2_TOF3_invM_y_bin2   -> SetXTitle("M_{K+K-}");
  h_TPC2_TOF3_invM_y_bin2   -> SetYTitle("Counts");
  TC_inv_y_bin2->cd();
  h_TPC2_TOF3_invM_y_bin2    -> Draw();
  h_mx_TPC2_TOF3_invM_y_bin2 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin2 -> SetFillStyle(3001);

  TCanvas * TC_inv_y_bin3 = new TCanvas("TC_inv_y_bin3","TC_inv_y_bin3",1280,720);
  h_TPC2_TOF3_invM_y_bin3   -> SetXTitle("M_{K+K-}");
  h_TPC2_TOF3_invM_y_bin3   -> SetYTitle("Counts");
  TC_inv_y_bin3->cd();
  h_TPC2_TOF3_invM_y_bin3    -> Draw();
  h_mx_TPC2_TOF3_invM_y_bin3 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin3 -> SetFillStyle(3001);

  TCanvas * TC_inv_y_bin4 = new TCanvas("TC_inv_y_bin4","TC_inv_y_bin4",1280,720);
  h_TPC2_TOF3_invM_y_bin4   -> SetXTitle("M_{K+K-}");
  h_TPC2_TOF3_invM_y_bin4   -> SetYTitle("Counts");
  TC_inv_y_bin4->cd();
  // h_TPC2_TOF3_invM_y_bin4 -> Rebin();
  h_TPC2_TOF3_invM_y_bin4    -> Draw();
  h_mx_TPC2_TOF3_invM_y_bin4 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin4 -> SetFillStyle(3001);

  //  h_mx_TPC2_TOF3_invM -> Draw("sames");

  //right bg norm
  double d_r_area     = a_d_int_range[3]-a_d_int_range[2];
  double d_r_prim_int = h_TPC2_TOF3_invM -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_mx_int   = h_mx_TPC2_TOF3_invM -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_norm     = (d_r_mx_int != 0.0) ? d_r_prim_int/d_r_mx_int : 1.0;

  double d_r_prim_int_y_bin2 = h_TPC2_TOF3_invM_y_bin2 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_mx_int_y_bin2   = h_mx_TPC2_TOF3_invM_y_bin2 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_norm_y_bin2     = (d_r_mx_int_y_bin2 != 0.0) ? d_r_prim_int_y_bin2/d_r_mx_int_y_bin2 : 1.0;

  double d_r_prim_int_y_bin3 = h_TPC2_TOF3_invM_y_bin3 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_mx_int_y_bin3   = h_mx_TPC2_TOF3_invM_y_bin3 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_norm_y_bin3     = (d_r_mx_int_y_bin3 != 0.0) ? d_r_prim_int_y_bin3/d_r_mx_int_y_bin3 : 1.0;

  double d_r_prim_int_y_bin4 = h_TPC2_TOF3_invM_y_bin4 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_mx_int_y_bin4   = h_mx_TPC2_TOF3_invM_y_bin4 -> Integral(a_iBin_range[2],a_iBin_range[3]);
  double d_r_norm_y_bin4     = (d_r_mx_int_y_bin4 != 0.0) ? d_r_prim_int_y_bin4/d_r_mx_int_y_bin4 : 1.0;

  cout<<" R: "<<d_r_area<<" : "<<d_r_prim_int<<" : "<<d_r_mx_int<<endl;
  //left bg norm
  double d_l_area =  a_d_int_range[1]-a_d_int_range[0];
  double d_l_prim_int = h_TPC2_TOF3_invM -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_mx_int   = h_mx_TPC2_TOF3_invM -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_norm     = (d_l_mx_int != 0.0) ? d_l_prim_int/d_l_mx_int : 1.0;

  double d_l_prim_int_y_bin2 = h_TPC2_TOF3_invM_y_bin2 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_mx_int_y_bin2   = h_mx_TPC2_TOF3_invM_y_bin2 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_norm_y_bin2     = (d_l_mx_int != 0.0) ? d_l_prim_int_y_bin2/d_l_mx_int_y_bin2 : 1.0;

  double d_l_prim_int_y_bin3 = h_TPC2_TOF3_invM_y_bin3 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_mx_int_y_bin3   = h_mx_TPC2_TOF3_invM_y_bin3 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_norm_y_bin3     = (d_l_mx_int != 0.0) ? d_l_prim_int_y_bin3/d_l_mx_int_y_bin3 : 1.0;

  double d_l_prim_int_y_bin4 = h_TPC2_TOF3_invM_y_bin4 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_mx_int_y_bin4   = h_mx_TPC2_TOF3_invM_y_bin4 -> Integral(a_iBin_range[0],a_iBin_range[1]);
  double d_l_norm_y_bin4     = (d_l_mx_int != 0.0) ? d_l_prim_int_y_bin4/d_l_mx_int_y_bin4 : 1.0;



  cout<<" L: "<<d_l_area<<" : "<<d_l_prim_int<<" : "<<d_l_mx_int<<endl;
  cout<<" l: "<<d_l_norm<<" r: "<<d_r_norm<<endl;
  double d_norm = ((d_r_norm/d_r_area) + (d_l_norm/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));

  double d_norm_y_bin2 = ((d_r_norm_y_bin2/d_r_area) + (d_l_norm_y_bin2/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
  double d_norm_y_bin3 = ((d_r_norm_y_bin3/d_r_area) + (d_l_norm_y_bin3/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));
  double d_norm_y_bin4 = ((d_r_norm_y_bin4/d_r_area) + (d_l_norm_y_bin4/d_l_area))/( (1.0/d_r_area) + (1.0/d_l_area));

  cout<<" d_norm = "<<d_norm<<endl;
  cout<<" d_norm_y_bin2 = "<<d_norm_y_bin2<<endl;
  cout<<" d_norm_y_bin3 = "<<d_norm_y_bin3<<endl;
  cout<<" d_norm_y_bin4 = "<<d_norm_y_bin4<<endl;


  h_mx_TPC2_TOF3_invM -> Sumw2();
  h_mx_TPC2_TOF3_invM -> Scale(d_norm);

  h_mx_TPC2_TOF3_invM_y_bin2 -> Sumw2();
  h_mx_TPC2_TOF3_invM_y_bin2 -> Scale(d_norm_y_bin2);

  h_mx_TPC2_TOF3_invM_y_bin3 -> Sumw2();
  h_mx_TPC2_TOF3_invM_y_bin3 -> Scale(d_norm_y_bin3);

  h_mx_TPC2_TOF3_invM_y_bin4 -> Sumw2();
  h_mx_TPC2_TOF3_invM_y_bin4 -> Scale(d_norm_y_bin4);

  // h_mx_TPC2_TOF3_invM_y_bin4 -> Rebin();


  TC_inv->cd();
  h_mx_TPC2_TOF3_invM -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM -> Draw("HISTsames");

  TC_inv_y_bin2->cd();
  h_mx_TPC2_TOF3_invM_y_bin2 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin2 -> Draw("HISTsames");

  TC_inv_y_bin3->cd();
  h_mx_TPC2_TOF3_invM_y_bin3 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin3 -> Draw("HISTsames");

  TC_inv_y_bin4->cd();
  h_mx_TPC2_TOF3_invM_y_bin4 -> SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin4 -> Draw("HISTsames");

  TLine * a_TL_int[4];
  for(int k = 0; k < 4; k++)
    {
      a_TL_int[k] = new TLine(a_d_int_range[k],0.0001,a_d_int_range[k],1.2*h_mx_TPC2_TOF3_invM -> GetMaximum());
      a_TL_int[k] -> SetLineColor(kRed);
      TC_inv->cd();
      a_TL_int[k] -> Draw();

      TC_inv_y_bin2->cd();
      a_TL_int[k] -> Draw();

      TC_inv_y_bin3->cd();
      a_TL_int[k] -> Draw();

      TC_inv_y_bin4->cd();
      a_TL_int[k] -> Draw();

    }//for(int k = 0; k < 4; k++)
  //--- END norm lines

  TH1D * h_sub = (TH1D*) h_TPC2_TOF3_invM -> Clone("h_sub");
  h_sub -> Reset();
  h_sub -> Sumw2();

  TH1D * h_sub_y_bin2= (TH1D*) h_TPC2_TOF3_invM_y_bin2 -> Clone("h_sub_y_bin2");
  h_sub_y_bin2 -> Reset();
  h_sub_y_bin2 -> Sumw2();

  TH1D * h_sub_y_bin3 = (TH1D*) h_TPC2_TOF3_invM_y_bin3 -> Clone("h_sub_y_bin3");
  h_sub_y_bin3 -> Reset();
  h_sub_y_bin3 -> Sumw2();

  TH1D * h_sub_y_bin4 = (TH1D*) h_TPC2_TOF3_invM_y_bin4 -> Clone("h_sub_y_bin4");
  h_sub_y_bin4 -> Reset();
  h_sub_y_bin4 -> Sumw2();
  for(int ijk = 1; ijk < h_sub->GetNbinsX()+1; ijk++)
    {
      double d_center   = h_TPC2_TOF3_invM -> GetBinCenter(ijk);
      double d_prim     = h_TPC2_TOF3_invM -> GetBinContent(ijk);
      double d_prim_err = h_TPC2_TOF3_invM -> GetBinError(ijk);

      double d_mx       = h_mx_TPC2_TOF3_invM -> GetBinContent(ijk);
      double d_mx_err   = h_mx_TPC2_TOF3_invM -> GetBinError(ijk);

      double d_sub      = d_prim - d_mx;
      double d_sub_err  = sqrt(d_prim_err*d_prim_err+d_mx_err*d_mx_err);

      //scale by 2 beacuse we are going to rebin by 2 later
      h_sub -> SetBinContent(ijk,d_sub);///2.);
      h_sub -> SetBinError(ijk,d_sub_err);///2.);

      double d_center_y_bin2   = h_TPC2_TOF3_invM_y_bin2 -> GetBinCenter(ijk);
      double d_prim_y_bin2     = h_TPC2_TOF3_invM_y_bin2 -> GetBinContent(ijk);
      double d_prim_err_y_bin2 = h_TPC2_TOF3_invM_y_bin2 -> GetBinError(ijk);

      double d_mx_y_bin2       = h_mx_TPC2_TOF3_invM_y_bin2 -> GetBinContent(ijk);
      double d_mx_err_y_bin2   = h_mx_TPC2_TOF3_invM_y_bin2 -> GetBinError(ijk);

      double d_sub_y_bin2      = d_prim_y_bin2 - d_mx_y_bin2;
      double d_sub_err_y_bin2  = sqrt(d_prim_err_y_bin2*d_prim_err_y_bin2+d_mx_err_y_bin2*d_mx_err_y_bin2);

      //scale by 2 beacuse we are going to rebin by 2 later
      h_sub_y_bin2 -> SetBinContent(ijk,d_sub_y_bin2);///2.);
      h_sub_y_bin2 -> SetBinError(ijk,d_sub_err_y_bin2);///2.);

      double d_center_y_bin3   = h_TPC2_TOF3_invM_y_bin3 -> GetBinCenter(ijk);
      double d_prim_y_bin3     = h_TPC2_TOF3_invM_y_bin3 -> GetBinContent(ijk);
      double d_prim_err_y_bin3 = h_TPC2_TOF3_invM_y_bin3 -> GetBinError(ijk);

      double d_mx_y_bin3       = h_mx_TPC2_TOF3_invM_y_bin3 -> GetBinContent(ijk);
      double d_mx_err_y_bin3   = h_mx_TPC2_TOF3_invM_y_bin3 -> GetBinError(ijk);

      double d_sub_y_bin3      = d_prim_y_bin3 - d_mx_y_bin3;
      double d_sub_err_y_bin3  = sqrt(d_prim_err_y_bin3*d_prim_err_y_bin3+d_mx_err_y_bin3*d_mx_err_y_bin3);

      //scale by 2 beacuse we are going to rebin by 2 later
      h_sub_y_bin3 -> SetBinContent(ijk,d_sub_y_bin3);///2.);
      h_sub_y_bin3 -> SetBinError(ijk,d_sub_err_y_bin3);///2.);

      double d_center_y_bin4   = h_TPC2_TOF3_invM_y_bin4 -> GetBinCenter(ijk);
      double d_prim_y_bin4     = h_TPC2_TOF3_invM_y_bin4 -> GetBinContent(ijk);
      double d_prim_err_y_bin4 = h_TPC2_TOF3_invM_y_bin4 -> GetBinError(ijk);

      double d_mx_y_bin4       = h_mx_TPC2_TOF3_invM_y_bin4 -> GetBinContent(ijk);
      double d_mx_err_y_bin4   = h_mx_TPC2_TOF3_invM_y_bin4 -> GetBinError(ijk);

      double d_sub_y_bin4      = d_prim_y_bin4 - d_mx_y_bin4;
      double d_sub_err_y_bin4  = sqrt(d_prim_err_y_bin4*d_prim_err_y_bin4+d_mx_err_y_bin4*d_mx_err_y_bin4);

      //scale by 2 beacuse we are going to rebin by 2 later
      h_sub_y_bin4 -> SetBinContent(ijk,d_sub_y_bin4);///2.);
      h_sub_y_bin4 -> SetBinError(ijk,d_sub_err_y_bin4);///2.);

    }//for(int ijk = 1; ijk < h_sub->GetNbinsX()+1; ijk++)

  TCanvas * TC_sub = new TCanvas("TC_sub","TC_sub",1280,720);
  // h_sub -> Rebin();
  TC_sub->cd();
  h_sub -> Draw("E");
  gStyle->SetOptFit(1111);

  TCanvas * TC_sub_y_bin2 = new TCanvas("TC_sub_y_bin2","TC_sub_y_bin2",1280,720);
  // h_sub_y_bin2 -> Rebin();
  TC_sub_y_bin2->cd();
  h_sub_y_bin2 -> Draw("E");
  gStyle->SetOptFit(1111);

  TCanvas * TC_sub_y_bin3 = new TCanvas("TC_sub_y_bin3","TC_sub_y_bin3",1280,720);
  // h_sub_y_bin3 -> Rebin();
  TC_sub_y_bin3->cd();
  h_sub_y_bin3 -> Draw("E");
  gStyle->SetOptFit(1111);

  TCanvas * TC_sub_y_bin4 = new TCanvas("TC_sub_y_bin4","TC_sub_y_bin4",1280,720);
  // h_sub_y_bin4 -> Rebin();
  TC_sub_y_bin4->cd();
  h_sub_y_bin4 -> Draw("E");
  gStyle->SetOptFit(1111);

  //Fit function
  //fit to Gauss plus constant
  TFormula * tform_single = new TFormula("tform_single","gaus(0)+[3]");
  TF1 * tf1_phi = new TF1("polygauss_single",tform_single->GetExpFormula(),1.,1.06);

  //fit to a simple gauss first to seed
  TF1 * tf1_gauss = new TF1("tf1_gauss","gaus",0.9,1.1);
  TC_sub->cd();
  h_sub -> Fit(tf1_gauss,"0","R",1.01,1.03);

  double d_prefit_p0    = tf1_gauss -> GetParameter(0);
  double d_prefit_mean  = tf1_gauss -> GetParameter(1);
  double d_prefit_sigma = tf1_gauss -> GetParameter(2);

  cout << "d_prefit_p0 = "<< d_prefit_p0 << endl;
  cout << "d_prefit_mean = "<< d_prefit_mean << endl;
  cout << "d_prefit_sigma = "<< d_prefit_sigma << endl;

  tf1_phi -> SetParameter(1,d_prefit_mean);
  tf1_phi -> SetParameter(2,d_prefit_sigma);// test
  tf1_phi -> SetParLimits(1,d_prefit_mean-d_prefit_sigma,d_prefit_mean+d_prefit_sigma);
  tf1_phi -> SetParLimits(2,0.66*d_prefit_sigma,1.5*d_prefit_sigma);

  tf1_phi -> SetLineColor(kBlue);
  h_sub   -> Fit(tf1_phi,"E+","R",0.98,1.1);
  h_sub   -> SetTitle(Form("%s After BG Subtraction",h_sub->GetTitle() ));
  h_sub   -> SetXTitle("M_{K+K-}");
  h_sub   -> SetYTitle("Counts");

  double d_p0    = tf1_phi -> GetParameter(0);
  double d_mean  = tf1_phi -> GetParameter(1);
  double d_sigma = tf1_phi -> GetParameter(2);
  double d_p3    = tf1_phi -> GetParameter(3);

  cout << "d_p0 = "<< d_p0 << endl;
  cout << "d_mean = "<< d_mean << endl;
  cout << "d_sigma = "<< d_sigma << endl;
  cout << "d_p3 = "<< d_p3 << endl;

  int iBin_3sigint_low = h_sub -> FindFixBin(d_mean - (3*d_sigma));
  int iBin_3sigint_hi  = h_sub -> FindFixBin(d_mean + (3*d_sigma));
  double d_3sig_integral_error;
  double d_3sig_integral = h_sub -> IntegralAndError(iBin_3sigint_low,iBin_3sigint_hi,d_3sig_integral_error,"");
  double d_dNdy_err = 1.0/sqrt(d_3sig_integral);
  double d_Y_mT_err = 1.0/sqrt(d_3sig_integral);

  //  TPaveText * ptext = new TPaveText(0.65,0.5,1.,0.95,"NDCARC");
  TPaveText * ptext = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext -> AddText(Form("Mean: %.4f",d_mean));
  ptext -> AddText(Form("Sigma: %.4f",d_sigma));
  ptext -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral));
  ptext -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error));

  TC_sub->cd();
  ptext -> Draw("same");

  // fit BG within 3 sigma by 2nd order polynomial
  // TFormula * tform_pol = new TFormula("tform_pol","pol2(0)");
  // TF1 * tf1_bg = new TF1("polynomial_bg",tform_pol->GetExpFormula(),(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  //
  // //fit to a simple polynomial first to seed
  // TF1 * tf1_pol = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  // TC_inv->cd();
  //
  // h_mx_TPC2_TOF3_invM->Fit(tf1_pol,"E","R",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  // double d_prefit_linear_bg  = tf1_pol -> GetParameter(1);
  // double d_prefit_quardra_bg = tf1_pol -> GetParameter(2);
  //
  // tf1_bg -> SetParameter(1,d_prefit_linear_bg);
  // tf1_bg -> SetParameter(2,d_prefit_quardra_bg);// test
  // // tf1_bg -> SetParLimits(1,);
  // // tf1_bg -> SetParLimits(2,);
  //
  // tf1_bg -> SetLineColor(kBlue);
  // h_mx_TPC2_TOF3_invM   -> Fit(tf1_bg,"E+","R",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  // h_mx_TPC2_TOF3_invM   -> SetTitle(Form("%s With Fit in 3#sigma",h_mx_TPC2_TOF3_invM->GetTitle() ));
  // h_mx_TPC2_TOF3_invM   -> SetXTitle("M_{K+K-}");
  // h_mx_TPC2_TOF3_invM   -> SetYTitle("Counts");
  //
  // double d_const_bg   = tf1_bg -> GetParameter(0);
  // double d_linear_bg  = tf1_bg -> GetParameter(1);
  // double d_quardra_bg = tf1_bg -> GetParameter(2);
  // double d_shift_bg = tf1_bg -> GetParameter(3);
  //
  // cout << "d_const_bg = " << d_const_bg << endl;
  // cout << "d_linear_bg = " << d_linear_bg << endl;
  // cout << "d_quardra_bg = " << d_quardra_bg << endl;
  // cout << "d_shift_bg = " << d_shift_bg << endl;
  //
  // TPaveText * ptext_bg = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  // ptext_bg -> AddText(Form("Constant: %.4f",d_const_bg));
  // ptext_bg -> AddText(Form("Linear: %.4f",d_linear_bg));
  // ptext_bg -> AddText(Form("Qudratic: %.4f",d_quardra_bg));
  // TC_inv->cd();
  // ptext_bg -> Draw("same");

  // END fit BG within 3 sigma by 2nd order polynomial

  //Fit function
  //fit to Gauss plus constant
  TFormula * tform_single_y_bin2 = new TFormula("tform_single_y_bin2","gaus(0)+[3]");
  TF1 * tf1_phi_y_bin2 = new TF1("polygauss_single_y_bin2",tform_single_y_bin2->GetExpFormula(),1.,1.06);

  //fit to a simple gauss first to seed
  TF1 * tf1_gauss_y_bin2 = new TF1("tf1_gauss_y_bin2","gaus",0.9,1.1);
  TC_sub_y_bin2->cd();
  h_sub_y_bin2 -> Fit(tf1_gauss_y_bin2,"0","R",1.01,1.03);

  double d_prefit_p0_y_bin2    = tf1_gauss_y_bin2 -> GetParameter(0);
  double d_prefit_mean_y_bin2  = tf1_gauss_y_bin2 -> GetParameter(1);
  double d_prefit_sigma_y_bin2 = tf1_gauss_y_bin2 -> GetParameter(2);

  cout << "d_prefit_p0_y_bin2 = "<< d_prefit_p0_y_bin2 << endl;
  cout << "d_prefit_mean_y_bin2 = "<< d_prefit_mean_y_bin2 << endl;
  cout << "d_prefit_sigma_y_bin2 = "<< d_prefit_sigma_y_bin2 << endl;

  tf1_phi_y_bin2 -> SetParameter(1,d_prefit_mean_y_bin2);
  tf1_phi_y_bin2 -> SetParameter(2,d_prefit_sigma_y_bin2);
  tf1_phi_y_bin2 -> SetParLimits(1,d_prefit_mean_y_bin2-d_prefit_sigma_y_bin2,d_prefit_mean_y_bin2+d_prefit_sigma_y_bin2);
  tf1_phi_y_bin2 -> SetParLimits(2,0.66*d_prefit_sigma_y_bin2,1.5*d_prefit_sigma_y_bin2);

  tf1_phi_y_bin2 -> SetLineColor(kBlue);
  h_sub_y_bin2   -> Fit(tf1_phi_y_bin2,"E+","R",0.98,1.1);
  h_sub_y_bin2   -> SetTitle(Form("%s After BG Subtraction",h_sub_y_bin2->GetTitle() ));
  h_sub_y_bin2   -> SetXTitle("M_{K+K-}");
  h_sub_y_bin2   -> SetYTitle("Counts");

  double d_p0_y_bin2    = tf1_phi_y_bin2 -> GetParameter(0);
  double d_mean_y_bin2  = tf1_phi_y_bin2 -> GetParameter(1);
  double d_sigma_y_bin2 = tf1_phi_y_bin2 -> GetParameter(2);
  double d_p3_y_bin2    = tf1_phi_y_bin2 -> GetParameter(3);

  cout << "d_p0_y_bin2 = "<< d_p0_y_bin2 << endl;
  cout << "d_mean_y_bin2 = "<< d_mean_y_bin2 << endl;
  cout << "d_sigma_y_bin2 = "<< d_sigma_y_bin2 << endl;
  cout << "d_p3_y_bin2 = "<< d_p3_y_bin2 << endl;

  int iBin_3sigint_low_y_bin2 = h_sub_y_bin2 -> FindFixBin(d_mean_y_bin2 - (3*d_sigma_y_bin2));
  int iBin_3sigint_hi_y_bin2  = h_sub_y_bin2 -> FindFixBin(d_mean_y_bin2 + (3*d_sigma_y_bin2));
  double d_3sig_integral_error_y_bin2;

  double d_3sig_integral_y_bin2 = h_sub_y_bin2 -> IntegralAndError(iBin_3sigint_low_y_bin2,iBin_3sigint_hi_y_bin2,d_3sig_integral_error_y_bin2,"");//integralerror
  double d_dNdy_err_y_bin2 = 1.0/sqrt(d_3sig_integral_y_bin2);
  double d_Y_mT_err_y_bin2 = 1.0/sqrt(d_3sig_integral_y_bin2);

  //  TPaveText * ptext = new TPaveText(0.65,0.5,1.,0.95,"NDCARC");
  TPaveText * ptext_y_bin2 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin2 -> AddText(Form("Mean: %.4f",d_mean_y_bin2));
  ptext_y_bin2 -> AddText(Form("Sigma: %.4f",d_sigma_y_bin2));
  ptext_y_bin2 -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral_y_bin2));
  ptext_y_bin2 -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error_y_bin2));

  TC_sub_y_bin2->cd();
  ptext_y_bin2 -> Draw("same");

  //Fit function
  //fit to Gauss plus constant
  TFormula * tform_single_y_bin3 = new TFormula("tform_single_y_bin3","gaus(0)+[3]");
  TF1 * tf1_phi_y_bin3 = new TF1("polygauss_single_y_bin3",tform_single_y_bin3->GetExpFormula(),1.,1.06);

  //fit to a simple gauss first to seed
  TF1 * tf1_gauss_y_bin3 = new TF1("tf1_gauss_y_bin3","gaus",0.9,1.1);
  TC_sub_y_bin3->cd();
  h_sub_y_bin3 -> Fit(tf1_gauss_y_bin3,"0","R",1.01,1.03);

  double d_prefit_p0_y_bin3    = tf1_gauss_y_bin3 -> GetParameter(0);
  double d_prefit_mean_y_bin3  = tf1_gauss_y_bin3 -> GetParameter(1);
  double d_prefit_sigma_y_bin3 = tf1_gauss_y_bin3 -> GetParameter(2);

  cout << "d_prefit_p0_y_bin3 = "<< d_prefit_p0_y_bin3 << endl;
  cout << "d_prefit_mean_y_bin3 = "<< d_prefit_mean_y_bin3 << endl;
  cout << "d_prefit_sigma_y_bin3 = "<< d_prefit_sigma_y_bin3 << endl;

  tf1_phi_y_bin3 -> SetParameter(1,d_prefit_mean_y_bin3);
  tf1_phi_y_bin3 -> SetParameter(2,d_prefit_sigma_y_bin3);
  tf1_phi_y_bin3 -> SetParLimits(1,d_prefit_mean_y_bin3-d_prefit_sigma_y_bin3,d_prefit_mean_y_bin3+d_prefit_sigma_y_bin3);
  tf1_phi_y_bin3 -> SetParLimits(2,0.66*d_prefit_sigma_y_bin3,1.5*d_prefit_sigma_y_bin3);

  tf1_phi_y_bin3 -> SetLineColor(kBlue);
  h_sub_y_bin3   -> Fit(tf1_phi_y_bin3,"E+","R",0.98,1.1);
  h_sub_y_bin3   -> SetTitle(Form("%s After BG Subtraction",h_sub_y_bin3->GetTitle() ));
  h_sub_y_bin3   -> SetXTitle("M_{K+K-}");
  h_sub_y_bin3   -> SetYTitle("Counts");

  double d_p0_y_bin3    = tf1_phi_y_bin3    -> GetParameter(0);
  double d_mean_y_bin3  = tf1_phi_y_bin3 -> GetParameter(1);
  double d_sigma_y_bin3 = tf1_phi_y_bin3 -> GetParameter(2);
  double d_p3_y_bin3    = tf1_phi_y_bin3    -> GetParameter(3);


  cout << "d_p0_y_bin3 = "<< d_p0_y_bin3 << endl;
  cout << "d_mean_y_bin3 = "<< d_mean_y_bin3 << endl;
  cout << "d_sigma_y_bin3 = "<< d_sigma_y_bin3 << endl;
  cout << "d_p3_y_bin3 = "<< d_p3_y_bin3 << endl;

  int iBin_3sigint_low_y_bin3 = h_sub_y_bin3 -> FindFixBin(d_mean_y_bin3 - (3*d_sigma_y_bin3));
  int iBin_3sigint_hi_y_bin3  = h_sub_y_bin3 -> FindFixBin(d_mean_y_bin3 + (3*d_sigma_y_bin3));
  double d_3sig_integral_error_y_bin3;

  double d_3sig_integral_y_bin3 = h_sub_y_bin3 -> IntegralAndError(iBin_3sigint_low_y_bin3,iBin_3sigint_hi_y_bin3,d_3sig_integral_error_y_bin3,"");
  double d_dNdy_err_y_bin3 = 1.0/sqrt(d_3sig_integral_y_bin3);
  double d_Y_mT_err_y_bin3 = 1.0/sqrt(d_3sig_integral_y_bin3);

  //  TPaveText * ptext = new TPaveText(0.65,0.5,1.,0.95,"NDCARC");
  TPaveText * ptext_y_bin3 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin3 -> AddText(Form("Mean: %.4f",d_mean_y_bin3));
  ptext_y_bin3 -> AddText(Form("Sigma: %.4f",d_sigma_y_bin3));
  ptext_y_bin3 -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral_y_bin3));
  ptext_y_bin3 -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error_y_bin3));

  TC_sub_y_bin3->cd();
  ptext_y_bin3 -> Draw("same");

  //Fit function
  //fit to Gauss plus constant
  TFormula * tform_single_y_bin4 = new TFormula("tform_single_y_bin4","gaus(0)+[3]");
  TF1 * tf1_phi_y_bin4 = new TF1("polygauss_single_y_bin4",tform_single_y_bin4->GetExpFormula(),1.,1.06);

  //fit to a simple gauss first to seed
  TF1 * tf1_gauss_y_bin4 = new TF1("tf1_gauss_y_bin4","gaus",0.9,1.1);
  TC_sub_y_bin4->cd();
  h_sub_y_bin4 -> Fit(tf1_gauss_y_bin4,"0","R",1.01,1.03);

  double d_prefit_p0_y_bin4    = tf1_gauss_y_bin4 -> GetParameter(0);
  double d_prefit_mean_y_bin4  = tf1_gauss_y_bin4 -> GetParameter(1);
  double d_prefit_sigma_y_bin4 = tf1_gauss_y_bin4 -> GetParameter(2);

  cout << "d_prefit_p0_y_bin4 = "<< d_prefit_p0_y_bin4 << endl;
  cout << "d_prefit_mean_y_bin4 = "<< d_prefit_mean_y_bin4 << endl;
  cout << "d_prefit_sigma_y_bin4 = "<< d_prefit_sigma_y_bin4 << endl;

  tf1_phi_y_bin4 -> SetParameter(1,d_prefit_mean_y_bin4);
  tf1_phi_y_bin4 -> SetParameter(2,d_prefit_sigma_y_bin4);
  tf1_phi_y_bin4 -> SetParLimits(1,d_prefit_mean_y_bin4-d_prefit_sigma_y_bin4,d_prefit_mean_y_bin4+d_prefit_sigma_y_bin4);
  tf1_phi_y_bin4 -> SetParLimits(2,0.66*d_prefit_sigma_y_bin4,1.5*d_prefit_sigma_y_bin4);

  tf1_phi_y_bin4 -> SetLineColor(kBlue);
  h_sub_y_bin4   -> Fit(tf1_phi_y_bin4,"E+","R",0.98,1.1);
  h_sub_y_bin4   -> SetTitle(Form("%s After BG Subtraction",h_sub_y_bin4->GetTitle() ));
  h_sub_y_bin4   -> SetXTitle("M_{K+K-}");
  h_sub_y_bin4   -> SetYTitle("Counts");

  double d_p0_y_bin4    = tf1_phi_y_bin4 -> GetParameter(0);
  double d_mean_y_bin4  = tf1_phi_y_bin4 -> GetParameter(1);
  double d_sigma_y_bin4 = tf1_phi_y_bin4 -> GetParameter(2);
  double d_p3_y_bin4    = tf1_phi_y_bin4 -> GetParameter(3);

  cout << "d_p0_y_bin4 = "<< d_p0_y_bin4 << endl;
  cout << "d_mean_y_bin4 = "<< d_mean_y_bin4 << endl;
  cout << "d_sigma_y_bin4 = "<< d_sigma_y_bin4 << endl;
  cout << "d_p3_y_bin4 = "<< d_p3_y_bin4 << endl;

  int iBin_3sigint_low_y_bin4 = h_sub_y_bin4 -> FindFixBin(d_mean_y_bin4 - (3*d_sigma_y_bin4));
  int iBin_3sigint_hi_y_bin4  = h_sub_y_bin4 -> FindFixBin(d_mean_y_bin4 + (3*d_sigma_y_bin4));
  double d_3sig_integral_error_y_bin4;

  double d_3sig_integral_y_bin4 = h_sub_y_bin4 -> IntegralAndError(iBin_3sigint_low_y_bin4,iBin_3sigint_hi_y_bin4,d_3sig_integral_error_y_bin4,"");
  double d_dNdy_err_y_bin4 = 1.0/sqrt(d_3sig_integral_y_bin4);
  double d_Y_mT_err_y_bin4 = 1.0/sqrt(d_3sig_integral_y_bin4);

  //  TPaveText * ptext = new TPaveText(0.65,0.5,1.,0.95,"NDCARC");
  TPaveText * ptext_y_bin4 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin4 -> AddText(Form("Mean: %.4f",d_mean_y_bin4));
  ptext_y_bin4 -> AddText(Form("Sigma: %.4f",d_sigma_y_bin4));
  ptext_y_bin4 -> AddText(Form("3#sigma Int: %.4f",d_3sig_integral_y_bin4));
  ptext_y_bin4 -> AddText(Form("3#sigma Int Err: %.4f",d_3sig_integral_error_y_bin4));

  TC_sub_y_bin4->cd();
  ptext_y_bin4 -> Draw("same");




  TFile * tf_out = new TFile("Phi_out_single_tot_eff_corr_invM.root","RECREATE");
  h_TPC2_TOF3_invM_y_bin2 -> Write();
  h_mx_TPC2_TOF3_invM_y_bin2 -> Write();
  h_sub_y_bin2->Write();

  h_TPC2_TOF3_invM_y_bin3 -> Write();
  h_mx_TPC2_TOF3_invM_y_bin3 -> Write();
  h_sub_y_bin3->Write();

  h_TPC2_TOF3_invM_y_bin4 -> Write();
  h_mx_TPC2_TOF3_invM_y_bin4 -> Write();
  h_sub_y_bin4->Write();

  h_TPC2_TOF3_invM -> Write();
  h_mx_TPC2_TOF3_invM -> Write();
  h_sub->Write();

  tf_out -> Close();

  double d_int_nevt      = (d_N_evt != 0.0)    ? d_3sig_integral/d_N_evt : 0.0;
  d_dNdy_err        = (d_N_evt != 0.0)    ? d_dNdy_err/d_N_evt : 0.0;
  d_Y_mT_err        = (d_N_evt != 0.0)    ? d_Y_mT_err/d_N_evt : 0.0;

  // from phiDecay.C
  double d_int_effacc    = 0.21289;//0.182376;
  double d_int_eff       = d_int_nevt / d_int_effacc;
  d_dNdy_err             = d_dNdy_err / d_int_effacc;
  d_Y_mT_err             = d_Y_mT_err / d_int_effacc;

  // delta eta
  double d_delta_eta     = fabs(0.0 -1.47);
  double d_delta_y       = 2.0;
  double d_delta_mT      = 2.0;
  double d_int_dmTdeta   = d_int_eff/(d_delta_y*d_delta_mT);
  double d_2pi_dmTdeta   = d_int_dmTdeta/(2.*TMath::Pi());
  d_Y_mT_err             = d_Y_mT_err/(d_delta_y*d_delta_mT);
  d_Y_mT_err             = d_Y_mT_err/(2.*TMath::Pi());

  double d_int_evt_delta_y = d_int_eff/d_delta_y;
  double d_dNdy                   = d_int_evt_delta_y;
  d_dNdy_err               = d_dNdy_err/d_delta_y;

  // TC_inv   -> SaveAs("/home/jbryslaw/MPC_EX/PHI_UPDATE_08_07_2018/h_inv.pdf");
  // TC_sub   -> SaveAs("/home/jbryslaw/MPC_EX/PHI_UPDATE_08_07_2018/h_sub.pdf");

  TCanvas * TC_dNdy= new TCanvas("TC_dNdy","TC_dNdy",480,480);
  TH1D * h_dNdy = new TH1D("h_dNdy","#phi dN/dy BG subtracted and EFF corrected",100,-2,2);

  TFile * tf_qa = new TFile("h_K_PID_QA_tot.root","READ");
  TH1D * h_K_pl_pT = (TH1D*)tf_qa -> Get("h_K_pl_pT");
  TH1D * h_K_min_pT = (TH1D*)tf_qa -> Get("h_K_min_pT");

  int iBin_low = h_K_min_pT -> FindFixBin(0.0);
  int iBin_hi = h_K_min_pT -> FindFixBin(3.0);

  double d_K_pl_int = h_K_pl_pT -> Integral(iBin_low,iBin_hi);
  double d_K_pl_int_err = (d_K_pl_int != 0.0) ?  1.0/d_K_pl_int : 0.0;
  double d_K_min_int = h_K_min_pT -> Integral(iBin_low,iBin_hi);
  double d_K_min_int_err = (d_K_min_int != 0.0) ?  1.0/d_K_min_int : 0.0;
  double d_phi_int_err = (d_3sig_integral != 0.0) ? 1.0/d_3sig_integral : 0.0;

  cout<<" Kpl: "<<d_K_pl_int/(2.*d_N_evt)<<endl;

  // new TCanvas();
  // double d_K_pl_R = d_3sig_integral/d_K_pl_int;
  // double d_K_pl_R_err = sqrt(d_phi_int_err*d_phi_int_err+d_K_pl_int_err*d_K_pl_int_err);

  // double d_K_min_R = d_3sig_integral/d_K_min_int;
  // double d_K_min_R_err = sqrt(d_phi_int_err*d_phi_int_err+d_K_min_int_err*d_K_min_int_err);
  // TH1D * h_K_pl_R = new TH1D("h_K_pl_R","#phi / K+",100,0.0,3.0);
  // TH1D * h_K_min_R = new TH1D("h_K_min_R","#phi / K-",100,0.0,3.0);

  // h_K_pl_R -> Fill(1.5,d_K_pl_R);
  // h_K_min_R -> Fill(1.5,d_K_min_R);

  // int iBin_k_pl = h_K_pl_R -> FindFixBin(1.5);
  // h_K_pl_R -> SetBinError(iBin_k_pl,d_K_pl_R_err);

  //   int iBin_k_min = h_K_min_R -> FindFixBin(1.5);
  // h_K_min_R -> SetBinError(iBin_k_min,d_K_min_R_err);


  // h_K_min_R -> SetMarkerStyle(21);
  // h_K_pl_R -> SetMarkerStyle(21);
  // h_K_min_R -> SetMarkerColor(kRed);
  // //h_K_pl_R -> SetMarkerColor(21);

  // h_K_pl_R -> SetMaximum(d_K_min_R*1.1);
  // h_K_pl_R -> SetMinimum(d_K_pl_R*0.9);
  // h_K_pl_R -> Draw();
  // h_K_min_R -> Draw("sames");
  // gPad -> BuildLegend();

  // cout<< d_K_min_R <<" "<<d_K_pl_R<<endl;

  TH1D * h_phi_TPC2_TOF3_tight_y = (TH1D*) tf_qa -> Get("h_phi_TPC2_TOF3_tight_y");
  double d_ymean = h_phi_TPC2_TOF3_tight_y -> GetMean();
  double d_ycm = -1.57;

  h_dNdy -> Fill(d_ymean - d_ycm,d_dNdy);
  int iBin_ymcm = h_dNdy -> FindFixBin(d_ymean - d_ycm);
  if(d_dNdy != 0.0)  h_dNdy -> SetBinError(iBin_ymcm, d_dNdy_err );
  h_dNdy -> SetMarkerStyle(20);
  h_dNdy -> SetXTitle("y-y_{cm}");
  h_dNdy -> SetYTitle("dN/dy");

  h_dNdy -> SetMinimum(d_dNdy *0.9);
  h_dNdy -> SetMaximum(d_dNdy* 1.1);
  TC_dNdy->cd();
  h_dNdy -> Draw("E");

  //  TC_dNdy -> SaveAs("/home/jbryslaw/MPC_EX/PHI_UPDATE_08_07_2018/TC_dNdy.pdf");

  //mT
  double d_Y_mT = d_2pi_dmTdeta;

  TCanvas * TC_dmT = new TCanvas("TC_dmT","TC_dmT",480,480);
  TH1D * h_mT = new TH1D("h_mT","#phi m_{T} BG subtracted and EFF corrected",100,0.0,2.0);
  h_mT -> Fill(1.0,d_Y_mT);
  int iBin_mT1 = h_mT -> FindFixBin(1.0);
  if(d_Y_mT != 0.0)  h_mT -> SetBinError(iBin_mT1, d_Y_mT_err );
  h_mT -> SetMarkerStyle(20);
  h_mT -> SetXTitle("m_{T}-m{0} (GeV/c^{2})");
  h_mT -> SetYTitle("(1/2#pim_{T})d^{2}N/dydm_{T} (GeV/c^{2})^{-2}");

  h_mT -> SetMinimum(d_Y_mT *0.9);
  h_mT -> SetMaximum(d_Y_mT* 1.1);
  TC_dmT->cd();
  h_mT -> Draw("E");

  //  TC_dmT -> SaveAs("/home/jbryslaw/MPC_EX/PHI_UPDATE_08_07_2018/TC_dmT.pdf");
  // h_out -> SetName("h_out");
  // h_out -> SaveAs("h_out.root");

  return;
}//void out_single()
