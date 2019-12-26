//==============================================================================
// A Macro to build event plane
// author: Ding Chen
// date : Aug 6, 2019
//
//==============================================================================

// C++ headers
#include <iostream>
#include <cstdio>
#include <vector>
#include <set>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"

// PicoDst headers
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoHelix.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoTrackCovMatrix.h"



//==================================== Main Funcition ==========================
void EventPlaneBuilder(const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_EventPlane")
{
  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);

  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !picoReader->chain() ) {
      std::cout << "No chain has been found." << std::endl;
  }

  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;

  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Number of events to read: " << events2read << std::endl;

  outFile.Append(".picoDst.result.root");

  TFile *outputFile = new TFile(outFile,"recreate");

  // Read information from first run
  TFile *f1 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/files/Event_Plane/EventPlane_4.5_Run1_pTcut.root","read");
  if(f1->IsOpen()) cout<<"EventPlane_4.5_Run1 loaded successfully!"<<endl;
  else {cout<<"No EventPlane_4.5_Run1, quit!"<<endl; return;}

  TH1D * h_1 = (TH1D*)f1->Get("h_Q1_y_sub1");
  TH1D * h_2 = (TH1D*)f1->Get("h_Q1_x_sub1");
  TH1D * h_3 = (TH1D*)f1->Get("h_Q1_y_sub2");
  TH1D * h_4 = (TH1D*)f1->Get("h_Q1_x_sub2");

  double d_Q1_y_sub1_m = h_1->GetMean();
  double d_Q1_x_sub1_m = h_2->GetMean();

  double d_Q1_y_sub2_m = h_3->GetMean();
  double d_Q1_x_sub2_m = h_4->GetMean();

  // End Read information from first run

  // Read information from second run
  TFile *f2 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/files/Event_Plane/EventPlane_4.5_Run2_pTcut.root","read");
  if(f2->IsOpen()) cout<<"EventPlane_4.5_Run2 loaded successfully!"<<endl;
  else {cout<<"No EventPlane_4.5_Run2, quit!"<<endl; return;}

  TH1D * h_isin_rec1s[20];
  TH1D * h_icos_rec1s[20];

  TH1D * h_isin_rec2s[20];
  TH1D * h_icos_rec2s[20];

  double d_isin_rec1s[20];
  double d_icos_rec1s[20];

  double d_isin_rec2s[20];
  double d_icos_rec2s[20];

  for(int i=0;i<20;i++)
  {
    h_isin_rec1s[i] = (TH1D*)f2->Get(Form("h_sin_iPsi1_rec1[%d]" ,i));
    h_icos_rec1s[i] = (TH1D*)f2->Get(Form("h_cos_iPsi1_rec1[%d]" ,i));

    h_isin_rec2s[i] = (TH1D*)f2->Get(Form("h_sin_iPsi1_rec2[%d]" ,i));
    h_icos_rec2s[i] = (TH1D*)f2->Get(Form("h_cos_iPsi1_rec2[%d]" ,i));

    d_isin_rec1s[i] = h_isin_rec1s[i]->GetMean();
    d_icos_rec1s[i] = h_icos_rec1s[i]->GetMean();

    d_isin_rec2s[i] = h_isin_rec2s[i]->GetMean();
    d_icos_rec2s[i] = h_icos_rec2s[i]->GetMean();
    // std::cout << Form("d_isin_rec1s[%d]" ,i) << d_isin_rec1s[i] << std::endl;
    // std::cout << Form("d_icos_rec1s[%d]" ,i) << d_icos_rec1s[i] << std::endl;
    // std::cout << Form("d_isin_rec2s[%d]" ,i) << d_isin_rec2s[i] << std::endl;
    // std::cout << Form("d_icos_rec2s[%d]" ,i) << d_icos_rec2s[i] << std::endl;

  }

  // END Read information from second run

  // Read Event Plane Resolution
  TFile *f3 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/files/Event_Plane/EventPlaneResolution_4.5_Run1_pTcut.root","read");
  if(f3->IsOpen()) cout<<"EventPlane_4.5_Run1 loaded successfully!"<<endl;
  else {cout<<"No EventPlane_4.5_Run2, quit!"<<endl; return;}

  TH1D * h_EPR_1 = (TH1D*)f3->Get("h_cos_Psi1_flat_diff");

  double d_cos_Psi1_flat_diff_avg = h_EPR_1->GetMean();
  double d_resolution1_sub        = TMath::Sqrt(d_cos_Psi1_flat_diff_avg);
  double d_resolution1_full_apx   = d_resolution1_sub*TMath::Sqrt(2);

  // END Read Event Plane Resolution

  double  d_zvtx  = -9999.0;
  std::vector <unsigned int> triggerIDs;
  // Event cut set up

  // Histogramning

  // Test Plots
  TH1D *  h_Q1_y_diff1   = new TH1D("h_Q1_y_diff1","h_Q1_y_diff1",100001,-500.005,500.005);
  TH1D *  h_Q1_x_diff1   = new TH1D("h_Q1_x_diff1","h_Q1_x_diff1",100001,-500.005,500.005);

  TH1D *  h_Q1_y_diff2   = new TH1D("h_Q1_y_diff2","h_Q1_y_diff2",100001,-500.005,500.005);
  TH1D *  h_Q1_x_diff2   = new TH1D("h_Q1_x_diff2","h_Q1_x_diff2",100001,-500.005,500.005);

  // End Test Plots
  TH1D *  h_eta         = new TH1D("h_eta","h_eta",300,-2.5,0.5);
  TH1D *  h_eta0        = new TH1D("h_eta0","h_eta0",300,-2.5,0.5);
  TH1D *  h_eta1        = new TH1D("h_eta1","h_eta1",300,-2.5,0.5);
  TH1D *  h_Psi1_sub1   = new TH1D("h_Psi1_sub1","h_Psi1_sub1",633,-3.165,3.165);
  TH1D *  h_Psi1_sub2   = new TH1D("h_Psi1_sub2","h_Psi1_sub2",633,-3.165,3.165);

  TH1D *  h_Psi1_rec1   = new TH1D("h_Psi1_rec1","h_Psi1_rec1",633,-3.165,3.165);
  TH1D *  h_Psi1_rec2   = new TH1D("h_Psi1_rec2","h_Psi1_rec2",633,-3.165,3.165);

  TH1D *  h_Psi1_flat1   = new TH1D("h_Psi1_flat1","h_Psi1_flat1",633,-3.165,3.165);
  TH1D *  h_Psi1_flat2   = new TH1D("h_Psi1_flat2","h_Psi1_flat2",633,-3.165,3.165);

  // TH1D *  h_cos_Psi1_flat_diff = new TH1D("h_cos_Psi1_flat_diff","h_cos_Psi1_flat_diff",1001,-5.005,5.005);

  TH1D *  h_v1_obs      = new TH1D("h_v1_obs","h_v1_obs",1001,-5.005,5.005);
  TH1D *  h_v1_obs_sub1 = new TH1D("h_v1_obs_sub1","h_v1_obs_sub1",1001,-5.005,5.005);
  TH1D *  h_v1_obs_sub2 = new TH1D("h_v1_obs_sub2","h_v1_obs_sub2",1001,-5.005,5.005);

  TH1D *  h_v1          = new TH1D("h_v1","h_v1",1001,-5.005,5.005);
  TH1D *  h_v1_sub1     = new TH1D("h_v1_sub1","h_v1_sub1",1001,-5.005,5.005);
  TH1D *  h_v1_sub2     = new TH1D("h_v1_sub2","h_v1_sub2",1001,-5.005,5.005);
  // Directed flow

  TH2D *  h2_PRO_plus_v1_vs_y_sub1 = new TH2D("h2_PRO_plus_v1_vs_y_sub1","h2_PRO_plus_v1_vs_y_sub1",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_plus_v1_vs_y_sub2 = new TH2D("h2_PRO_plus_v1_vs_y_sub2","h2_PRO_plus_v1_vs_y_sub2",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_plus_v1_vs_y      = new TH2D("h2_PRO_plus_v1_vs_y","h2_PRO_plus_v1_vs_y",1001,-5.005,5.005,1001,-5.005,5.005);

  TH2D *  h2_PRO_plus_v1_vs_pT_sub1 = new TH2D("h2_PRO_plus_v1_vs_pT_sub1","h2_PRO_plus_v1_vs_pT_sub1",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_plus_v1_vs_pT_sub2 = new TH2D("h2_PRO_plus_v1_vs_pT_sub2","h2_PRO_plus_v1_vs_pT_sub2",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_plus_v1_vs_pT      = new TH2D("h2_PRO_plus_v1_vs_pT","h2_PRO_plus_v1_vs_pT",1001,-5.005,5.005,1001,-5.005,5.005);

  TH2D *  h2_PRO_minus_v1_vs_y_sub1 = new TH2D("h2_PRO_minus_v1_vs_y_sub1","h2_PRO_minus_v1_vs_y_sub1",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_y_sub2 = new TH2D("h2_PRO_minus_v1_vs_y_sub2","h2_PRO_minus_v1_vs_y_sub2",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_y      = new TH2D("h2_PRO_minus_v1_vs_y","h2_PRO_minus_v1_vs_y",1001,-5.005,5.005,1001,-5.005,5.005);

  TH2D *  h2_PRO_minus_v1_vs_pT_sub1 = new TH2D("h2_PRO_minus_v1_vs_pT_sub1","h2_PRO_minus_v1_vs_pT_sub1",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_pT_sub2 = new TH2D("h2_PRO_minus_v1_vs_pT_sub2","h2_PRO_minus_v1_vs_pT_sub2",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_pT      = new TH2D("h2_PRO_minus_v1_vs_pT","h2_PRO_minus_v1_vs_pT",1001,-5.005,5.005,1001,-5.005,5.005);

  TProfile * tp1 = new TProfile("tp1","tp1",15,-3,0,-5.0,5.0);
  TProfile * tp1_sub1 = new TProfile("tp1_sub1","tp1_sub1",15,-3,0,-5.0,5.0);
  TProfile * tp1_sub2 = new TProfile("tp1_sub2","tp1_sub2",15,-3,0,-5.0,5.0);

  TProfile * tp1_PRO_v1_vs_y = new TProfile("tp1_PRO_v1_vs_y","tp1_PRO_v1_vs_y",15,-3,0,-5.0,5.0);
  TProfile * tp1_PRO_v1_vs_y_sub1 = new TProfile("tp1_PRO_v1_vs_y_sub1","tp1_PRO_v1_vs_y_sub1",15,-3,0,-5.0,5.0);
  TProfile * tp1_PRO_v1_vs_y_sub2 = new TProfile("tp1_PRO_v1_vs_y_sub2","tp1_PRO_v1_vs_y_sub2",15,-3,0,-5.0,5.0);

  TProfile * tp1_PRO_v1_vs_pT = new TProfile("tp1_PRO_v1_vs_pT","tp1_PRO_v1_vs_pT",50,0,5,-5.0,5.0);
  TProfile * tp1_PRO_v1_vs_pT_sub1 = new TProfile("tp1_PRO_v1_vs_pT_sub1","tp1_PRO_v1_vs_pT_sub1",50,0,5,-5.0,5.0);
  TProfile * tp1_PRO_v1_vs_pT_sub2 = new TProfile("tp1_PRO_v1_vs_pT_sub2","tp1_PRO_v1_vs_pT_sub2",50,0,5,-5.0,5.0);

  TProfile * tp1_PI_plus_v1_vs_y = new TProfile("tp1_PI_plus_v1_vs_y","tp1_PI_plus_v1_vs_y",15,-3,0,-5.0,5.0);
  TProfile * tp1_PI_plus_v1_vs_y_sub1 = new TProfile("tp1_PI_plus_v1_vs_y_sub1","tp1_PI_plus_v1_vs_y_sub1",15,-3,0,-5.0,5.0);
  TProfile * tp1_PI_plus_v1_vs_y_sub2 = new TProfile("tp1_PI_plus_v1_vs_y_sub2","tp1_PI_plus_v1_vs_y_sub2",15,-3,0,-5.0,5.0);

  TProfile * tp1_PI_plus_v1_vs_pT = new TProfile("tp1_PI_plus_v1_vs_pT","tp1_PI_plus_v1_vs_pT",50,0,5,-5.0,5.0);
  TProfile * tp1_PI_plus_v1_vs_pT_sub1 = new TProfile("tp1_PI_plus_v1_vs_pT_sub1","tp1_PI_plus_v1_vs_pT_sub1",50,0,5,-5.0,5.0);
  TProfile * tp1_PI_plus_v1_vs_pT_sub2 = new TProfile("tp1_PI_plus_v1_vs_pT_sub2","tp1_PI_plus_v1_vs_pT_sub2",50,0,5,-5.0,5.0);

  TProfile * tp1_PI_minus_v1_vs_y = new TProfile("tp1_PI_minus_v1_vs_y","tp1_PI_minus_v1_vs_y",15,-3,0,-5.0,5.0);
  TProfile * tp1_PI_minus_v1_vs_y_sub1 = new TProfile("tp1_PI_minus_v1_vs_y_sub1","tp1_PI_minus_v1_vs_y_sub1",15,-3,0,-5.0,5.0);
  TProfile * tp1_PI_minus_v1_vs_y_sub2 = new TProfile("tp1_PI_minus_v1_vs_y_sub2","tp1_PI_minus_v1_vs_y_sub2",15,-3,0,-5.0,5.0);

  TProfile * tp1_PI_minus_v1_vs_pT = new TProfile("tp1_PI_minus_v1_vs_pT","tp1_PI_minus_v1_vs_pT",50,0,5,-5.0,5.0);
  TProfile * tp1_PI_minus_v1_vs_pT_sub1 = new TProfile("tp1_PI_minus_v1_vs_pT_sub1","tp1_PI_minus_v1_vs_pT_sub1",50,0,5,-5.0,5.0);
  TProfile * tp1_PI_minus_v1_vs_pT_sub2 = new TProfile("tp1_PI_minus_v1_vs_pT_sub2","tp1_PI_minus_v1_vs_pT_sub2",50,0,5,-5.0,5.0);

  // PID Plots
  TH1D *  h_PRO_plus_mT_Diff = new TH1D("h_PRO_plus_mT_Diff","h_PRO_plus_mT_Diff",1000,0.0,5.0);
  TH1D *  h_PRO_plus_pT = new TH1D("h_PRO_plus_pT","Proton+ pT",1000,0.0,5.0);
  TH1D *  h_PRO_plus_y = new TH1D("h_PRO_plus_y","Proton+ y",1000,-5.0,0.0);
  TH2D *  h2_PRO_plus_pT_vs_y = new TH2D("h2_PRO_plus_pT_vs_y","h2_PRO_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);

  TH1D *  h_PI_plus_mT_Diff  = new TH1D("h_PI_plus_mT_Diff","h_PI_plus_mT_Diff",1000,0.0,5.0);
  TH1D *  h_PI_plus_pT  = new TH1D("h_PI_plus_pT","Pion+ pT",1000,0.0,5.0);
  TH1D *  h_PI_plus_y  = new TH1D("h_PI_plus_y","Pion+ y",1000,-5.0,0.0);
  TH2D *  h2_PI_plus_pT_vs_y  = new TH2D("h2_PI_plus_pT_vs_y","h2_PI_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);

  TH1D *  h_PI_minus_mT_Diff  = new TH1D("h_PI_minus_mT_Diff","h_PI_minus_mT_Diff",1000,0.0,5.0);
  TH1D *  h_PI_minus_pT  = new TH1D("h_PI_minus_pT","Pion+ pT",1000,0.0,5.0);
  TH1D *  h_PI_minus_y  = new TH1D("h_PI_minus_y","Pion+ y",1000,-5.0,0.0);
  TH2D *  h2_PI_minus_pT_vs_y  = new TH2D("h2_PI_minus_pT_vs_y","h2_PI_minus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);

  TH2D *  h2_PRO_plus_eta_vs_y  = new TH2D("h2_PRO_plus_eta_vs_y","h2_PRO_plus_eta_vs_y",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PRO_plus_eta_vs_y_sub1  = new TH2D("h2_PRO_plus_eta_vs_y_sub1","h2_PRO_plus_eta_vs_y_sub1",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PRO_plus_eta_vs_y_sub2  = new TH2D("h2_PRO_plus_eta_vs_y_sub2","h2_PRO_plus_eta_vs_y_sub2",1000,-5.0,0.0,1000,-5.0,0.0);

  TH2D *  h2_PI_plus_eta_vs_y  = new TH2D("h2_PI_plus_eta_vs_y","h2_PI_plus_eta_vs_y",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PI_plus_eta_vs_y_sub1  = new TH2D("h2_PI_plus_eta_vs_y_sub1","h2_PI_plus_eta_vs_y_sub1",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PI_plus_eta_vs_y_sub2  = new TH2D("h2_PI_plus_eta_vs_y_sub2","h2_PI_plus_eta_vs_y_sub2",1000,-5.0,0.0,1000,-5.0,0.0);

  TH2D *  h2_PI_minus_eta_vs_y  = new TH2D("h2_PI_minus_eta_vs_y","h2_PI_minus_eta_vs_y",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PI_minus_eta_vs_y_sub1  = new TH2D("h2_PI_minus_eta_vs_y_sub1","h2_PI_minus_eta_vs_y_sub1",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PI_minus_eta_vs_y_sub2  = new TH2D("h2_PI_minus_eta_vs_y_sub2","h2_PI_minus_eta_vs_y_sub2",1000,-5.0,0.0,1000,-5.0,0.0);

  // Plots for Recentering Run
  // TH1D *  h_Q1_y_sub1   = new TH1D("h_Q1_y_sub1","h_Q1_y_sub1",100001,-500.005,500.005);
  // TH1D *  h_Q1_x_sub1   = new TH1D("h_Q1_x_sub1","h_Q1_x_sub1",100001,-500.005,500.005);
  //
  // TH1D *  h_Q1_y_sub2   = new TH1D("h_Q1_y_sub2","h_Q1_y_sub2",100001,-500.005,500.005);
  // TH1D *  h_Q1_x_sub2   = new TH1D("h_Q1_x_sub2","h_Q1_x_sub2",100001,-500.005,500.005);

  // Plots for Flattening Run

  // TH1D *  h_sin_iPsi1_rec1[20];//   = new TH1D("h_sin_iPsi1_rec1","h_sin_iPsi1_rec1",1001,-5.005,5.005);
  // TH1D *  h_cos_iPsi1_rec1[20];//   = new TH1D("h_cos_iPsi1_rec1","h_cos_iPsi1_rec1",1001,-5.005,5.005);
  //
  // TH1D *  h_sin_iPsi1_rec2[20];//   = new TH1D("h_sin_iPsi1_rec2","h_sin_iPsi1_rec2",1001,-5.005,5.005);
  // TH1D *  h_cos_iPsi1_rec2[20];//   = new TH1D("h_cos_iPsi1_rec2","h_cos_iPsi1_rec2",1001,-5.005,5.005);

  // for(int i=0;i<20;i++)
  // {
  //   h_sin_iPsi1_rec1[i] = new TH1D(Form("h_sin_iPsi1_rec1[%d]" ,i),Form("h_sin_iPsi1_rec1[%d]" ,i),1001,-5.005,5.005);
  //   h_cos_iPsi1_rec1[i] = new TH1D(Form("h_cos_iPsi1_rec1[%d]" ,i),Form("h_cos_iPsi1_rec1[%d]" ,i),1001,-5.005,5.005);
  //
  //   h_sin_iPsi1_rec2[i] = new TH1D(Form("h_sin_iPsi1_rec2[%d]" ,i),Form("h_sin_iPsi1_rec2[%d]" ,i),1001,-5.005,5.005);
  //   h_cos_iPsi1_rec2[i] = new TH1D(Form("h_cos_iPsi1_rec2[%d]" ,i),Form("h_cos_iPsi1_rec2[%d]" ,i),1001,-5.005,5.005);
  //
  // }


  // END Histogramning

  //============================ Event Loop ====================================
  for(Long64_t iEvent = 0; iEvent < events2read; iEvent++)
  {
    std::cout << "Working on event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);

    if( !readEvent )
    {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    StPicoDst *dst = picoReader->picoDst();
    // Retrieve picoDst

    StPicoEvent *event = dst->event();
    // Retrieve event information

    Int_t runId_Dst   = event->runId();
    Int_t eventId_Dst = event->eventId();
    // std::cout << " runId from picoDst = "   <<  runId_Dst << std::endl;
    // std::cout << " eventId from picoDst = " <<  eventId_Dst << std::endl;

    if( !event )
    {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    //============================ Trigger Selection ===========================
    triggerIDs.clear();
    triggerIDs       = event->triggerIds();
    bool b_bad_trig = true;

    // loop for the trigger ids and see if any == 1
    for(unsigned int i=0; i < triggerIDs.size(); i++)
      {
        if(triggerIDs[i] == 1) b_bad_trig = false;
      }

    //=========================== End Trigger Slection =========================

    TVector3 pVtx     = event->primaryVertex();
    // Primary Vertex

    //=========================== Z-VTX Selection ==============================
    TVector3 v3D_vtx  = event->primaryVertex();
    d_zvtx = pVtx.z();
    bool b_bad_zvtx   = ((d_zvtx < 210.0) || (d_zvtx > 212.0));
    //======================== END Z-VTX Selection =============================

    bool b_bad_evt  = b_bad_zvtx || b_bad_trig ;
    if(b_bad_evt) continue;
    // Bad Event Cut

    // Event Plane Set Up

    double d_cos_Psi1_flat_diff  = -999;
    double d_v1_obs_sub1         = -999;
    double d_v1_obs_sub2         = -999;
    double d_v1_obs_sub1_avg     = -999;
    double d_v1_obs_sub2_avg     = -999;
    double d_v1_sub1             = -999;
    double d_v1_sub2             = -999;

    // Event plane resolution for sub-events

    double d_Psi1_flat1 = 0;
    double d_Psi1_flat2 = 0;

    double d_Psi1_rec1_shift = 0;
    double d_Psi1_rec2_shift = 0;

    double d_Psi1_rec1 = 0;
    double d_Psi1_rec2 = 0;

    double d_Psi1_sub1 = 0;
    double d_Psi1_sub2 = 0;

    double d_Q1_x_sub1 = 0;
    double d_Q1_y_sub1 = 0;

    double d_Q1_x_sub2 = 0;
    double d_Q1_y_sub2 = 0;

    double d_Q1_x_rec1 = 0;
    double d_Q1_y_rec1 = 0;

    double d_Q1_x_rec2 = 0;
    double d_Q1_y_rec2 = 0;




    // double d_sin_iPsi1_rec1[20];
    // double d_cos_iPsi1_rec1[20];
    // double d_sin_iPsi1_rec2[20];
    // double d_cos_iPsi1_rec2[20];



    // End Evnet Plane Set Up

    //======================== Track Analysis ==================================
    Int_t nTracks = dst->numberOfTracks();
    //Number of tracks in this event
    double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut // Mon Jul  3 09:16:46 EDT 2017
    if(nTracks <= 10) continue;
    //====================== EP Track loop =====================================
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);
      // Retrieve i-th Pico Track

      if(!picoTrack) continue;
      // PicoTrack Cut

      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      bool b_bad_track    = b_bad_dEdx || b_bad_tracking;
      unsigned short nHits = picoTrack->nHits();
      if(nHits <= 10) continue;

      if(b_bad_track) continue;
      // Track-Level Cut

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      StPicoBTofPidTraits *trait        = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      double d_tofBeta                  = -999;
      if(trait) d_tofBeta               = trait->btofBeta();

      double d_eta       = picoTrack->pMom().Eta();
      double d_Phi       = picoTrack->pMom().Phi();
      double d_pT        = picoTrack->pPt();
      if(d_pT > 2.1) continue;

      double d_Q1_x_itrk = d_pT * TMath::Cos(d_Phi);
      double d_Q1_y_itrk = d_pT * TMath::Sin(d_Phi);






      if(d_eta<-1.26)
      {
        h_eta ->Fill(d_eta);
        h_eta0->Fill(d_eta);

        d_Q1_x_sub1 += d_Q1_x_itrk;
        d_Q1_y_sub1 += d_Q1_y_itrk;
      }
      else if(d_eta>-1.08)
      {
        h_eta ->Fill(d_eta);
        h_eta1->Fill(d_eta);

        d_Q1_x_sub2 += d_Q1_x_itrk;
        d_Q1_y_sub2 += d_Q1_y_itrk;
      }

      if(d_tofBeta == -999) continue;
      // TOF Beta Cut
    }

    //================== END EP Track loop =====================================

    if(d_Q1_x_sub1 != 0 || d_Q1_y_sub1 != 0)
    {
      d_Psi1_sub1 = TMath::ATan2(d_Q1_y_sub1,d_Q1_x_sub1);
      d_Psi1_rec1 = TMath::ATan2((d_Q1_y_sub1-d_Q1_y_sub1_m),(d_Q1_x_sub1-d_Q1_x_sub1_m));

      d_Psi1_rec1_shift = 0;
      for(int i=0; i<20;i++)
      {
        d_Psi1_rec1_shift += (2*(- (d_isin_rec1s[i])*TMath::Cos((i+1)*d_Psi1_rec1)
                                 + (d_icos_rec1s[i])*TMath::Sin((i+1)*d_Psi1_rec1))/(i+1));
      }
      d_Psi1_flat1 = d_Psi1_rec1 + d_Psi1_rec1_shift;

      h_Psi1_flat1->Fill(d_Psi1_flat1);
      h_Psi1_rec1->Fill(d_Psi1_rec1);
      h_Psi1_sub1->Fill(d_Psi1_sub1);

      h_Q1_y_diff1->Fill(d_Q1_y_sub1-d_Q1_y_sub1_m);
      h_Q1_x_diff1->Fill(d_Q1_x_sub1-d_Q1_x_sub1_m);

      // h_Q1_y_sub1->Fill(d_Q1_y_sub1);
      // h_Q1_x_sub1->Fill(d_Q1_x_sub1);

      // 20 orders of shifts for flatten run


      // for(int i=0; i<20;i++)
      // {
      //   d_sin_iPsi1_rec1[i] = TMath::Sin((i+1)*d_Psi1_rec1);
      //   d_cos_iPsi1_rec1[i] = TMath::Cos((i+1)*d_Psi1_rec1);
      //
      //   h_sin_iPsi1_rec1[i] -> Fill(d_sin_iPsi1_rec1[i]);
      //   h_cos_iPsi1_rec1[i] -> Fill(d_cos_iPsi1_rec1[i]);
      //
      // }

    }
    if(d_Q1_x_sub2 != 0 || d_Q1_y_sub2 != 0)
    {
      d_Psi1_sub2 = TMath::ATan2(d_Q1_y_sub2,d_Q1_x_sub2);
      d_Psi1_rec2 = TMath::ATan2((d_Q1_y_sub2-d_Q1_y_sub2_m),(d_Q1_x_sub2-d_Q1_x_sub2_m));

      d_Psi1_rec2_shift = 0;
      for(int i=0; i<20;i++)
      {
        d_Psi1_rec2_shift += (2*(- (d_isin_rec2s[i])*TMath::Cos((i+1)*d_Psi1_rec2)
                                 + (d_icos_rec2s[i])*TMath::Sin((i+1)*d_Psi1_rec2))/(i+1));
      }
      d_Psi1_flat2 = d_Psi1_rec2 + d_Psi1_rec2_shift;


      h_Psi1_flat2->Fill(d_Psi1_flat2);
      h_Psi1_rec2->Fill(d_Psi1_rec2);
      h_Psi1_sub2->Fill(d_Psi1_sub2);

      h_Q1_y_diff2->Fill(d_Q1_y_sub2-d_Q1_y_sub2_m);
      h_Q1_x_diff2->Fill(d_Q1_x_sub2-d_Q1_x_sub2_m);
      // h_Q1_y_sub2->Fill(d_Q1_y_sub2);
      // h_Q1_x_sub2->Fill(d_Q1_x_sub2);

      // 20 orders of shifts for flatten run

      // for(int i=0; i<20;i++)
      // {
      //   d_sin_iPsi1_rec2[i] = TMath::Sin((i+1)*d_Psi1_rec2);
      //   d_cos_iPsi1_rec2[i] = TMath::Cos((i+1)*d_Psi1_rec2);
      //
      //
      //   h_sin_iPsi1_rec2[i] -> Fill(d_sin_iPsi1_rec2[i]);
      //   h_cos_iPsi1_rec2[i] -> Fill(d_cos_iPsi1_rec2[i]);
      //
      // }
    }

    d_cos_Psi1_flat_diff = TMath::Cos(d_Psi1_flat1-d_Psi1_flat2);
    // if(d_cos_Psi1_flat_diff != -999) h_cos_Psi1_flat_diff->Fill(d_cos_Psi1_flat_diff);
    // Subevents Event Plane Resolution

    //=================== Observed Flow Track Loop =============================
    // FLow Track Set Up
    d_v1_obs_sub1    = 0;
    d_v1_obs_sub2    = 0;
    Int_t  nSub1trks = 0;
    Int_t  nSub2trks = 0;

    double d_PI_m    = 0.13957018;
    double d_PRO_m   = 0.9382720813;
    double d_K_m     = 0.493677;
    double d_E_m     = 0.0005109989;
    // END Flow Track Set Up

    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);
      // Retrieve i-th Pico Track

      if(!picoTrack) continue;
      // PicoTrack Cut

      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      bool b_bad_track    = b_bad_dEdx || b_bad_tracking;
      if(b_bad_track) continue;
      // Track-Level Cut

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      StPicoBTofPidTraits *trait        = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      double d_tofBeta                  = -999;
      if(trait) d_tofBeta               = trait->btofBeta();
      double d_tofBeta0                 = -999;
      if(trait) d_tofBeta0              = trait->btofBeta();

      double d_eta         = picoTrack->pMom().Eta();
      double d_Phi         = picoTrack->pMom().Phi();
      double d_pT          = picoTrack->pPt();


      // if(d_eta<-1.26)
      // {
      //   double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat1);
      //
      //   d_v1_obs_sub1 += d_v1_obs_itrk1;
      //   nSub1trks     ++;
      //
      // }
      // else if(d_eta>-1.08)
      // {
      //   double d_v1_obs_itrk2 = TMath::Cos(d_Phi - d_Psi1_flat2);
      //
      //   d_v1_obs_sub2 += d_v1_obs_itrk2;
      //   nSub2trks     ++;
      // }


      // d_v1_obs_sub1_avg = d_v1_obs_sub1/nSub1trks;
      // d_v1_obs_sub2_avg = d_v1_obs_sub2/nSub2trks;
      // Has some Problem here, should be put outside the track loop!!!
      // d_v1_sub1 = d_v1_obs_sub1_avg/d_resolution1_full_apx;
      // d_v1_sub2 = d_v1_obs_sub2_avg/d_resolution1_full_apx;


      // if(d_v1_obs_sub1_avg != -999)
      // {
      //   h_v1_obs_sub1->Fill(d_v1_obs_sub1_avg);
      //   h_v1_obs     ->Fill(d_v1_obs_sub1_avg);
      //
      //   h_v1_sub1    ->Fill(d_v1_sub1);
      //   h_v1         ->Fill(d_v1_sub1);
      //
      // }
      // if(d_v1_obs_sub2_avg != -999)
      // {
      //   h_v1_obs_sub2->Fill(d_v1_obs_sub2_avg);
      //   h_v1_obs     ->Fill(d_v1_obs_sub2_avg);
      //
      //   h_v1_sub2    ->Fill(d_v1_sub2);
      //   h_v1         ->Fill(d_v1_sub2);
      // }
      //=========================== PID ========================================
      if(d_tofBeta == -999) continue;
      // TOF Beta Cut
      double d_TPCnSigmaElectron = fabs(picoTrack->nSigmaElectron());
      double d_TPCnSigmaPion     = fabs(picoTrack->nSigmaPion());
      double d_TPCnSigmaProton   = fabs(picoTrack->nSigmaProton());
      double d_TPCnSigmaKaon     = fabs(picoTrack->nSigmaKaon());

      bool b_PI  = fabs(d_TPCnSigmaPion) < d_SigmaCutLevel;
      bool b_PRO = fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
      bool b_K   = fabs(d_TPCnSigmaKaon) < d_SigmaCutLevel;
      bool b_E   = fabs(d_TPCnSigmaElectron)< d_SigmaCutLevel;

      if( b_E
         && (d_TPCnSigmaElectron < d_TPCnSigmaPion)
         && (d_TPCnSigmaElectron < d_TPCnSigmaProton)
         && (d_TPCnSigmaElectron < d_TPCnSigmaKaon) )
        { b_E = true; b_PI = false; b_PRO = false; b_K = false;}

      if( b_PI
         && (d_TPCnSigmaPion < d_TPCnSigmaElectron)
         && (d_TPCnSigmaPion < d_TPCnSigmaKaon)
         && (d_TPCnSigmaPion < d_TPCnSigmaProton) )
         { b_E = false; b_PI = true; b_PRO = false; b_K = false;}

      if( b_PRO
         && (d_TPCnSigmaProton < d_TPCnSigmaElectron)
         && (d_TPCnSigmaProton < d_TPCnSigmaPion)
         && (d_TPCnSigmaProton < d_TPCnSigmaKaon) )
         { b_E = false; b_PI = false; b_PRO = true; b_K = false;}

      if( b_K
         && (d_TPCnSigmaKaon < d_TPCnSigmaElectron)
         && (d_TPCnSigmaKaon < d_TPCnSigmaProton)
         && (d_TPCnSigmaKaon < d_TPCnSigmaPion) )
         { b_E = false; b_PI = false; b_PRO = false; b_K = true;}

      double d_charge  = picoTrack->charge();

      double d_px0        = picoTrack->pMom().x();
      double d_py0        = picoTrack->pMom().y();
      double d_pz0        = picoTrack->pMom().z();
      double d_pT0        = picoTrack->pPt();
      double d_mom0       = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

      double mass2      = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

      double d_E_PRO      = sqrt(d_PRO_m*d_PRO_m + d_mom0*d_mom0);
      double d_E_PI       = sqrt(d_PI_m*d_PI_m + d_mom0*d_mom0);
      double d_E_K        = sqrt(d_K_m*d_K_m + d_mom0*d_mom0);

      double d_mT_PRO      = sqrt(d_pT0*d_pT0 + d_PRO_m*d_PRO_m);
      double d_mT_PI       = sqrt(d_pT0*d_pT0 + d_PI_m*d_PI_m);
      double d_mT_K        = sqrt(d_pT0*d_pT0 + d_K_m*d_K_m);
      double d_mT_E        = sqrt(d_pT0*d_pT0 + d_E_m*d_E_m);


      double d_y_PRO      = 0.5*TMath::Log((d_E_PRO + d_pz0)/(d_E_PRO - d_pz0));
      double d_y_PI       = 0.5*TMath::Log((d_E_PI + d_pz0)/(d_E_PI - d_pz0));
      double d_y_K        = 0.5*TMath::Log((d_E_K + d_pz0)/(d_E_K - d_pz0));

      if(d_charge > 0.0 && d_eta<=-1.17)
        {
          // if( b_E
          //   && (fabs(d_TPCnSigmaElectron) < 3.0)
          //   && (mass2 < 0.005)
          //   && (mass2 > -0.008)
          //   && (d_pT0 < 0.3)){
          //
          // }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35)){
          //
          // }
          if(
            // b_PI
            // && (fabs(d_TPCnSigmaPion) < 2.0)
            // && (mass2 > -0.03)
            // && (mass2 < 0.06)
            // && ((d_pT0 > 0.4) || (mass2 > 0.008))

            (trait==NULL && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0))
                ||(trait && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 && (mass2 > -0.15) && (mass2 < 0.12))

            // (d_pT0<1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) )
            // ||(d_pT0>=1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 &&  ((d_pT0 > 0.4) || (mass2 > 0.008)) && (mass2 > -0.03) && (mass2 < 0.06))
          ){
            double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

            d_v1_sub1 = d_v1_obs_itrk1/d_resolution1_full_apx;
            tp1_PI_plus_v1_vs_y_sub1->Fill(d_y_PI,d_v1_sub1);
            tp1_PI_plus_v1_vs_y->Fill(d_y_PI,d_v1_sub1);

            tp1_PI_plus_v1_vs_pT_sub1->Fill(d_pT0,d_v1_sub1);
            tp1_PI_plus_v1_vs_pT->Fill(d_pT0,d_v1_sub1);

            h_PI_plus_pT -> Fill(d_pT0);
            h_PI_plus_y  -> Fill(d_y_PI);
            h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

            h2_PI_plus_eta_vs_y->Fill(d_y_PI,d_eta);
            h2_PI_plus_eta_vs_y_sub1->Fill(d_y_PI,d_eta);
          }
          if(
            // b_PRO
            // && (fabs(d_TPCnSigmaProton) < 2.0)
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)

            (trait==NULL && d_pT0>0.4 && d_pT0<2.0 && (fabs(d_TPCnSigmaProton) < 2.0))
                ||(trait && d_pT0>0.4 && d_pT0<2.0 && (fabs(d_TPCnSigmaProton) < 2.0) && d_tofBeta0>0 && (mass2 > 0.4) && (mass2 < 1.4))

            // (d_pT0<1.0 && (b_PRO && (fabs(d_TPCnSigmaProton) < 2.0)))
            // ||(d_pT0>=1.0 && (b_PRO && (fabs(d_TPCnSigmaProton) < 2.0)) && d_tofBeta0>0 && (mass2 > 0.7) && (mass2 < 1.1))
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)
            )
            {
              double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

              d_v1_sub1 = d_v1_obs_itrk1/d_resolution1_full_apx;
              tp1_PRO_v1_vs_y_sub1->Fill(d_y_PRO,d_v1_sub1);
              tp1_PRO_v1_vs_y->Fill(d_y_PRO,d_v1_sub1);

              tp1_PRO_v1_vs_pT_sub1->Fill(d_pT0,d_v1_sub1);
              tp1_PRO_v1_vs_pT->Fill(d_pT0,d_v1_sub1);

              h_PRO_plus_pT -> Fill(d_pT0);
              h_PRO_plus_y  -> Fill(d_y_PRO);
              h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
              h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));

              h2_PRO_plus_eta_vs_y->Fill(d_y_PRO,d_eta);
              h2_PRO_plus_eta_vs_y_sub1->Fill(d_y_PRO,d_eta);

              // h2_PRO_plus_v1_vs_y_sub1  -> Fill(d_y_PRO,d_v1_sub1);
              // h2_PRO_plus_v1_vs_pT_sub1 -> Fill(d_pT0,d_v1_sub1);
              // h2_PRO_plus_v1_vs_y       -> Fill(d_y_PRO,d_v1_sub1);
              // h2_PRO_plus_v1_vs_pT      -> Fill(d_pT0,d_v1_sub1);
          }
        }//if(d_charge > 0.0)
      else if(d_eta<=-1.17)
        {
          // if( b_E
          //   && (fabs(d_TPCnSigmaElectron) < 3.0)
          //   && (mass2 < 0.005)
          //   && (mass2 > -0.008)
          //   && (d_pT0 < 0.3)){
          //
          // }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35)){
          //
          // }
          if(
            // b_PI
            // && (fabs(d_TPCnSigmaPion) < 2.0)
            // && (mass2 > -0.03)
            // && (mass2 < 0.06)
            // && ((d_pT0 > 0.4) || (mass2 > 0.008))

            (trait==NULL && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0))
                ||(trait && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 && (mass2 > -0.15) && (mass2 < 0.12))

            //(d_pT0<1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) )
            //||(d_pT0>=1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 &&  ((d_pT0 > 0.4) || (mass2 > 0.008)) && (mass2 > -0.03) && (mass2 < 0.06))
          ){
            double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

            d_v1_sub1 = d_v1_obs_itrk1/d_resolution1_full_apx;
            tp1_PI_minus_v1_vs_y_sub1->Fill(d_y_PI,d_v1_sub1);
            tp1_PI_minus_v1_vs_y->Fill(d_y_PI,d_v1_sub1);

            tp1_PI_minus_v1_vs_pT_sub1->Fill(d_pT0,d_v1_sub1);
            tp1_PI_minus_v1_vs_pT->Fill(d_pT0,d_v1_sub1);

            h_PI_minus_pT -> Fill(d_pT0);
            h_PI_minus_y  -> Fill(d_y_PI);
            h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

            h2_PI_minus_eta_vs_y->Fill(d_y_PI,d_eta);
            h2_PI_minus_eta_vs_y_sub1->Fill(d_y_PI,d_eta);
          }
          if( b_PRO
            && (fabs(d_TPCnSigmaProton) < 3.0)
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)
            )
            {

              // h2_PRO_minus_v1_vs_y_sub1  -> Fill(d_y_PRO,d_v1_sub1);
              // h2_PRO_minus_v1_vs_pT_sub1 -> Fill(d_pT0,d_v1_sub1);
              // h2_PRO_minus_v1_vs_y       -> Fill(d_y_PRO,d_v1_sub1);
              // h2_PRO_minus_v1_vs_pT      -> Fill(d_pT0,d_v1_sub1);
          }
        }// sub 1

      if(d_charge > 0.0 && d_eta>-1.17)
        {
          // if( b_E
          //   && (fabs(d_TPCnSigmaElectron) < 3.0)
          //   && (mass2 < 0.005)
          //   && (mass2 > -0.008)
          //   && (d_pT0 < 0.3)){
          //
          // }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35)){
          //
          // }
          if(
            // b_PI
            // && (fabs(d_TPCnSigmaPion) < 2.0)
            // && (mass2 > -0.03)
            // && (mass2 < 0.06)
            // && ((d_pT0 > 0.4) || (mass2 > 0.008))

            (trait==NULL && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0))
                ||(trait && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 && (mass2 > -0.15) && (mass2 < 0.12))

            //  (d_pT0<1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) )
            // ||(d_pT0>=1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 &&  ((d_pT0 > 0.4) || (mass2 > 0.008)) && (mass2 > -0.03) && (mass2 < 0.06))
          ){
            double d_v1_obs_itrk2 = TMath::Cos(d_Phi - d_Psi1_flat1);

            d_v1_sub2 = d_v1_obs_itrk2/d_resolution1_full_apx;
            tp1_PI_plus_v1_vs_y_sub2->Fill(d_y_PI,d_v1_sub2);
            tp1_PI_plus_v1_vs_y->Fill(d_y_PI,d_v1_sub2);

            tp1_PI_plus_v1_vs_pT_sub2->Fill(d_pT0,d_v1_sub2);
            tp1_PI_plus_v1_vs_pT->Fill(d_pT0,d_v1_sub2);

            h_PI_plus_pT -> Fill(d_pT0);
            h_PI_plus_y  -> Fill(d_y_PI);
            h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

            h2_PI_plus_eta_vs_y->Fill(d_y_PI,d_eta);
            h2_PI_plus_eta_vs_y_sub2->Fill(d_y_PI,d_eta);
          }
          if(
            // b_PRO
            // && (fabs(d_TPCnSigmaProton) < 2.0)
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)

            (trait==NULL && d_pT0>0.4 && d_pT0<2.0 && (fabs(d_TPCnSigmaProton) < 2.0))
                ||(trait && d_pT0>0.4 && d_pT0<2.0 && (fabs(d_TPCnSigmaProton) < 2.0) && d_tofBeta0>0 && (mass2 > 0.4) && (mass2 < 1.4))

            // (d_pT0<1.0 && (b_PRO && (fabs(d_TPCnSigmaProton) < 2.0)))
            // ||(d_pT0>=1.0 && (b_PRO && (fabs(d_TPCnSigmaProton) < 2.0)) && d_tofBeta0>0 && (mass2 > 0.7) && (mass2 < 1.1))
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)
            )
            {
              double d_v1_obs_itrk2 = TMath::Cos(d_Phi - d_Psi1_flat1);

              d_v1_sub2 = d_v1_obs_itrk2/d_resolution1_full_apx;
              tp1_PRO_v1_vs_y_sub2->Fill(d_y_PRO,d_v1_sub2);
              tp1_PRO_v1_vs_y->Fill(d_y_PRO,d_v1_sub2);

              tp1_PRO_v1_vs_pT_sub2->Fill(d_pT0,d_v1_sub2);
              tp1_PRO_v1_vs_pT->Fill(d_pT0,d_v1_sub2);

              h_PRO_plus_pT -> Fill(d_pT0);
              h_PRO_plus_y  -> Fill(d_y_PRO);
              h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
              h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));

              h2_PRO_plus_eta_vs_y->Fill(d_y_PRO,d_eta);
              h2_PRO_plus_eta_vs_y_sub2->Fill(d_y_PRO,d_eta);

              // h2_PRO_plus_v1_vs_y_sub2  -> Fill(d_y_PRO,d_v1_sub2);
              // h2_PRO_plus_v1_vs_pT_sub2 -> Fill(d_pT0,d_v1_sub2);
              // h2_PRO_plus_v1_vs_y       -> Fill(d_y_PRO,d_v1_sub2);
              // h2_PRO_plus_v1_vs_pT      -> Fill(d_pT0,d_v1_sub2);
          }
        }//if(d_charge > 0.0)
      else if(d_eta>-1.17)
        {
          // if( b_E
          //   && (fabs(d_TPCnSigmaElectron) < 3.0)
          //   && (mass2 < 0.005)
          //   && (mass2 > -0.008)
          //   && (d_pT0 < 0.3)){
          //
          // }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35)){
          //
          // }
          if(
            // b_PI
            // && (fabs(d_TPCnSigmaPion) < 2.0)
            // && (mass2 > -0.03)
            // && (mass2 < 0.06)
            // && ((d_pT0 > 0.4) || (mass2 > 0.008))

            (trait==NULL && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0))
                ||(trait && d_pT0>0.2 && d_pT0<1.6 && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 && (mass2 > -0.15) && (mass2 < 0.12))

            // (d_pT0<1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) )
            // ||(d_pT0>=1.0 && b_PI && (fabs(d_TPCnSigmaPion) < 2.0) && d_tofBeta0>0 &&  ((d_pT0 > 0.4) || (mass2 > 0.008)) && (mass2 > -0.03) && (mass2 < 0.06))
          ){
            double d_v1_obs_itrk2 = TMath::Cos(d_Phi - d_Psi1_flat1);

            d_v1_sub2 = d_v1_obs_itrk2/d_resolution1_full_apx;
            tp1_PI_minus_v1_vs_y_sub2->Fill(d_y_PI,d_v1_sub2);
            tp1_PI_minus_v1_vs_y->Fill(d_y_PI,d_v1_sub2);

            tp1_PI_minus_v1_vs_pT_sub2->Fill(d_pT0,d_v1_sub2);
            tp1_PI_minus_v1_vs_pT->Fill(d_pT0,d_v1_sub2);

            h_PI_minus_pT -> Fill(d_pT0);
            h_PI_minus_y  -> Fill(d_y_PI);
            h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

            h2_PI_minus_eta_vs_y->Fill(d_y_PI,d_eta);
            h2_PI_minus_eta_vs_y_sub2->Fill(d_y_PI,d_eta);
          }
          if( b_PRO
            && (fabs(d_TPCnSigmaProton) < 3.0)
            // && (mass2 > 0.7)
            // && (mass2 < 1.1)
            )
            {
              // h2_PRO_minus_v1_vs_y_sub2  -> Fill(d_y_PRO,d_v1_sub2);
              // h2_PRO_minus_v1_vs_pT_sub2 -> Fill(d_pT0,d_v1_sub2);
              // h2_PRO_minus_v1_vs_y       -> Fill(d_y_PRO,d_v1_sub2);
              // h2_PRO_minus_v1_vs_pT      -> Fill(d_pT0,d_v1_sub2);
          }
        }//sub 2
        //======================= END PID ========================================
        if(d_tofBeta == -999) continue;
        // TOF Beta Cut

    }

    //================ END Observed Flow Track Loop ============================



    //==================== END Track Analysis ==================================

  }
  //======================== END Event Loop ====================================

  outputFile->cd();

  h_eta ->SetTitle("Pseudo Rapidity");
  h_eta ->GetXaxis()->SetTitle("eta");

  h_eta0->SetTitle("Subevents_1 Pseudo Rapidity");
  h_eta0->GetXaxis()->SetTitle("eta");

  h_eta1->SetTitle("Subevents_2 Pseudo Rapidity");
  h_eta1->GetXaxis()->SetTitle("eta");

  h_Psi1_sub1->SetTitle("Subevents_1 Psi_{1}");
  h_Psi1_sub1->GetXaxis()->SetTitle("Psi_{1}");
  h_Psi1_sub1->GetYaxis()->SetTitle("N_{events}");

  h_Psi1_sub2->SetTitle("Subevents_2 Psi_{1}");
  h_Psi1_sub2->GetXaxis()->SetTitle("Psi_{1}");
  h_Psi1_sub2->GetYaxis()->SetTitle("N_{events}");

  h_Psi1_rec1->SetTitle("Subevents_1 Recentered Psi_{1}");
  h_Psi1_rec1->GetXaxis()->SetTitle("Psi_{1rec}");
  h_Psi1_rec1->GetYaxis()->SetTitle("N_{events}");


  h_Psi1_rec2->SetTitle("Subevents_2 Recentered Psi_{1}");
  h_Psi1_rec2->GetXaxis()->SetTitle("Psi_{1rec}");
  h_Psi1_rec2->GetYaxis()->SetTitle("N_{events}");

  h_Psi1_flat1->SetTitle("Subevents_1 Flattened Psi_{1}");
  h_Psi1_flat1->GetXaxis()->SetTitle("Psi_{1flat}");
  h_Psi1_flat1->GetYaxis()->SetTitle("N_{events}");


  h_Psi1_flat2->SetTitle("Subevents_2 Flattened Psi_{1}");
  h_Psi1_flat2->GetXaxis()->SetTitle("Psi_{1flat}");
  h_Psi1_flat2->GetYaxis()->SetTitle("N_{events}");

  tp1_PRO_v1_vs_y_sub1->GetXaxis()->SetTitle("y");
  tp1_PRO_v1_vs_y_sub1->GetYaxis()->SetTitle("v1");
  tp1_PRO_v1_vs_y_sub2->GetXaxis()->SetTitle("y");
  tp1_PRO_v1_vs_y_sub2->GetYaxis()->SetTitle("v1");
  tp1_PRO_v1_vs_y->GetXaxis()->SetTitle("y");
  tp1_PRO_v1_vs_y->GetYaxis()->SetTitle("v1");

  tp1_PRO_v1_vs_pT_sub1->GetXaxis()->SetTitle("pT");
  tp1_PRO_v1_vs_pT_sub1->GetYaxis()->SetTitle("v1");
  tp1_PRO_v1_vs_pT_sub2->GetXaxis()->SetTitle("pT");
  tp1_PRO_v1_vs_pT_sub2->GetYaxis()->SetTitle("v1");
  tp1_PRO_v1_vs_pT->GetXaxis()->SetTitle("pT");
  tp1_PRO_v1_vs_pT->GetYaxis()->SetTitle("v1");

  h2_PRO_plus_eta_vs_y->GetXaxis()->SetTitle("y");
  h2_PRO_plus_eta_vs_y_sub1->GetXaxis()->SetTitle("y");
  h2_PRO_plus_eta_vs_y_sub2->GetXaxis()->SetTitle("y");
  h2_PI_plus_eta_vs_y->GetXaxis()->SetTitle("y");
  h2_PI_plus_eta_vs_y_sub1->GetXaxis()->SetTitle("y");
  h2_PI_plus_eta_vs_y_sub2->GetXaxis()->SetTitle("y");
  h2_PI_minus_eta_vs_y->GetXaxis()->SetTitle("y");
  h2_PI_minus_eta_vs_y_sub1->GetXaxis()->SetTitle("y");
  h2_PI_minus_eta_vs_y_sub2->GetXaxis()->SetTitle("y");

  h2_PRO_plus_eta_vs_y->GetYaxis()->SetTitle("eta");
  h2_PRO_plus_eta_vs_y_sub1->GetYaxis()->SetTitle("eta");
  h2_PRO_plus_eta_vs_y_sub2->GetYaxis()->SetTitle("eta");
  h2_PI_plus_eta_vs_y->GetYaxis()->SetTitle("eta");
  h2_PI_plus_eta_vs_y_sub1->GetYaxis()->SetTitle("eta");
  h2_PI_plus_eta_vs_y_sub2->GetYaxis()->SetTitle("eta");
  h2_PI_minus_eta_vs_y->GetYaxis()->SetTitle("eta");
  h2_PI_minus_eta_vs_y_sub1->GetYaxis()->SetTitle("eta");
  h2_PI_minus_eta_vs_y_sub2->GetYaxis()->SetTitle("eta");

  // h_cos_Psi1_flat_diff->SetTitle("Subevents correlations");
  // h_cos_Psi1_flat_diff->GetXaxis()->SetTitle("Cos(Psi_{1flat}^1-Psi_{1flat}^2)");
  // h_cos_Psi1_flat_diff->GetYaxis()->SetTitle("N_{events}");

  // h_Q1_y_sub1->SetTitle("Sub_1 Q1_y");
  // h_Q1_y_sub1->GetXaxis()->SetTitle("Q1_y");
  // h_Q1_x_sub1->SetTitle("Sub_1 Q1_x");
  // h_Q1_x_sub1->GetXaxis()->SetTitle("Q1_x");
  //
  // h_Q1_y_sub2->SetTitle("Sub_2 Q1_y");
  // h_Q1_y_sub2->GetXaxis()->SetTitle("Q1_y");
  // h_Q1_x_sub2->SetTitle("Sub_2 Q1_x");
  // h_Q1_x_sub2->GetXaxis()->SetTitle("Q1_x");

  h_eta ->Write();
  h_eta0->Write();
  h_eta1->Write();
  h_Psi1_sub1->Write();
  h_Psi1_sub2->Write();

  // h_Q1_y_sub1->Write();
  // h_Q1_x_sub1->Write();
  //
  // h_Q1_y_sub2->Write();
  // h_Q1_x_sub2->Write();

  h_Psi1_rec1->Write();
  h_Psi1_rec2->Write();

  h_Psi1_flat1->Write();
  h_Psi1_flat2->Write();
  //
  //
  tp1_PRO_v1_vs_y->Write();
  tp1_PRO_v1_vs_y_sub1->Write();
  tp1_PRO_v1_vs_y_sub2->Write();

  tp1_PRO_v1_vs_pT->Write();
  tp1_PRO_v1_vs_pT_sub1->Write();
  tp1_PRO_v1_vs_pT_sub2->Write();

  tp1_PI_plus_v1_vs_y->Write();
  tp1_PI_plus_v1_vs_y_sub1->Write();
  tp1_PI_plus_v1_vs_y_sub2->Write();

  tp1_PI_plus_v1_vs_pT->Write();
  tp1_PI_plus_v1_vs_pT_sub1->Write();
  tp1_PI_plus_v1_vs_pT_sub2->Write();

  tp1_PI_minus_v1_vs_y->Write();
  tp1_PI_minus_v1_vs_y_sub1->Write();
  tp1_PI_minus_v1_vs_y_sub2->Write();

  tp1_PI_minus_v1_vs_pT->Write();
  tp1_PI_minus_v1_vs_pT_sub1->Write();
  tp1_PI_minus_v1_vs_pT_sub2->Write();
  //
  h_PRO_plus_mT_Diff->Write();
  h_PRO_plus_pT->Write();
  h_PRO_plus_y->Write();
  h2_PRO_plus_pT_vs_y->Write();

  h_PI_plus_mT_Diff->Write();
  h_PI_plus_pT->Write();
  h_PI_plus_y->Write();
  h2_PI_plus_pT_vs_y->Write();

  h_PI_minus_mT_Diff->Write();
  h_PI_minus_pT->Write();
  h_PI_minus_y->Write();
  h2_PI_minus_pT_vs_y->Write();

  h2_PRO_plus_eta_vs_y->Write();
  h2_PRO_plus_eta_vs_y_sub1->Write();
  h2_PRO_plus_eta_vs_y_sub2->Write();
  h2_PI_plus_eta_vs_y->Write();
  h2_PI_plus_eta_vs_y_sub1->Write();
  h2_PI_plus_eta_vs_y_sub2->Write();
  h2_PI_minus_eta_vs_y->Write();
  h2_PI_minus_eta_vs_y_sub1->Write();
  h2_PI_minus_eta_vs_y_sub2->Write();

  // h_cos_Psi1_flat_diff->Write();

  // h_v1_obs_sub1->Write();
  // h_v1_obs_sub2->Write();
  // h_v1_obs     ->Write();
  //
  // h_v1_sub1->Write();
  // h_v1_sub2->Write();
  // h_v1     ->Write();

  // h2_PRO_plus_v1_vs_y_sub1->Write();
  // h2_PRO_plus_v1_vs_y_sub2->Write();
  // h2_PRO_plus_v1_vs_y->Write();
  //
  // h2_PRO_plus_v1_vs_pT_sub1->Write();
  // h2_PRO_plus_v1_vs_pT_sub2->Write();
  // h2_PRO_plus_v1_vs_pT->Write();
  //
  // h2_PRO_minus_v1_vs_y_sub1->Write();
  // h2_PRO_minus_v1_vs_y_sub2->Write();
  // h2_PRO_minus_v1_vs_y->Write();
  //
  // h2_PRO_minus_v1_vs_pT_sub1->Write();
  // h2_PRO_minus_v1_vs_pT_sub2->Write();
  // h2_PRO_minus_v1_vs_pT->Write();

  // h_Q1_y_diff1->Write();
  // h_Q1_x_diff1->Write();
  //
  // h_Q1_y_diff2->Write();
  // h_Q1_x_diff2->Write();

  // for(int i=0; i<20;i++)
  // {
  //   h_sin_iPsi1_rec1[i]->Write();
  //   h_cos_iPsi1_rec1[i]->Write();
  //
  //   h_sin_iPsi1_rec2[i]->Write();
  //   h_cos_iPsi1_rec2[i]->Write();
  //
  // }
}
//================================ END Main Funcition ==========================
