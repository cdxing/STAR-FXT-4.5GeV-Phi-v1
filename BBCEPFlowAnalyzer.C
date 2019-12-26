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
void BBCEPFlowAnalyzer(const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
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

  // Read BBC EP from Yang's Event Plane Root file
  TFile *f1 = new TFile("/star/data01/pwg/rseto/eventPlanesPicoDST/eventPlanes.picoDst.result1.root","read");
	if(f1->IsOpen()) cout<<"Event planes loaded successfully!"<<endl;
	else {cout<<"No event planes, quit!"<<endl; return;}
	TTree *t1 = (TTree*)f1->Get("EPTree");
	if(!t1) return;
	Int_t nEntries = t1->GetEntries();
	cout<<nEntries<<" total events found!"<<endl;

	// TFile *o1 = new TFile("testPlots.root","recreate");
	Float_t PI = TMath::Pi();
	TH1D *h1 = new TH1D("runId","Run Id",10,16140031,16140041);
	TH1D *h2 = new TH1D("eventId","Event Id",1e3,0,7.e5);
	TH1D *h3 = new TH1D("mult","Multiplicity",250,0,250);
	TH1D *h4 = new TH1D("bbcEast","#psi^{BBC}_{East} [Radian]",500,-0.5*PI,2.5*PI);
	TH1D *h5 = new TH1D("tpcEast","#psi^{TPC}_{East} [Radian]",500,-0.5*PI,2.5*PI);
	TH1D *h6 = new TH1D("tpcWest","#psi^{TPC}_{West} [Radian]",500,-0.5*PI,2.5*PI);
	TH1D *h7 = new TH1D("trackId","Track Id",4000,0,4000);
	TH1D *h8 = new TH1D("tpcEastInd","#psi^{TPC}_{East} for each track [Radian]",500,-0.5*PI,2.5*PI);
	TH1D *h9 = new TH1D("tpcWestInd","#psi^{TPC}_{West} for each track [Radian]",500,-0.5*PI,2.5*PI);

	TLeaf *leaf_runId = t1->GetLeaf("RunId");
	TLeaf *leaf_eventId = t1->GetLeaf("EventId");
	TLeaf *leaf_mult = t1->GetLeaf("Multiplicity");
	TLeaf *leaf_trackId = t1->GetLeaf("TrackId");
	TLeaf *leaf_TPCeastEP = t1->GetLeaf("TPCeastEP");
	TLeaf *leaf_TPCwestEP = t1->GetLeaf("TPCeastEP");
	TLeaf *leaf_TPCeastIndEP = t1->GetLeaf("TPCeastIndEP");
	TLeaf *leaf_TPCwestIndEP = t1->GetLeaf("TPCeastIndEP");
	TLeaf *leaf_BBCeastEP = t1->GetLeaf("BBCeastEP");
  // END Read BBC EP from Yang's Event Plane Root file

  // Read a map of # of entry vs runId vs event Id
  TFile *f2 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/eventPlanesPicoDST/testPlots.root","read");
	if(f1->IsOpen()) cout<<"ID map loaded successfully!"<<endl;
	else {cout<<"No ID maps, quit!"<<endl; return;}
	TH2D *h2_id = (TH2D*)f2->Get("h2_1");
	if(!h2_id) return;
  // END Read a map of # of entry vs runId vs event Id




  double  d_zvtx  = -9999.0;
  std::vector <unsigned int> triggerIDs;
  // Event cut set up

  // Histogramning

  // End Test Plots

  // Directed flow

  TH2D *  h2_PRO_plus_v1_vs_y      = new TH2D("h2_PRO_plus_v1_vs_y","h2_PRO_plus_v1_vs_y",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_plus_v1_vs_pT      = new TH2D("h2_PRO_plus_v1_vs_pT","h2_PRO_plus_v1_vs_pT",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_y      = new TH2D("h2_PRO_minus_v1_vs_y","h2_PRO_minus_v1_vs_y",1001,-5.005,5.005,1001,-5.005,5.005);
  TH2D *  h2_PRO_minus_v1_vs_pT      = new TH2D("h2_PRO_minus_v1_vs_pT","h2_PRO_minus_v1_vs_pT",1001,-5.005,5.005,1001,-5.005,5.005);

  TProfile * tp1_PRO_v1_vs_y_cent[7];
  TProfile * tp1_PRO_v1_vs_pT_cent[7];
  TProfile * tp1_PI_plus_v1_vs_y_cent[7];
  TProfile * tp1_PI_plus_v1_vs_pT_cent[7];
  TProfile * tp1_PI_minus_v1_vs_y_cent[7];
  TProfile * tp1_PI_minus_v1_vs_pT_cent[7];

  for(Int_t i=0;i<7;i++)
  {
   tp1_PRO_v1_vs_y_cent[i] = new TProfile(Form("tp1_PRO_v1_vs_y_cent[%d]",i),Form("tp1_PRO_v1_vs_y_cent[%d]",i),15,-3,0,-5.0,5.0);
   tp1_PRO_v1_vs_pT_cent[i] = new TProfile(Form("tp1_PRO_v1_vs_pT_cent[%d]",i),Form("tp1_PRO_v1_vs_pT_cent[%d]",i),50,0,5,-5.0,5.0);
   tp1_PI_plus_v1_vs_y_cent[i] = new TProfile(Form("tp1_PI_plus_v1_vs_y_cent[%d]",i),Form("tp1_PI_plus_v1_vs_y_cent[%d]",i),15,-3,0,-5.0,5.0);
   tp1_PI_plus_v1_vs_pT_cent[i] = new TProfile(Form("tp1_PI_plus_v1_vs_pT_cent[%d]",i),Form("tp1_PI_plus_v1_vs_pT_cent[%d]",i),50,0,5,-5.0,5.0);
   tp1_PI_minus_v1_vs_y_cent[i] = new TProfile(Form("tp1_PI_minus_v1_vs_y_cent[%d]",i),Form("tp1_PI_minus_v1_vs_y_cent[%d]",i),15,-3,0,-5.0,5.0);
   tp1_PI_minus_v1_vs_pT_cent[i] = new TProfile(Form("tp1_PI_minus_v1_vs_pT_cent[%d]",i),Form("tp1_PI_minus_v1_vs_pT_cent[%d]",i),50,0,5,-5.0,5.0);
  }
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
  TH2D *  h2_PI_plus_eta_vs_y  = new TH2D("h2_PI_plus_eta_vs_y","h2_PI_plus_eta_vs_y",1000,-5.0,0.0,1000,-5.0,0.0);
  TH2D *  h2_PI_minus_eta_vs_y  = new TH2D("h2_PI_minus_eta_vs_y","h2_PI_minus_eta_vs_y",1000,-5.0,0.0,1000,-5.0,0.0);


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
    if( !event )
    {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    Int_t runId_Dst   = event->runId();
    Int_t eventId_Dst = event->eventId();
    Int_t entrynumber = h2_id->GetBinContent(runId_Dst-16140030,eventId_Dst+1);
    Int_t runId = -999;
    Int_t eventId = -999;

    double BBCeastEP = -999;
    // for(Int_t ievent = 0; ievent < nEntries; ievent++)
    // {
      t1->GetEntry(entrynumber);
      runId = leaf_runId->GetValue(0);
      eventId = leaf_eventId->GetValue(0);

      if(runId != runId_Dst) continue;
      if(eventId != eventId_Dst) continue;

      if(runId == runId_Dst && eventId == eventId_Dst)
      {
        BBCeastEP = leaf_BBCeastEP->GetValue(0);
      }

    // }
    if(runId != runId_Dst || eventId != eventId_Dst || runId == -999 || eventId == -999 || BBCeastEP == -999) continue;
    std::cout << "ID Matched" << std::endl;


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

    // double d_v1             = -999;
    // END Event Plane Set Up

    //======================== Track Analysis ==================================
    Int_t nTracks = dst->numberOfTracks();
    //Number of tracks in this event
    double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut // Mon Jul  3 09:16:46 EDT 2017
    // if(nTracks <= 10) continue;

    //=================== Centrality Flow Track Loop ===========================
    bool b_cent[7] = {false};
    int nGoodTracks = 0;
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);

      if(!picoTrack) continue;

      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      bool b_bad_track    = b_bad_dEdx || b_bad_tracking ;
      if(b_bad_track) continue;

      //============= Primary Track Cut ===================
      if(!picoTrack->isPrimary()) continue;

      StPicoBTofPidTraits *trait = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      double d_tofBeta0        = -999;
      if(trait) d_tofBeta0 = trait->btofBeta();
      nGoodTracks++;
      if(d_tofBeta0 == -999) continue;
    }
    bool b_pileup   = (nGoodTracks >  240);
    b_cent[0]       = (nGoodTracks >= 153);
    b_cent[1]       = (nGoodTracks >= 121 && nGoodTracks < 153);
    b_cent[2]       = (nGoodTracks >= 97  && nGoodTracks < 121);
    b_cent[3]       = (nGoodTracks >= 77  && nGoodTracks < 97);
    b_cent[4]       = (nGoodTracks >= 61  && nGoodTracks < 77);
    b_cent[5]       = (nGoodTracks >= 48  && nGoodTracks < 61);
    b_cent[6]       = (nGoodTracks < 48);
    bool b_low_mult = (nGoodTracks <= 25);


    //=============== END Centrality Flow Track Loop ===========================


    //=================== Observed Flow Track Loop =============================
    // FLow Track Set Up

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

      for(Int_t i=0; i<7;i++)
      {
        if(b_cent[i])
        {
          if(d_charge > 0.0 )
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
                // double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

                double d_v1 = TMath::Cos(d_Phi - BBCeastEP);

                tp1_PI_plus_v1_vs_y_cent[i]->Fill(d_y_PI,d_v1);
                tp1_PI_plus_v1_vs_pT_cent[i]->Fill(d_pT0,d_v1);

                h_PI_plus_pT -> Fill(d_pT0);
                h_PI_plus_y  -> Fill(d_y_PI);
                h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

                h2_PI_plus_eta_vs_y->Fill(d_y_PI,d_eta);
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
                  // double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

                  double d_v1 = TMath::Cos(d_Phi - BBCeastEP);//d_v1_obs_itrk1/d_resolution1_full_apx;


                  tp1_PRO_v1_vs_y_cent[i]->Fill(d_y_PRO,d_v1);
                  tp1_PRO_v1_vs_pT_cent[i]->Fill(d_pT0,d_v1);

                  h_PRO_plus_pT -> Fill(d_pT0);
                  h_PRO_plus_y  -> Fill(d_y_PRO);
                  h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
                  h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));

                  h2_PRO_plus_eta_vs_y->Fill(d_y_PRO,d_eta);
              }
            }//if(d_charge > 0.0)
          else
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
                // double d_v1_obs_itrk1 = TMath::Cos(d_Phi - d_Psi1_flat2);

                double d_v1 = TMath::Cos(d_Phi - BBCeastEP);//d_v1_obs_itrk1/d_resolution1_full_apx;


                tp1_PI_minus_v1_vs_y_cent[i]->Fill(d_y_PI,d_v1);
                tp1_PI_minus_v1_vs_pT_cent[i]->Fill(d_pT0,d_v1);

                h_PI_minus_pT -> Fill(d_pT0);
                h_PI_minus_y  -> Fill(d_y_PI);
                h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));

                h2_PI_minus_eta_vs_y->Fill(d_y_PI,d_eta);
              }
              if( b_PRO
                && (fabs(d_TPCnSigmaProton) < 3.0)
                // && (mass2 > 0.7)
                // && (mass2 < 1.1)
                )
                {


              }
            }
        }
      }


        //======================= END PID ========================================
        // if(d_tofBeta == -999) continue;
        // TOF Beta Cut

    }

    //================ END Observed Flow Track Loop ============================



    //==================== END Track Analysis ==================================

  }
  //======================== END Event Loop ====================================

  outputFile->cd();

  // tp1_PRO_v1_vs_y->GetXaxis()->SetTitle("y");
  // tp1_PRO_v1_vs_y->GetYaxis()->SetTitle("v1");
  //
  // tp1_PRO_v1_vs_pT->GetXaxis()->SetTitle("pT");
  // tp1_PRO_v1_vs_pT->GetYaxis()->SetTitle("v1");

  h2_PRO_plus_eta_vs_y->GetXaxis()->SetTitle("y");
  h2_PI_plus_eta_vs_y->GetXaxis()->SetTitle("y");
  h2_PI_minus_eta_vs_y->GetXaxis()->SetTitle("y");

  h2_PRO_plus_eta_vs_y->GetYaxis()->SetTitle("eta");
  h2_PI_plus_eta_vs_y->GetYaxis()->SetTitle("eta");
  h2_PI_minus_eta_vs_y->GetYaxis()->SetTitle("eta");

  for(Int_t i=0;i<7;i++)
  {
  tp1_PRO_v1_vs_y_cent[i] ->Write();
  tp1_PRO_v1_vs_pT_cent[i] ->Write();
  tp1_PI_plus_v1_vs_y_cent[i] ->Write();
  tp1_PI_plus_v1_vs_pT_cent[i] ->Write();
  tp1_PI_minus_v1_vs_y_cent[i] ->Write();
  tp1_PI_minus_v1_vs_pT_cent[i] ->Write();
  }

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
  h2_PI_plus_eta_vs_y->Write();
  h2_PI_minus_eta_vs_y->Write();

}
//================================ END Main Funcition ==========================
