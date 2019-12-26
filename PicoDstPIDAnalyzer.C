/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 *
 * \brief Incoporated Jason's FXT 4.5 Phi analysis code
 *
 * \author Yang Wu
 * \date May 15, 2019
 *
 */

// This is needed for calling standalone classes (not needed on RACF)
//#define _VANILLA_ROOT_

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

//using namespace std;

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
//#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
//R__LOAD_LIBRARY(StPicoEvent/libStPicoDst)
//#endif

// inFile - is a name of name.picoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.picoDst.root files


class st_track : public TObject
{
public:
    st_track(){ reset(); }
    virtual ~st_track(){}
    bool b_pos_charge;
    /* bool b_bad_TOF_match; */
    /* bool b_bad_ToF; */

    bool b_PI;
    bool b_PRO;
    bool b_K;

    /* unsigned int i_event; */
    /* unsigned int nGoodTracks; */
    unsigned int runNumber;
    unsigned int eventNumber;

    double MagField;

    float px;
    float py;
    float pz;

    float x_vtx;
    float y_vtx;
    float z_vtx;



    StPicoPhysicalHelix trackhelix;

    void reset()
    {
      b_pos_charge = 0;

      b_PI  = 0;
      b_PRO = 0;
      b_K   = 0;

      /* i_event     = 0; */

      runNumber   = 0;
      eventNumber = 0;

      MagField = 0.0;

      px = 0.0;
      py = 0.0;
      pz = 0.0;

      x_vtx = 0.0;
      y_vtx = 0.0;
      z_vtx = 0.0;

    }//void reset()


private:
    ClassDef(st_track,1);
    //    ClassDef(st_track,1);
};//  class st_track : public TObject

//st_K
class st_K : public TObject
{
public:
    st_K(){ reset(); }
    virtual ~st_K(){}
    bool b_pos_charge;
    bool b_bad_TOF;
    /* bool b_bad_TOF_match; */
    /* bool b_bad_ToF; */

    unsigned int runNumber;
    unsigned int eventNumber;
    unsigned int nGoodTracks;
    unsigned int nTrkvsCuts;


    float px;
    float py;
    float pz;

    float x_vtx;
    float y_vtx;
    float z_vtx;

    float tofBeta;
    float dEdxFit;
    float d_dEdxTru;
    /* float eta; */
    /* float phi; */
    float DCA_r;
    /* float d_zvtx; */

    /* unsigned int nhits0; */

    float d_TPCnSigmaKaon;
    float d_TOFnSigmaKaon;

    void reset()
    {
      b_pos_charge = 0;
      b_bad_TOF = 0;
      /* i_event     = 0; */

      runNumber   = 0;
      eventNumber = 0;
      nGoodTracks = 0;
      nTrkvsCuts  = 0;

      px = 0.0;
      py = 0.0;
      pz = 0.0;

      x_vtx = 0.0;
      y_vtx = 0.0;
      z_vtx = 0.0;

      d_TPCnSigmaKaon = 0.0;
      d_TOFnSigmaKaon = 0.0;

      DCA_r = -9999.0;
    }//void reset()

private:
    ClassDef(st_K,1);
    //    ClassDef(st_K,1);
};//  class st_K : public TObject
  //ENd DEF st_K


//_________________
void PicoDstPIDAnalyzer(const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test"
                      ) {
  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

    //gROOT->ProcessLine( ".L StPicoEvent/libStPicoDst.so" );

    //gSystem->Load(".sl73_gcc485/lib/libStPicoEvent.so");

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  //Long64_t events2read = picoReader->chain()->GetEntries();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);
  //picoReader->SetStatus("BTowHit",0);
  //picoReader->SetStatus("EmcTrigger",0);
  //picoReader->SetStatus("TrackCovMatrix",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

    if( !picoReader->chain() ) {
        std::cout << "No chain has been found." << std::endl;
    }
    Long64_t eventsInTree = picoReader->tree()->GetEntries();
    std::cout << "eventsInTree: "  << eventsInTree << std::endl;
    Long64_t events2read = picoReader->chain()->GetEntries();

    std::cout << "Number of events to read: " << events2read
    << std::endl;

    outFile.Append(".picoDst.result.root");

    TFile *outputFile = new TFile(outFile,"recreate");

    // Histogramming
    // Histogramming

    double  d_zvtx  = -9999.0;
    unsigned int  i_event = 0;//-1;

    unsigned int  runNumber        = 0;
    unsigned int  eventNumber      = 0;

    std::vector <unsigned int> triggerIDs;

    //Create the outout Root File
    //TFile * outFile = new TFile(fileName.Data(), "RECREATE");
    //outFile = new TFile(d_OutputFileName.c_str(), "RECREATE");
    //TFile * outFile = new TFile("testOutput.root", "RECREATE");

    st_K Kaoninfo;
    st_track trackinfo;

    TTree *  t_K = new TTree("t_K","t_K");
    t_K -> Branch("Kaoninfo",&Kaoninfo);

    TTree *  t_tracks   = new TTree("t_tracks","t_tracks");
    t_tracks   -> Branch("trackinfo",&trackinfo);



    double d_PI_m2   = 0.019479835;
    double d_PRO_m2  = 0.880354;
    double d_K_m2    = 0.24371698;

    TH2D *  h2_m2_QA_pq = new TH2D("h2_m2_QA_pq","h2_m2_QA_pq",5000,-5.0,5.0,5000,-0.2,5.0);
    TH2D *  h2_m2_QA_pq_1 = new TH2D("h2_m2_QA_pq_1","h2_m2_QA_pq_1",5000,-5.0,5.0,5000,-0.2,5.0);
    TH2D *  h2_m2_QA_pT = new TH2D("h2_m2_QA_pT","h2_m2_QA_pT",5000,-5.0,5.0,5000,-2.0,5.0);
    TH2D *  h2_m2_QA_pT_1 = new TH2D("h2_m2_QA_pT_1","h2_m2_QA_pT_1",5000,-5.0,5.0,5000,-2.0,5.0);

    TH1D *  h_E_plus_mT_Diff   = new TH1D("h_E_plus_mT_Diff","h_E_plus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_K_plus_mT_Diff   = new TH1D("h_K_plus_mT_Diff","h_K_plus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_PI_plus_mT_Diff  = new TH1D("h_PI_plus_mT_Diff","h_PI_plus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_PRO_plus_mT_Diff = new TH1D("h_PRO_plus_mT_Diff","h_PRO_plus_mT_Diff",1000,0.0,5.0);

    TH1D *  h_E_minus_mT_Diff   = new TH1D("h_E_minus_mT_Diff","h_E_minus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_K_minus_mT_Diff   = new TH1D("h_K_minus_mT_Diff","h_K_minus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_PI_minus_mT_Diff  = new TH1D("h_PI_minus_mT_Diff","h_PI_minus_mT_Diff",1000,0.0,5.0);
    TH1D *  h_PRO_minus_mT_Diff = new TH1D("h_PRO_minus_mT_Diff","h_PRO_minus_mT_Diff",1000,0.0,5.0);

    TH1D *  h_E_plus_pT   = new TH1D("h_E_plus_pT","Electron+ pT",1000,0.0,5.0);
    TH1D *  h_K_plus_pT   = new TH1D("h_K_plus_pT","Kaon+ pT",1000,0.0,5.0);
    TH1D *  h_PI_plus_pT  = new TH1D("h_PI_plus_pT","Pion+ pT",1000,0.0,5.0);
    TH1D *  h_PRO_plus_pT = new TH1D("h_PRO_plus_pT","Proton+ pT",1000,0.0,5.0);

    TH1D *  h_E_minus_pT   = new TH1D("h_E_minus_pT","Electron- pT",1000,0.0,5.0);
    TH1D *  h_K_minus_pT   = new TH1D("h_K_minus_pT","Kaon- pT",1000,0.0,5.0);
    TH1D *  h_PI_minus_pT  = new TH1D("h_PI_minus_pT","Pion- pT",1000,0.0,5.0);
    TH1D *  h_PRO_minus_pT = new TH1D("h_PRO_minus_pT","Proton- pT",1000,0.0,5.0);

    TH1D *  h_E_plus_y   = new TH1D("h_E_plus_y","Electron+ y",1000,-5.0,0.0);
    TH1D *  h_K_plus_y   = new TH1D("h_K_plus_y","Kaon+ y",1000,-5.0,0.0);
    TH1D *  h_PI_plus_y  = new TH1D("h_PI_plus_y","Pion+ y",1000,-5.0,0.0);
    TH1D *  h_PRO_plus_y = new TH1D("h_PRO_plus_y","Proton+ y",1000,-5.0,0.0);

    TH1D *  h_E_minus_y   = new TH1D("h_E_minus_y","Electron- y",1000,-5.0,0.0);
    TH1D *  h_K_minus_y   = new TH1D("h_K_minus_y","Kaon- y",1000,-5.0,0.0);
    TH1D *  h_PI_minus_y  = new TH1D("h_PI_minus_y","Pion- y",1000,-5.0,0.0);
    TH1D *  h_PRO_minus_y = new TH1D("h_PRO_minus_y","Proton- y",1000,-5.0,0.0);

    TH2D *  h2_E_plus_pT_vs_y   = new TH2D("h2_E_plus_pT_vs_y","h2_E_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_K_plus_pT_vs_y   = new TH2D("h2_K_plus_pT_vs_y","h2_K_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PI_plus_pT_vs_y  = new TH2D("h2_PI_plus_pT_vs_y","h2_PI_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PRO_plus_pT_vs_y = new TH2D("h2_PRO_plus_pT_vs_y","h2_PRO_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);

    TH2D *  h2_E_minus_pT_vs_y   = new TH2D("h2_E_minus_pT_vs_y","h2_E_minus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_K_minus_pT_vs_y   = new TH2D("h2_K_minus_pT_vs_y","h2_K_minus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PI_minus_pT_vs_y  = new TH2D("h2_PI_minus_pT_vs_y","h2_PI_minus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PRO_minus_pT_vs_y = new TH2D("h2_PRO_minus_pT_vs_y","h2_PRO_minus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);


    // //Set the Output File Name
    //TString fileName(d_OutputFileName);
    //fileName.Append(mFileIndex);
    //fileName.Append(".root");
    //  fileName.Prepend(mOutDir);

    // Loop over events
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    //for(Long64_t iEvent=0; iEvent<10; iEvent++) {

        std::cout << "Working on event #[" << (iEvent+1)
        << "/" << events2read << "]" << std::endl;

        Bool_t readEvent = picoReader->readPicoEvent(iEvent);
        if( !readEvent ) {
            std::cout << "Something went wrong, Master! Nothing to analyze..."
            << std::endl;
            break;
        }

        // Retrieve picoDst
        StPicoDst *dst = picoReader->picoDst();

        // Retrieve event information
        StPicoEvent *event = dst->event();
        if( !event ) {
            std::cout << "Something went wrong, Master! Event is hiding from me..."
            << std::endl;
            break;
        }
        //hRefMult->Fill( event->refMult() );


        const Float_t B = event->bField(); //cout<<"Magnetic field is "<<B<<endl;

        int i_muEvent_num = event -> eventId();

        runNumber        = event->runId();
        eventNumber      = event->eventId();

        double d_MagField = event->bField();

        //=========================================================
        //          Trigger Selection
        //=========================================================
        // loop for the trigger ids and see if any == 1
        triggerIDs.clear();
        triggerIDs       = event->triggerIds();
        bool b_bad_trig = true;
        //  cout<<" SIZE : "<<triggerIDs.size()<<endl;
        for(int i=0; i < triggerIDs.size(); i++)
          {
            //      cout<<" trigger: "<<triggerIDs[i]<<endl;
            if(triggerIDs[i] == 1) b_bad_trig = false;
          }//for(int i; i < triggerIDs.size(); i++)
           //bool b_good_fxt_triger = muEvent->triggerIdCollection().nominal().isTrigger(1);
           //=========================================================
           //      END Trigger Selection
           //=========================================================


        TVector3 pVtx = event->primaryVertex();
        //hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
        //hVtxZ->Fill( event->primaryVertex().Z() );

        //bool b_0vtx     = pVtx;

        //------------ Z-VTX Selection ----------------------------
        TVector3 v3D_vtx = event->primaryVertex();
        d_zvtx = pVtx.z();
        bool b_bad_zvtx        = ((d_zvtx < 210.0) || (d_zvtx > 212.0));
        //-------- END Z-VTX Selection ----------------------------


        bool b_bad_evt  = /*b_0vtx ||*/ /*b_low_mult || b_pileup ||*/ b_bad_zvtx /*|| (!b_good_fxt_triger)*/ || b_bad_trig /*|| b_bad_vtxID*/;
        if(b_bad_evt) continue; // Events skipped are not counted in the event count
        // Track analysis
        Int_t nTracks = dst->numberOfTracks();
        Int_t nMatrices = dst->numberOfTrackCovMatrices();
        // if(nTracks != nMatrices) {
        //     //std::cout << "Number of tracks and matrices do not match!" << std::endl;
        // }
        //std::cout << "Number of tracks in event: " << nTracks << std::endl;

        bool b_cent_01  = false; // 0  < centrality <= 5%
        bool b_cent_02  = false; // 5  < centrality <= 10%
        bool b_cent_03  = false; // 10 < centrality <= 15%
        bool b_cent_04  = false; // 15 < centrality <= 20%
        bool b_cent_05  = false; // 20 < centrality <= 25%
        bool b_cent_06  = false; // 25 < centrality <= 30%
        bool b_cent_07  = false; // centrality > 30%

        // # of cuts = 0
        // ievtcut[0] = ievtcut[0] + 1;
        // itrkcut[0] = itrkcut[0] + nTracks;
        // std::cout << "Event(s) with 0 cut(s)" << ievtcut[0] << std::endl;
        // std::cout << "Track(s) with 0 cut(s)" << itrkcut[0] << std::endl;

        // Track loop
        int nGoodTracks = 0;
        int nTrkvsCuts  = 0;
        h2_m2_QA_pq_1   ->Reset();
        h2_m2_QA_pT_1   ->Reset();
        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

            // Retrieve i-th pico track
            StPicoTrack *picoTrack = dst->track(iTrk);

            if(!picoTrack) continue;


            //============== === Track Level Cuts ==============
            bool b_bad_dEdx     = false;
            bool b_bad_tracking = false;

            b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
            b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
            //bool b_not_enough_hits = false;//((double)picoTrack->nHitsFit()) < 15;
            bool b_bad_track    = b_bad_dEdx || b_bad_tracking /*|| b_not_enough_hits*/;
            if(b_bad_track) continue;



            //============== END Track Level Cuts ==============
            //============= Primary Track Cut ===================
            if(!picoTrack->isPrimary()) continue;

            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
            double d_tofBeta0        = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();

            double d_px0      = picoTrack->gMom().x();
            double d_py0      = picoTrack->gMom().y();
            double d_pz0      = picoTrack->gMom().z();
            double d_pT0      = picoTrack->gPt();
            double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
            double mass2      = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

            nGoodTracks++;
            if(d_tofBeta0 == -999) continue;

            h2_m2_QA_pq_1   ->Fill(d_mom0/(picoTrack->charge()),mass2);
            h2_m2_QA_pT_1   ->Fill(d_pT0/(picoTrack->charge()),mass2);
            nTrkvsCuts++;


        }




        bool b_pileup   = (nGoodTracks >  240);
        b_cent_01       = (nGoodTracks >= 153);
        b_cent_02       = (nGoodTracks >= 121 && nGoodTracks < 153);
        b_cent_03       = (nGoodTracks >= 97  && nGoodTracks < 121);
        b_cent_04       = (nGoodTracks >= 77  && nGoodTracks < 97);
        b_cent_05       = (nGoodTracks >= 61  && nGoodTracks < 77);
        b_cent_06       = (nGoodTracks >= 48  && nGoodTracks < 61);
        b_cent_07       = (nGoodTracks < 48);
        bool b_low_mult = (nGoodTracks <= 25);
        //take the low multiplicity cut out for now
        b_low_mult = false;

        // bool b_bad_evt  = /*b_0vtx ||*/ /*b_low_mult || b_pileup ||*/ b_bad_zvtx /*|| (!b_good_fxt_triger)*/ || b_bad_trig /*|| b_bad_vtxID*/;
        // if(b_bad_evt) continue; // Events skipped are not counted in the event count


        h2_m2_QA_pq   ->Add(h2_m2_QA_pq_1);
        h2_m2_QA_pT   ->Add(h2_m2_QA_pT_1);

        //=========================================================
        //      Track Cut Settings
        //=========================================================
        int i_Nhits_min = 9;  //minimum number of hits
                              // double d_pT_min = 0.2;//0.05;
                              // double d_mom_min = 0.2;//0.05;

        //Mon Jul  3 09:17:54 EDT 2017 from phi paper
        double d_pT_min = 0.1;//0.05;
        double d_pT_max = 10.0;//0.05;
        double d_mom_min = 0.1;//0.05;
        double d_mom_max = 10.0;// 1.0;//0.05;
                                //  double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut
        double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut // Mon Jul  3 09:16:46 EDT 2017

        //--- Not being used Fri Jun 16 10:48:30 EDT 2017
        double d_PRO_daughter_DCA = 0.3; // trk must be greater than this dca
        double d_PI_daughter_DCA  = 1.0; // trk must be greater than this dca
                                         //---


        // Lambda: 1.0, K0s 1.2
        //  double d_cut_dca_daughters    = 2.0; //distance between daughters, must be less than this
        double d_cut_dca_daughters_lam    = 1.0;//2.0; //distance between daughters, must be less than this
        double d_cut_dca_daughters_k0s    = 1.2; //distance between daughters, must be less than this


        // Lambda: 1.5 Anti Lambda: 2.0 K0s: 1.0
        //double d_cut_dca_mother       = 5.0; // must be less than this
        double d_cut_dca_mother_lam       = 1.5;//5.0; // must be less than this
        double d_cut_dca_mother_k0s       = 1.0; // must be less than this

        // Lambda: 3.0, K0s: 2.0
        //  double d_cut_mother_decay_length = 3.0; // must be greater than this
        double d_cut_mother_decay_length_lam = 3.0; // must be greater than this
        double d_cut_mother_decay_length_k0s = 2.0; // must be greater than this
        double d_cut_mother_decay_length_RHO = 0.5; // must be LESS    than this
        double d_cut_mother_decay_length_PHI = 0.5; // must be LESS    than this
        //=========================================================
        //  END Track Cut Settings
        //=========================================================

        //=========================================================
        //      Preliminary Track Loop
        //=========================================================
        vector<StPicoTrack *> v_tracks;
        vector<StPicoTrack *> v_tracks_pl;
        vector<StPicoTrack *> v_tracks_mi;
        int index = 0;
        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

            // Retrieve i-th pico track
            StPicoTrack *picoTrack = dst->track(iTrk);

            if(picoTrack == NULL) continue;
            //short flag = picoTrack->flag();
            //if(flag <=0 || flag >=1000 )continue;
            if(picoTrack->isPrimary()) continue;
            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
            //bool tofAvailability = trait;
            unsigned short nHits = picoTrack->nHits();
            if(nHits < i_Nhits_min) continue;
            double d_TPCnSigmaPion   = fabs(picoTrack->nSigmaPion());
            double d_TPCnSigmaProton = fabs(picoTrack->nSigmaProton());
            double d_TPCnSigmaKaon   = fabs(picoTrack->nSigmaKaon());

            //double d_pidProbPion     = picoTrack -> pidProbPion();
            //double d_pidProbProton   = picoTrack -> pidProbProton();
            //double d_pidProbKaon     = picoTrack -> pidProbKaon();

            //double d_TOFnSigmaPion   = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaPion());
            //double d_TOFnSigmaProton = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaProton());
            //double d_TOFnSigmaKaon   = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaKaon());

            double d_tofBeta0        = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();
            //double dEdxFit           = (dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->dEdxFit();
            //double d_dEdxTru         = (dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->dEdxTruncated();

            // PID Cuts
            bool b_PI  = fabs(d_TPCnSigmaPion) < d_SigmaCutLevel;
            bool b_PRO = fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
            //      bool b_K   = fabs(d_TPCnSigmaKaon) < d_SigmaCutLevel;
            bool b_K   = fabs(d_TPCnSigmaKaon) < 4.0;

            //bool b_K_tof   = fabs(d_TOFnSigmaKaon) < d_SigmaCutLevel;
            // tof cuts
            //      if(!b_K_tof) b_K = false;

            //pick the most likely particle from the TPC nsigmas
            if( b_PI
               && (d_TPCnSigmaPion < d_TPCnSigmaProton)
               && (d_TPCnSigmaPion < d_TPCnSigmaKaon) )
              { b_PI = true; b_PRO = false;}// b_K = false;}

            if( b_PRO
               && (d_TPCnSigmaProton < d_TPCnSigmaPion)
               && (d_TPCnSigmaProton < d_TPCnSigmaKaon) )
              { b_PRO = true; b_PI = false;}// b_K = false;}

            if( b_K
               && (d_TPCnSigmaKaon < d_TPCnSigmaProton)
               && (d_TPCnSigmaKaon < d_TPCnSigmaPion) )
              { b_K = true; b_PRO = false; b_PI = false;}


            // if(d_pidProbPion   < 0.8) b_PI  = false;
            // if(d_pidProbProton < 0.8) b_PRO = false;
            // if(d_pidProbKaon   < 0.8) b_K   = false;

            //      StBTofPidTraits  trk_tofPidTraits = gbltrk_i->btofPidTraits();


            double d_charge  = picoTrack->charge();
            //Anti lambda's too? ///////      if(d_charge <= 0) b_PRO = false;

            // cout<<" d_TPCnSigmaPion: "<<fabs(d_TPCnSigmaPion)
            //       <<" d_TPCnSigmaProton: "<<fabs(d_TPCnSigmaProton)
            //       <<" d_TPCnSigmaKaon: "<<fabs(d_TPCnSigmaKaon)
            //       <<" NOR : "<<(!(b_PI||b_PRO||b_K))<<endl;

            if(!(b_PI||b_PRO||b_K)) continue;

            if(!b_K) continue;

            // if( (b_PI&&b_PRO)
            //       ||(b_K&&b_PRO)
            //       ||(b_PI&&b_K) ) continue;

            double d_px0      = picoTrack->gMom().x();
            double d_py0      = picoTrack->gMom().y();
            double d_pz0      = picoTrack->gMom().z();
            double d_pT0      = picoTrack->gPt();
            double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

            // if( (d_pT0<d_pT_min) || (d_pT0 > d_pT_max)) continue;
            // if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max)) continue;

            if( d_pT0<d_pT_min) continue;
            if( d_mom0<d_mom_min) continue;






            // }//      if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max))

            //============== ===  Track QA Cuts  ==============
            bool b_bad_dEdx     = false;
            bool b_bad_tracking = false;

            //TOF Cuts
            bool b_bad_TOF_match = false;
            bool b_bad_ToF       = false;

            b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
            b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);

            // //TOF Cuts
            // ruins eta distro?
            //    b_bad_TOF_match = (gbltrk_i->btofPidTraits().matchFlag() < 1);
            //      b_bad_ToF       = (gbltrk_i->btofPidTraits().timeOfFlight() <= 0.0);
            bool b_bad_track     = b_bad_dEdx || b_bad_tracking || b_bad_TOF_match || b_bad_ToF;
            // trackinfo.b_bad_TOF_match = (track->btofPidTraits().matchFlag() < 1);
            // trackinfo.b_bad_ToF       = (track->btofPidTraits().timeOfFlight() <= 0.0);

            if(b_bad_track) continue;
            //============== END  Track QA Cuts  ==============

            StPicoPhysicalHelix trackhelix = picoTrack->helix(B);
            double helixpathl = trackhelix.pathLength(v3D_vtx, false);
            TVector3 v3D_dca = trackhelix.at(helixpathl)-v3D_vtx;
            double d_helix_DCA_r = v3D_dca.Mag();
            double d_DCA_r_cut = 1.0;

            if(b_K) d_DCA_r_cut = 1.0;

            if(d_helix_DCA_r > d_DCA_r_cut) continue;

            if(d_charge > 0.0) v_tracks_pl.push_back(picoTrack);
            else if(d_charge < 0.0) v_tracks_mi.push_back(picoTrack);
            v_tracks.push_back(picoTrack);


            trackinfo.reset();
            if(d_charge > 0.0)      trackinfo.b_pos_charge = true;
            else if(d_charge < 0.0) trackinfo.b_pos_charge = false;
            else continue;

            if(b_PI)  trackinfo.b_PI = true;
            if(b_PRO) trackinfo.b_PRO = true;
            if(b_K)   trackinfo.b_K = true;

            trackinfo.runNumber   = runNumber;
            trackinfo.eventNumber = eventNumber;

            trackinfo.MagField = d_MagField;

            trackinfo.px = d_px0;
            trackinfo.py = d_py0;
            trackinfo.pz = d_pz0;


            trackinfo.x_vtx = v3D_vtx.x();
            trackinfo.y_vtx = v3D_vtx.y();
            trackinfo.z_vtx = v3D_vtx.z();

            trackinfo.trackhelix = trackhelix;
            t_tracks -> Fill();
            index++;
          }

        //Fri Jul 28 20:42:22 EDT 2017
        //testing purposes, repeat for primary tracks
        vector<StPicoTrack *> v_pri_tracks;
        vector<StPicoTrack *> v_pri_tracks_pl;
        vector<StPicoTrack *> v_pri_tracks_mi;

        index = 0;
        double d_PI_m    = 0.13957018;
        double d_PRO_m   = 0.9382720813;
        double d_K_m     = 0.493677;
        double d_E_m     = 0.0005109989;
        // double d_PI_m2   = 0.019479835145;
        // double d_PRO_m2  = 0.880354498547;
        // double d_K_m2    = 0.243716980329;

        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

            // Retrieve i-th pico track
            StPicoTrack *picoTrack = dst->track(iTrk);

            if(picoTrack == NULL) continue;

            if(!picoTrack->isPrimary()) continue;
            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );

            unsigned short nHits = picoTrack->nHits();
            if(nHits < i_Nhits_min) continue;
            double d_TPCnSigmaElectron = fabs(picoTrack->nSigmaElectron());
            double d_TPCnSigmaPion     = fabs(picoTrack->nSigmaPion());
            double d_TPCnSigmaProton   = fabs(picoTrack->nSigmaProton());
            double d_TPCnSigmaKaon     = fabs(picoTrack->nSigmaKaon());



            double tofBeta           = -999;
            if(trait) tofBeta = trait->btofBeta();
            double d_tofBeta0 = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();


            // PID Cuts
            bool b_PI  = fabs(d_TPCnSigmaPion) < d_SigmaCutLevel;
            bool b_PRO = fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
            bool b_K   = fabs(d_TPCnSigmaKaon) < d_SigmaCutLevel;
            bool b_E   = fabs(d_TPCnSigmaElectron)< d_SigmaCutLevel;
            //      bool b_K   = fabs(d_TPCnSigmaKaon) < 0.5;


            //pick the most likely particle from the TPC nsigmas
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



            if(d_charge > 0.0)
              {
                if( b_E
                  && (fabs(d_TPCnSigmaElectron) < 3.0)
                  && (mass2 < 0.005)
                  && (mass2 > -0.008)
                  && (d_pT0 < 0.3)){
                  h_E_plus_pT -> Fill(d_pT0);
                  h_E_plus_y  -> Fill(d_y_K);
                  h2_E_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
                  h_E_plus_mT_Diff -> Fill((d_mT_E - d_E_m));
                }
                if( b_K
                  && (fabs(d_TPCnSigmaKaon) < 3.0)
                  && (mass2 > 0.15)
                  && (mass2 < 0.35)){
                  h_K_plus_pT -> Fill(d_pT0);
                  h_K_plus_y  -> Fill(d_y_K);
                  h2_K_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
                  h_K_plus_mT_Diff -> Fill((d_mT_K - d_K_m));
                }
                if( b_PI
                  && (fabs(d_TPCnSigmaPion) < 3.0)
                  && (mass2 > -0.03)
                  && (mass2 < 0.06)
                  && ((d_pT0 > 0.4) || (mass2 > 0.008))
                ){
                  h_PI_plus_pT -> Fill(d_pT0);
                  h_PI_plus_y  -> Fill(d_y_PI);
                  h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                  h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
                }
                if( b_PRO
                  && (fabs(d_TPCnSigmaProton) < 3.0)
                  && (mass2 > 0.7)
                  && (mass2 < 1.1))
                  {
                  h_PRO_plus_pT -> Fill(d_pT0);
                  h_PRO_plus_y  -> Fill(d_y_PRO);
                  h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
                  h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
                }
              }//if(d_charge > 0.0)
            else
              {
                if( b_E
                  && (fabs(d_TPCnSigmaElectron) < 3.0)
                  && (mass2 < 0.005)
                  && (mass2 > -0.008)
                  && (d_pT0 < 0.3)){
                  h_E_minus_pT -> Fill(d_pT0);
                  h_E_minus_y  -> Fill(d_y_K);
                  h2_E_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
                  h_E_minus_mT_Diff -> Fill((d_mT_E - d_E_m));
                }
                if( b_K
                  && (fabs(d_TPCnSigmaKaon) < 3.0)
                  && (mass2 > 0.15)
                  && (mass2 < 0.35)){
                  h_K_minus_pT -> Fill(d_pT0);
                  h_K_minus_y  -> Fill(d_y_K);
                  h2_K_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
                  h_K_minus_mT_Diff -> Fill((d_mT_K - d_K_m));
                }
                if( b_PI
                  && (fabs(d_TPCnSigmaPion) < 3.0)
                  && (mass2 > -0.03)
                  && (mass2 < 0.06)
                  && ((d_pT0 > 0.4) || (mass2 > 0.008))
                ){
                  h_PI_minus_pT -> Fill(d_pT0);
                  h_PI_minus_y  -> Fill(d_y_PI);
                  h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                  h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
                }
                if( b_PRO
                  && (fabs(d_TPCnSigmaProton) < 3.0)
                  && (mass2 > 0.7)
                  && (mass2 < 1.1))
                  {
                  h_PRO_minus_pT -> Fill(d_pT0);
                  h_PRO_minus_y  -> Fill(d_y_PRO);
                  h2_PRO_minus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
                  h_PRO_minus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
                }
              }//else (d_charge > 0.0)


            if(d_charge > 0.0)
              {
                if(b_PI||b_PRO) continue;
              }//if(d_charge > 0.0)
            else
              {
                if(!(b_PI||b_PRO||b_K)) continue;
                if(!b_K) continue;
              }//else (d_charge > 0.0)


            // if( (b_PI&&b_PRO)
            //       ||(b_K&&b_PRO)
            //       ||(b_PI&&b_K) ) continue;

            if( (d_pT0<d_pT_min) || (d_pT0 > d_pT_max)) continue;
            if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max)) continue;

            Kaoninfo.reset();

            // if( (d_mom0>d_mom_min) && (d_mom0 < 1.0))
            //     {
            // }//      if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max))

            //============== ===  Track QA Cuts  ==============
            bool b_bad_dEdx     = false;
            bool b_bad_tracking = false;

            //TOF Cuts
            bool b_bad_TOF_match = false;
            bool b_bad_ToF       = false;

            b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
            b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);

            // //TOF Cuts
            // ruins eta distro?
            if(trait) b_bad_TOF_match = (trait->btofMatchFlag() < 1);
            if(trait) b_bad_ToF       = (trait->btof() <= 0.0);
            bool b_bad_track     = b_bad_dEdx || b_bad_tracking; //|| b_bad_TOF_match || b_bad_ToF;
                                                                 // trackinfo.b_bad_TOF_match = (track->btofPidTraits().matchFlag() < 1);
                                                                 // trackinfo.b_bad_ToF       = (track->btofPidTraits().timeOfFlight() <= 0.0);

            if(b_bad_track) continue;
            //============== END  Track QA Cuts  ==============
            StPicoPhysicalHelix trackhelix = picoTrack->helix(B);
            double helixpathl = trackhelix.pathLength(v3D_vtx, false);
            TVector3 v3D_dca = trackhelix.at(helixpathl)-v3D_vtx;
            double d_helix_DCA_r = v3D_dca.Mag();
            double d_DCA_r_cut = 4.0;

            TVector3 v3D_obj_DCA = picoTrack->gDCA(pVtx);
            double d_obj_DCA = v3D_obj_DCA.Mag();


            if(d_helix_DCA_r > d_DCA_r_cut) continue;

            if(d_charge > 0.0) v_pri_tracks_pl.push_back(picoTrack);
            else if(d_charge < 0.0) v_pri_tracks_mi.push_back(picoTrack);
            v_tracks.push_back(picoTrack);

            //fill ttree
            if(d_charge < 0)      Kaoninfo.b_pos_charge = false;
            else                  Kaoninfo.b_pos_charge = true;
            Kaoninfo.runNumber   = runNumber;
            Kaoninfo.eventNumber = eventNumber;
            Kaoninfo.nGoodTracks = nGoodTracks;

            Kaoninfo.px = d_px0;
            Kaoninfo.py = d_py0;
            Kaoninfo.pz = d_pz0;

            Kaoninfo.tofBeta = tofBeta;
            //Kaoninfo.dEdxFit = dEdxFit;
            //Kaoninfo.d_dEdxTru = d_dEdxTru;

            Kaoninfo.x_vtx = v3D_vtx.x();
            Kaoninfo.y_vtx = v3D_vtx.y();
            Kaoninfo.z_vtx = v3D_vtx.z();

            Kaoninfo.d_TPCnSigmaKaon = d_TPCnSigmaKaon;
            //Kaoninfo.d_TOFnSigmaKaon = d_TOFnSigmaKaon;

            Kaoninfo.DCA_r = d_obj_DCA;

            Kaoninfo.b_bad_TOF = b_bad_ToF || b_bad_TOF_match;
            t_K -> Fill();
            index++;


          }
           // END testing purposes, repeat for primary tracks
           //  END Preliminary Track Loop
           //=========================================================

        i_event++;

    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    outputFile->cd();

    bool b_write_tree = false;
    if(b_write_tree)  t_tracks   -> Write();
    t_K ->Write();
    h2_m2_QA_pq->Write();
    h2_m2_QA_pT->Write();

    h_E_plus_mT_Diff-> Write();
    h_K_plus_mT_Diff-> Write();
    h_PI_plus_mT_Diff-> Write();
    h_PRO_plus_mT_Diff-> Write();

    h_E_minus_mT_Diff-> Write();
    h_K_minus_mT_Diff-> Write();
    h_PI_minus_mT_Diff-> Write();
    h_PRO_minus_mT_Diff-> Write();

    h_E_minus_pT   -> Write();
    h_K_minus_pT   -> Write();
    h_PI_minus_pT  -> Write();
    h_PRO_minus_pT -> Write();

    h_E_plus_pT   -> Write();
    h_K_plus_pT   -> Write();
    h_PI_plus_pT  -> Write();
    h_PRO_plus_pT -> Write();

    h_E_plus_y -> Write();
    h_K_plus_y -> Write();
    h_PI_plus_y -> Write();
    h_PRO_plus_y -> Write();

    h_E_minus_y -> Write();
    h_K_minus_y -> Write();
    h_PI_minus_y -> Write();
    h_PRO_minus_y -> Write();

    h2_E_plus_pT_vs_y -> Write();
    h2_K_plus_pT_vs_y -> Write();
    h2_PI_plus_pT_vs_y -> Write();
    h2_PRO_plus_pT_vs_y -> Write();

    h2_E_minus_pT_vs_y -> Write();
    h2_K_minus_pT_vs_y -> Write();
    h2_PI_minus_pT_vs_y -> Write();
    h2_PRO_minus_pT_vs_y -> Write();


    gROOT->GetListOfFiles()->Remove(outputFile);
    outputFile->Close();



  picoReader->Finish();

    //outputFile->Write();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    << std::endl;
}
