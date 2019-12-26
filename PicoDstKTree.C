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

    /* float pT; // is (*charge) */
    /* float pz; */
    //    float tofBeta;
    /* float dEdxFit; */
    //    float d_dEdxTru;
    /* float eta; */
    /* float phi; */
    /* float DCA_r; */
    /* float d_zvtx; */

    /* unsigned int nhits0; */

    /* float TPCnSigmaPion0; */
    /* float TPCnSigmaProton0; */
    /* float TPCnSigmaKaon0; */

    /* float d_pidProbPion0; */
    /* float d_pidProbProton0; */
    /* float d_pidProbKaon0; */

    /* float d_TOFnSigmaPion0; */
    /* float d_TOFnSigmaProton0; */
    /* float d_TOFnSigmaKaon0; */

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


      /* pair_mom = 0.0; */
      /* pair_pT = 0.0; */

      /* nhits0 = 0; */

      /* TPCnSigmaPion0 = 0.0; */
      /* TPCnSigmaProton0 = 0.0; */
      /* TPCnSigmaKaon0 = 0.0; */

      /* d_pidProbPion0 = 0.0; */
      /* d_pidProbProton0 = 0.0; */
      /* d_pidProbKaon0 = 0.0; */

      /* d_TOFnSigmaPion0 = 0.0; */
      /* d_TOFnSigmaProton0 = 0.0; */
      /* d_TOFnSigmaKaon0 = 0.0; */

      //      trackhelix0 = 0.0;

      /* nhits1 = 0; */

      /* TPCnSigmaPion1 = 0.0; */
      /* TPCnSigmaProton1 = 0.0; */
      /* TPCnSigmaKaon1 = 0.0; */

      /* d_pidProbPion1 = 0.0; */
      /* d_pidProbProton1 = 0.0; */
      /* d_pidProbKaon1 = 0.0; */

      /* d_TOFnSigmaPion1 = 0.0; */
      /* d_TOFnSigmaProton1 = 0.0; */
      /* d_TOFnSigmaKaon1 = 0.0; */

      //      trackhelix1 = 0.0;

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
void PicoDstKTree(const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
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

    int ievtcut[6] = {0};
    int itrkcut[6] = {0};

    double d_PI_m2   = 0.019479835;
    double d_PRO_m2  = 0.880354;
    double d_K_m2    = 0.24371698;

    TH1D *  h_evt_vs_cuts = new TH1D("h_evt_vs_cuts","h_evt_vs_cuts",11,-0.5,10.5);
    TH1D *  h_trk_vs_cuts = new TH1D("h_trk_vs_cuts","h_trk_vs_cuts",11,-0.5,10.5);

    TH1D *  h_evt      = new TH1D("h_evt","h_evt",1,0,1);
    TH1D *  h_vtx      = new TH1D("h_vtx","h_vtx",100,200,220);
    TH1D *  h_pT       = new TH1D("h_pT","h_pT",64,0.0,32.0);
    TH2D *  h2_dEdx_PI_pq = new TH2D("h2_dEdx_PI_pq","h2_dEdx_PI_pq",500,-2.0,2.0,500,0.6,3);
    TH2D *  h2_dEdx_PRO_pq = new TH2D("h2_dEdx_PRO_pq","h2_dEdx_PRO_pq",500,-2.0,2.0,500,0.6,3);
    TH2D *  h2_dEdx_K_pq = new TH2D("h2_dEdx_K_pq","h2_dEdx_K_pq",500,-2.0,2.0,500,0.6,3);
    TH2D *  h2_dEdx_All_pq = new TH2D("h2_dEdx_All_pq","h2_dEdx_All_pq",500,-2.0,2.0,2000,0.,10.);
    TH2D *  h2_dEdx_All_pq_1 = new TH2D("h2_dEdx_All_pq_1","h2_dEdx_All_pq_1",500,-2.0,2.0,2000,0.,10.);
    TH2D *  h2_m2_QA_pq = new TH2D("h2_m2_QA_pq","h2_m2_QA_pq",5000,-5.0,5.0,5000,-0.2,1.6);
    TH2D *  h2_m2_QA_pq_1 = new TH2D("h2_m2_QA_pq_1","h2_m2_QA_pq_1",5000,-5.0,5.0,5000,-0.2,1.6);
    TH2D *  h2_m2_QA_pT = new TH2D("h2_m2_QA_pT","h2_m2_QA_pT",5000,-5.0,5.0,5000,-2.0,1.6);
    TH2D *  h2_m2_QA_pT_1 = new TH2D("h2_m2_QA_pT_1","h2_m2_QA_pT_1",5000,-5.0,5.0,5000,-2.0,1.6);

    TH1D *  h_mult     = new TH1D("h_mult","h_mult",1600,0.0,1600.0);
    TH2D *  h2_mult     = new TH2D("h2_mult","h2_mult",400,0.0,400.0,100000,0.,100000.);//test

    TH1D *  h_DCA_r    = new TH1D("h_DCA_r","h_DCA_r",100,0.0,2.0);
    TH1D *  h_DCA_PRO    = new TH1D("h_DCA_PRO","h_DCA_PRO",100,0.0,2.0);
    TH1D *  h_DCA_PIM    = new TH1D("h_DCA_PIM","h_DCA_PIM",100,0.0,2.0);
    TH1D *  h_DCA_mother_r    = new TH1D("h_DCA_mother_r","h_DCA_mother_r",100,0.0,2.0);
    TH1D *  h_decay_length    = new TH1D("h_decay_length","h_decay_length",1000,0.0,20.0);

    TH1D *  h_K_DCA_r      = new TH1D("h_K_DCA_r","h_K_DCA_r",200,0.0,3.0);
    TH1D *  h_K_obj_DCA_r  = new TH1D("h_K_obj_DCA_r","h_K_obj_DCA_r",200,0.0,3.0);
    TH1D *  h_K_diff_DCA_r = new TH1D("h_K_diff_DCA_r","h_K_diff_DCA_r",400,-3.0,3.0);


    TH1D *  h_rapidity_pm = new TH1D("h_rapidity_pm","h_rapidity_pm",100,-2.0,2.0);
    TH2D *  h2_mT_rapidity_pm = new TH2D("h2_mT_rapidity_pm","h2_mT_rapidity_pm",100,0.0,1.0,20,-2.05,-0.05);

    TH1D *  h_nodecaylength_cut_inv_m_PHI    = new TH1D("h_nodecaylength_cut_inv_m_PHI","h_nodecaylength_cut_inv_m_PHI",1000,0.9,1.1);//.,2.);
    TH1D *  h_nodipangle_cut_inv_m_PHI    = new TH1D("h_nodipangle_cut_inv_m_PHI","h_nodipangle_cut_inv_m_PHI",1000,0.9,1.1);//.,2.);

    TH1D *  h_inv_m_PHI         = new TH1D("h_inv_m_PHI","h_inv_m_PHI",1000,0.9,1.1);
    TH1D *  h_prim_inv_m_PHI    = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",1000,0.9,1.1);
    //TH1D *  h_prim_inv_m_PHI    = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",1000,0.,2.1);


    TH1D *  h_inv_m_K0S         = new TH1D("h_inv_m_K0S","h_inv_m_K0S",1000,0.0,1.0);
    TH1D *  h_inv_m_RHO         = new TH1D("h_inv_m_RHO","h_inv_m_RHO",2000,0.45,1.2);
    //TH1D *  h_inv_m_LAMBDA = new TH1D("h_inv_m_LAMBDA","h_inv_m_LAMBDA",10000,0.0,10.0);
    TH1D *  h_inv_m_LAMBDA = new TH1D("h_inv_m_LAMBDA","h_inv_m_LAMBDA",1000,1.0,1.18);

    TH2D *  h2_pidprob_vs_nsigma_PI  = new TH2D("h2_pidprob_vs_nsigma_PI","h2_pidprob_vs_nsigma_PI",100,0.0,1.0,100,0.0,10.0);
    TH2D *  h2_pidprob_vs_nsigma_PRO = new TH2D("h2_pidprob_vs_nsigma_PRO","h2_pidprob_vs_nsigma_PRO",100,0.0,1.0,100,0.0,10.0);
    TH2D *  h2_pidprob_vs_nsigma_K   = new TH2D("h2_pidprob_vs_nsigma_K","h2_pidprob_vs_nsigma_K",100,0.0,1.0,100,0.0,10.0);

    TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_PI  = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_PI", "h2_TOF_nsigma_vs_TPC_nsigma_PI",100,0.0,10.0,100,0.0,10.0);
    TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_PRO = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_PRO","h2_TOF_nsigma_vs_TPC_nsigma_PRO",100,0.0,10.0,100,0.0,10.0);
    TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_K   = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_K",  "h2_TOF_nsigma_vs_TPC_nsigma_K",100,0.0,10.0,100,0.0,10.0);


    TH2D *  h2_TOF_nsigma_K_vs_PI  = new TH2D("h2_TOF_nsigma_K_vs_PI",   "h2_TOF_nsigma_K_vs_PI",100,0.0,10.0,100,0.0,10.0);
    TH2D *  h2_TPC_nsigma_K_vs_PI  = new TH2D("h2_TPC_nsigma_K_vs_PI",   "h2_TPC_nsigma_K_vs_PI",100,0.0,10.0,100,0.0,10.0);
    TH2D *  h2_TOF_nsigma_K_vs_PRO = new TH2D("h2_TOF_nsigma_K_vs_PRO",  "h2_TOF_nsigma_K_vs_PRO",100,0.0,10.0,100,0.0,10.0);
    TH2D *  h2_TPC_nsigma_K_vs_PRO = new TH2D("h2_TPC_nsigma_K_vs_PRO",  "h2_TPC_nsigma_K_vs_PRO",100,0.0,10.0,100,0.0,10.0);

    TH1D *  h_PHI_decay_length     = new TH1D("h_PHI_decay_length","h_PHI_decay_length",100,0,10.0);
    TH1D *  h_dip_angle            = new TH1D("h_dip_angle","h_dip_angle",1000,-4,4);

    TH1D *  h_K_plus_pT   = new TH1D("h_K_plus_pT","Kaon+ pT",1000,0.0,5.0);
    TH1D *  h_PI_plus_pT  = new TH1D("h_PI_plus_pT","Pion+ pT",1000,0.0,5.0);
    TH1D *  h_PRO_plus_pT = new TH1D("h_PRO_plus_pT","Proton+ pT",1000,0.0,5.0);

    TH1D *  h_K_minus_pT   = new TH1D("h_K_minus_pT","Kaon- pT",1000,0.0,5.0);
    TH1D *  h_PI_minus_pT  = new TH1D("h_PI_minus_pT","Pion- pT",1000,0.0,5.0);
    TH1D *  h_PRO_minus_pT = new TH1D("h_PRO_minus_pT","Proton- pT",1000,0.0,5.0);

    TH1D *  h_K_plus_y   = new TH1D("h_K_plus_y","Kaon+ y",1000,-5.0,0.0);
    TH1D *  h_PI_plus_y  = new TH1D("h_PI_plus_y","Pion+ y",1000,-5.0,0.0);
    TH1D *  h_PRO_plus_y = new TH1D("h_PRO_plus_y","Proton+ y",1000,-5.0,0.0);

    TH1D *  h_K_minus_y   = new TH1D("h_K_minus_y","Kaon- y",1000,-5.0,0.0);
    TH1D *  h_PI_minus_y  = new TH1D("h_PI_minus_y","Pion- y",1000,-5.0,0.0);
    TH1D *  h_PRO_minus_y = new TH1D("h_PRO_minus_y","Proton- y",1000,-5.0,0.0);

    TH2D *  h2_K_plus_pT_vs_y   = new TH2D("h2_K_plus_pT_vs_y","h2_K_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PI_plus_pT_vs_y  = new TH2D("h2_PI_plus_pT_vs_y","h2_PI_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);
    TH2D *  h2_PRO_plus_pT_vs_y = new TH2D("h2_PRO_plus_pT_vs_y","h2_PRO_plus_pT_vs_y",1000,-5.0,0.0,1000,0.0,5.0);

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
        h_vtx -> Fill(d_zvtx);
        bool b_bad_zvtx        = ((d_zvtx < 210.0) || (d_zvtx > 212.0));
        //-------- END Z-VTX Selection ----------------------------

        // Track analysis
        Int_t nTracks = dst->numberOfTracks();
        Int_t nMatrices = dst->numberOfTrackCovMatrices();
        if(nTracks != nMatrices) {
            //std::cout << "Number of tracks and matrices do not match!" << std::endl;
        }
        //std::cout << "Number of tracks in event: " << nTracks << std::endl;

        bool b_cent_01  = false; // 0  < centrality <= 5%
        bool b_cent_02  = false; // 5  < centrality <= 10%
        bool b_cent_03  = false; // 10 < centrality <= 15%
        bool b_cent_04  = false; // 15 < centrality <= 20%
        bool b_cent_05  = false; // 20 < centrality <= 25%
        bool b_cent_06  = false; // 25 < centrality <= 30%
        bool b_cent_07  = false; // centrality > 30%

        // # of cuts = 0
        ievtcut[0] = ievtcut[0] + 1;
        itrkcut[0] = itrkcut[0] + nTracks;
        // std::cout << "Event(s) with 0 cut(s)" << ievtcut[0] << std::endl;
        // std::cout << "Track(s) with 0 cut(s)" << itrkcut[0] << std::endl;

        // Track loop
        int nGoodTracks = 0;
        int nTrkvsCuts  = 0;
        h2_dEdx_All_pq_1->Reset();
        h2_m2_QA_pq_1   ->Reset();
        h2_m2_QA_pT_1   ->Reset();
        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

            // Retrieve i-th pico track
            StPicoTrack *picoTrack = dst->track(iTrk);

            if(!picoTrack) continue;

            // # of cuts = 1
            // h_evt_vs_cuts -> Fill(1);
            itrkcut[1] = itrkcut[1] + 1;

            //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;


            //StPicoPhysicalHelix helix = picoTrack->helix(B);
            //cout<<"Testing helix parameter: "<<helix.period()<<endl;

            //============== === Track Level Cuts ==============
            bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
            bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
            //bool b_not_enough_hits = false;//((double)picoTrack->nHitsFit()) < 15;
            bool b_bad_track    = b_bad_dEdx || b_bad_tracking /*|| b_not_enough_hits*/;
            if(b_bad_track) continue;

            // # of cuts = 2
            // h_evt_vs_cuts -> Fill(2);
            itrkcut[2] = itrkcut[2] + 1;

            //============== END Track Level Cuts ==============
            //============= Primary Track Cut ===================
            if(!picoTrack->isPrimary()) continue;

            // # of cuts = 3
            // h_evt_vs_cuts -> Fill(3);
            itrkcut[3] = itrkcut[3] + 1;

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

            h2_dEdx_All_pq_1->Fill(d_mom0/(picoTrack->charge()),picoTrack->dEdx());
            nGoodTracks++;
            if(d_tofBeta0 == -999) continue;
            // # of cuts = 4
            // h_evt_vs_cuts -> Fill(4);
            itrkcut[4] = itrkcut[4] + 1;

            h2_m2_QA_pq_1   ->Fill(d_mom0/(picoTrack->charge()),mass2);
            h2_m2_QA_pT_1   ->Fill(d_pT0/(picoTrack->charge()),mass2);
            nTrkvsCuts++;

            /*
            // Check if track has TOF signal
            if( picoTrack->isTofTrack() ) {
                // Retrieve corresponding trait
                StPicoBTofPidTraits *trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
                if( !trait ) {
                    std::cout << "O-oh... No BTofPidTrait # " << picoTrack->bTofPidTraitsIndex()
                    << " for track # " << iTrk << std::endl;
                    std::cout << "Check that you turned on the branch!" << std::endl;
                    continue;
                }
                // Fill beta
                //hTofBeta->Fill( trait->btofBeta() );
            } //if( isTofTrack() )
            */

        } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

        // h_mult->Fill(nGoodTracks);

        ievtcut[1] = ievtcut[1]+1;
        ievtcut[2] = ievtcut[2]+1;
        ievtcut[3] = ievtcut[3]+1;
        ievtcut[4] = ievtcut[4]+1;


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

        bool b_bad_evt  = /*b_0vtx ||*/ /*b_low_mult || b_pileup ||*/ b_bad_zvtx /*|| (!b_good_fxt_triger)*/ || b_bad_trig /*|| b_bad_vtxID*/;
        if(b_bad_evt) continue; // Events skipped are not counted in the event count

        // # of cuts = 5
        ievtcut[5] = ievtcut[5] + 1;
        itrkcut[5] = itrkcut[5] + nTrkvsCuts;

        h_mult->Fill(nGoodTracks);

        h2_dEdx_All_pq->Add(h2_dEdx_All_pq_1);
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



            // if( (d_mom0>d_mom_min) && (d_mom0 < 1.0))
            //     {
            //if(b_PI)  h2_pidprob_vs_nsigma_PI -> Fill(d_pidProbPion,d_TPCnSigmaPion);
            //if(b_PRO) h2_pidprob_vs_nsigma_PRO -> Fill(d_pidProbProton,d_TPCnSigmaProton);
            //if(b_K)   h2_pidprob_vs_nsigma_K -> Fill(d_pidProbKaon,d_TPCnSigmaKaon);

            // if(b_PI)  h2_TOF_nsigma_vs_TPC_nsigma_PI  -> Fill(d_TOFnSigmaPion,d_TPCnSigmaPion);
            // if(b_PRO) h2_TOF_nsigma_vs_TPC_nsigma_PRO -> Fill(d_TOFnSigmaProton,d_TPCnSigmaProton);
            // if(b_K)   h2_TOF_nsigma_vs_TPC_nsigma_K   -> Fill(d_TOFnSigmaKaon,d_TPCnSigmaKaon);

            //h2_TOF_nsigma_vs_TPC_nsigma_PI  -> Fill(d_TOFnSigmaPion,d_TPCnSigmaPion);
            //h2_TOF_nsigma_vs_TPC_nsigma_PRO -> Fill(d_TOFnSigmaProton,d_TPCnSigmaProton);
            //h2_TOF_nsigma_vs_TPC_nsigma_K   -> Fill(d_TOFnSigmaKaon,d_TPCnSigmaKaon);

            //h2_TOF_nsigma_K_vs_PI  -> Fill(d_TOFnSigmaKaon,d_TOFnSigmaPion);
            h2_TPC_nsigma_K_vs_PI  -> Fill(d_TPCnSigmaKaon,d_TPCnSigmaPion);
            //h2_TOF_nsigma_K_vs_PRO -> Fill(d_TOFnSigmaKaon,d_TOFnSigmaProton);
            h2_TPC_nsigma_K_vs_PRO -> Fill(d_TPCnSigmaKaon,d_TPCnSigmaProton);

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
            h_DCA_r -> Fill(d_helix_DCA_r);
            double d_DCA_r_cut = 1.0;
            if(b_PRO) {d_DCA_r_cut = d_PRO_daughter_DCA; h_DCA_PRO->Fill(d_helix_DCA_r);  }
            if(b_PI)  {d_DCA_r_cut = d_PI_daughter_DCA; h_DCA_PIM->Fill(d_helix_DCA_r);  }
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

            // trackinfo.tofBeta = tofBeta;
            // trackinfo.dEdxFit = dEdxFit;
            // trackinfo.d_dEdxTru = d_dEdxTru;


            trackinfo.x_vtx = v3D_vtx.x();
            trackinfo.y_vtx = v3D_vtx.y();
            trackinfo.z_vtx = v3D_vtx.z();

            trackinfo.trackhelix = trackhelix;
            t_tracks -> Fill();
            //      cout<<" track !"<<endl;
            index++;
          }//  while ( ( track = (StPicoTrack*)itr_globaltrk.Next() ) )

        //Fri Jul 28 20:42:22 EDT 2017
        //testing purposes, repeat for primary tracks
        vector<StPicoTrack *> v_pri_tracks;
        vector<StPicoTrack *> v_pri_tracks_pl;
        vector<StPicoTrack *> v_pri_tracks_mi;

        index = 0;
        double d_PI_m    = 0.13957018;
        double d_PRO_m   = 0.9382720813;
        double d_K_m     = 0.493677;
        // double d_PI_m2   = 0.019479835145;
        // double d_PRO_m2  = 0.880354498547;
        // double d_K_m2    = 0.243716980329;

        for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

            // Retrieve i-th pico track
            StPicoTrack *picoTrack = dst->track(iTrk);

            if(picoTrack == NULL) continue;
            //short flag = picoTrack->flag();
            //if(flag <=0 || flag >=1000 )continue;
            if(!picoTrack->isPrimary()) continue;
            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
            //bool tofAvailability = trait;
            unsigned short nHits = picoTrack->nHits();
            if(nHits < i_Nhits_min) continue;
            double d_TPCnSigmaPion   = fabs(picoTrack->nSigmaPion());
            double d_TPCnSigmaProton = fabs(picoTrack->nSigmaProton());
            double d_TPCnSigmaKaon   = fabs(picoTrack->nSigmaKaon());

            //double d_pidProbPion     = primtrk_i -> pidProbPion();
            //double d_pidProbProton   = primtrk_i -> pidProbProton();
            //double d_pidProbKaon     = primtrk_i -> pidProbKaon();

            //double d_TOFnSigmaPion   = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaPion());//Not recommended in PID
            //double d_TOFnSigmaProton = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaProton());//Not recommended in PID
            //double d_TOFnSigmaKaon   = fabs((dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->sigmaKaon());//Not recommended in PID

            double tofBeta           = -999;
            if(trait) tofBeta = trait->btofBeta();
            double d_tofBeta0 = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();
            //double dEdxFit           = (dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->dEdxFit();
            //double d_dEdxTru         = (dst->btofPidTraits(picoTrack->bTofPidTraitsIndex()))->dEdxTruncated();

            // PID Cuts
            bool b_PI  = fabs(d_TPCnSigmaPion) < d_SigmaCutLevel;
            bool b_PRO = fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
            bool b_K   = fabs(d_TPCnSigmaKaon) < d_SigmaCutLevel;
            //      bool b_K   = fabs(d_TPCnSigmaKaon) < 0.5;

            //bool b_K_tof   = fabs(d_TOFnSigmaKaon) < d_SigmaCutLevel;
            // tof cuts
            //      if(!b_K_tof) b_K = false;

            //pick the most likely particle from the TPC nsigmas
            if( b_PI
               && (d_TPCnSigmaPion < d_TPCnSigmaProton)
               && (d_TPCnSigmaPion < d_TPCnSigmaKaon) )
              { b_PI = true; b_PRO = false; b_K = false;}

            if( b_PRO
               && (d_TPCnSigmaProton < d_TPCnSigmaPion)
               && (d_TPCnSigmaProton < d_TPCnSigmaKaon) )
              { b_PRO = true; b_PI = false;b_K = false;}

            if( b_K
               && (d_TPCnSigmaKaon < d_TPCnSigmaProton)
               && (d_TPCnSigmaKaon < d_TPCnSigmaPion) )
              { b_K = true; b_PRO = false; b_PI = false;}


            // if(d_pidProbPion   < 0.8) b_PI  = false;
            // if(d_pidProbProton < 0.8) b_PRO = false;
            // if(d_pidProbKaon   < 0.8) b_K   = false;
            //      StBTofPidTraits  trk_tofPidTraits = primtrk_i->btofPidTraits();

            double d_charge  = picoTrack->charge();
            //Anti lambda's too? ///////      if(d_charge <= 0) b_PRO = false;

            // cout<<" d_TPCnSigmaPion: "<<fabs(d_TPCnSigmaPion)
            //       <<" d_TPCnSigmaProton: "<<fabs(d_TPCnSigmaProton)
            //       <<" d_TPCnSigmaKaon: "<<fabs(d_TPCnSigmaKaon)
            //       <<" NOR : "<<(!(b_PI||b_PRO||b_K))<<endl;

            /// relax cuts for k- --> k+ Fri Aug 10 15:01:45 EDT 2018

            double d_px0        = picoTrack->pMom().x();
            double d_py0        = picoTrack->pMom().y();
            double d_pz0        = picoTrack->pMom().z();
            double d_pT0        = picoTrack->pPt();
            double d_mom0       = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

            double mass2        = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

            double d_E_PRO      = sqrt(d_PRO_m*d_PRO_m + d_mom0*d_mom0);
            double d_E_PI       = sqrt(d_PI_m*d_PI_m + d_mom0*d_mom0);
            double d_E_K        = sqrt(d_K_m*d_K_m + d_mom0*d_mom0);

            double d_y_PRO      = 0.5*TMath::Log((d_E_PRO + d_pz0)/(d_E_PRO - d_pz0));
            double d_y_PI       = 0.5*TMath::Log((d_E_PI + d_pz0)/(d_E_PI - d_pz0));
            double d_y_K        = 0.5*TMath::Log((d_E_K + d_pz0)/(d_E_K - d_pz0));



            if(d_charge > 0.0)
              {
                if( b_K
                  && (fabs(d_TPCnSigmaKaon) < 3.0)
                  && (mass2 > 0.15)
                  && (mass2 < 0.35)){
                  h_K_plus_pT -> Fill(d_pT0);
                  h_K_plus_y  -> Fill(d_y_K);
                  h2_K_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
                }
                if( b_PI
                  && (fabs(d_TPCnSigmaPion) < 3.0)
                  && (mass2 > -0.03)
                  && (mass2 < 0.06)){
                  h_PI_plus_pT -> Fill(d_pT0);
                  h_PI_plus_y  -> Fill(d_y_PI);
                  h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                }
                if( b_PRO
                  && (fabs(d_TPCnSigmaProton) < 3.0)
                  && (mass2 > 0.7)
                  && (mass2 < 1.1))
                  {
                  h_PRO_plus_pT -> Fill(d_pT0);
                  h_PRO_plus_y  -> Fill(d_y_PRO);
                  h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
                }
              }//if(d_charge > 0.0)
            else
              {
                if( b_K
                  && (fabs(d_TPCnSigmaKaon) < 3.0)
                  && (mass2 > 0.15)
                  && (mass2 < 0.35)){
                  h_K_minus_pT -> Fill(d_pT0);
                  h_K_minus_y  -> Fill(d_y_K);
                  h2_K_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
                }
                if( b_PI
                  && (fabs(d_TPCnSigmaPion) < 3.0)
                  && (mass2 > -0.03)
                  && (mass2 < 0.06)){
                  h_PI_minus_pT -> Fill(d_pT0);
                  h_PI_minus_y  -> Fill(d_y_PI);
                  h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
                }
                if( b_PRO
                  && (fabs(d_TPCnSigmaProton) < 3.0)
                  && (mass2 > 0.7)
                  && (mass2 < 1.1))
                  {
                  h_PRO_minus_pT -> Fill(d_pT0);
                  h_PRO_minus_y  -> Fill(d_y_PRO);
                  h2_PRO_minus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
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

            h_K_DCA_r      -> Fill(d_helix_DCA_r);
            h_K_obj_DCA_r  -> Fill(d_obj_DCA);
            h_K_diff_DCA_r -> Fill(d_helix_DCA_r-d_obj_DCA);

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

            // if(primtrk_i == NULL) continue;
            // short flag = primtrk_i->flag();
            // if(flag <=0 || flag >=1000 )continue;
            // unsigned short nHits = primtrk_i->nHits();
            // if(nHits < i_Nhits_min) continue;
            // double d_TPCnSigmaPion   = fabs(primtrk_i->nSigmaPion());
            // double d_TPCnSigmaProton = fabs(primtrk_i->nSigmaProton());
            // double d_TPCnSigmaKaon   = fabs(primtrk_i->nSigmaKaon());

            // double d_pidProbPion     = primtrk_i -> pidProbPion();
            // double d_pidProbProton   = primtrk_i -> pidProbProton();
            // double d_pidProbKaon     = primtrk_i -> pidProbKaon();

            // double d_TOFnSigmaPion   = (primtrk_i->btofPidTraits()).sigmaPion();
            // double d_TOFnSigmaProton = (primtrk_i->btofPidTraits()).sigmaProton();
            // double d_TOFnSigmaKaon   = (primtrk_i->btofPidTraits()).sigmaKaon();

            // // PID Cuts
            // bool b_PI  = fabs(d_TPCnSigmaPion) < d_SigmaCutLevel;
            // bool b_PRO = fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
            // bool b_K   = fabs(d_TPCnSigmaKaon) < d_SigmaCutLevel;

            // bool b_K_tof   = fabs(d_TOFnSigmaKaon) < d_SigmaCutLevel;
            // // tof cuts
            // //      if(!b_K_tof) b_K = false;

            // //pick the most likely particle from the TPC nsigmas
            // if( b_PI
            //       && (d_TPCnSigmaPion < d_TPCnSigmaProton)
            //       && (d_TPCnSigmaPion < d_TPCnSigmaKaon) )
            //     { b_PI = true; b_PRO = false;}// b_K = false;}

            // if( b_PRO
            //       && (d_TPCnSigmaProton < d_TPCnSigmaPion)
            //       && (d_TPCnSigmaProton < d_TPCnSigmaKaon) )
            //     { b_PRO = true; b_PI = false;}// b_K = false;}

            // if( b_K
            //       && (d_TPCnSigmaKaon < d_TPCnSigmaProton)
            //       && (d_TPCnSigmaKaon < d_TPCnSigmaPion) )
            //     { b_K = true; b_PRO = false; b_PI = false;}

            // double d_charge  = primtrk_i->charge();

            // if(!(b_PI||b_PRO||b_K)) continue;

            // //just to Ks for now
            // if(!b_K) continue;

            // double d_px0      = primtrk_i->p().x();
            // double d_py0      = primtrk_i->p().y();
            // double d_pz0      = primtrk_i->p().z();
            // double d_pT0      = primtrk_i->pt();
            // double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

            // if( (d_pT0<d_pT_min) || (d_pT0 > d_pT_max)) continue;
            // if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max)) continue;

            // //============== ===  Track QA Cuts  ==============
            // bool b_bad_dEdx     = false;
            // bool b_bad_tracking = false;

            // //TOF Cuts
            // bool b_bad_TOF_match = false;
            // bool b_bad_ToF       = false;

            // b_bad_dEdx     = (primtrk_i->nHitsDedx() <= 0);
            // b_bad_tracking = (((double)primtrk_i->nHitsFit(kTpcId) / (double)primtrk_i->nHitsPoss(kTpcId)) < 0.52);

            // bool b_bad_track     = b_bad_dEdx || b_bad_tracking || b_bad_TOF_match || b_bad_ToF;
            // if(b_bad_track) continue;
            // //============== END  Track QA Cuts  ==============
            // StPhysicalHelixD trackhelix = primtrk_i->helix();
            // double helixpathl = trackhelix.pathLength(v3D_vtx, false);
            // TVector3 v3D_dca = trackhelix.at(helixpathl)-v3D_vtx;
            // double d_helix_DCA_r = v3D_dca.Mag();
            // double d_DCA_r_cut = 1.0;

            // if(d_helix_DCA_r < d_DCA_r_cut) continue;

            // if(d_charge > 0.0) v_pri_tracks_pl.push_back(primtrk_i);
            // else if(d_charge < 0.0) v_pri_tracks_mi.push_back(primtrk_i);
            // v_pri_tracks.push_back(primtrk_i);
          }//  while ( ( primtrk_i = (StPicoTrack*)itr_globaltrk.Next() ) )  //=========================================================
          // h_mult->Fill(index);
           // END testing purposes, repeat for primary tracks
           //  END Preliminary Track Loop
           //=========================================================
           // cout<<endl<<" v_tracks_pl.size(): "<< v_tracks_pl.size()<<" v_tracks_mi.size(): "<< v_tracks_mi.size()<<endl;
           // cout<<" v_pri_tracks_pl.size(): "<< v_pri_tracks_pl.size()<<" v_pri_tracks_mi.size(): "<< v_pri_tracks_mi.size()<<endl<<endl;

        //=========================================================
        //      Invariant Mass Nested Track Loop
        //=========================================================
        // double d_PI_m    = 0.13957018;
        // double d_PRO_m   = 0.9382720813;
        // double d_K_m     = 0.493677;
        if(v_tracks_pl.size() <= 0) continue;
        if(v_tracks_mi.size() <= 0) continue;
        vector<StPicoTrack *> v_tracks0 = v_tracks_pl;
        vector<StPicoTrack *> v_tracks1 = v_tracks_mi;

        for(int i = 0; i < v_tracks0.size();i++)
            //  for(int i = 0; i < v_tracks0.size()-1;i++)
          {
            StPicoTrack * picoTrack0 = v_tracks0[i];
            if(!picoTrack0) continue;

            double d_TPCnSigmaPion0   = fabs(picoTrack0->nSigmaPion());
            double d_TPCnSigmaProton0 = fabs(picoTrack0->nSigmaProton());
            double d_TPCnSigmaKaon0   = fabs(picoTrack0->nSigmaKaon());

            // PID Cuts
            bool b_PI0  = fabs(d_TPCnSigmaPion0) < d_SigmaCutLevel;
            bool b_PRO0 = fabs(d_TPCnSigmaProton0) < d_SigmaCutLevel;
            bool b_K0   = fabs(d_TPCnSigmaKaon0) < d_SigmaCutLevel;

            //pick the most likely particle from the TPC nsigmas
            if( b_PI0
               && (d_TPCnSigmaPion0 < d_TPCnSigmaProton0)
               && (d_TPCnSigmaPion0 < d_TPCnSigmaKaon0) )
              { b_PI0 = true; b_PRO0 = false;}// b_K0 = false;}

            if( b_PRO0
               && (d_TPCnSigmaProton0 < d_TPCnSigmaPion0)
               && (d_TPCnSigmaProton0 < d_TPCnSigmaKaon0) )
              { b_PRO0 = true; b_PI0 = false;}// b_K0 = false;}

            if( b_K0
               && (d_TPCnSigmaKaon0 < d_TPCnSigmaProton0)
               && (d_TPCnSigmaKaon0 < d_TPCnSigmaPion0) )
              { b_K0 = true; b_PRO0 = false; b_PI0 = false;}

            //double d_TOFnSigmaKaon0   = (v_tracks0[i]->btofPidTraits()).sigmaKaon();
            //bool b_K_tof0   = fabs(d_TOFnSigmaKaon0) < d_SigmaCutLevel;
            //tof cuts
            //      if(!b_K_tof0) b_K0 = false;

            double d_M0 = -9999.0;
            double d_px0      = picoTrack0->gMom().x();
            double d_py0      = picoTrack0->gMom().y();
            double d_pz0      = picoTrack0->gMom().z();
            double d_pT0      = picoTrack0->gPt();
            double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack0->isTofTrack()) trait = dst->btofPidTraits(picoTrack0->bTofPidTraitsIndex());
            double d_tofBeta0 = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();
            //double d_dEdxFit0 = picoTrack->probPidTraits().dEdxFit();
            //double d_dEdxTru0 = picoTrack->probPidTraits().dEdxTruncated();
            double d_charge0  = picoTrack0->charge();
            double d_mc0      = d_mom0/d_charge0;
            double d_eta0     = picoTrack0->gMom().Eta();
            double d_phi0     = picoTrack0->gMom().Phi();

            StPicoPhysicalHelix    trackhelix0 = picoTrack0->helix(B);

            if(b_PRO0)
              {
                d_M0 = d_PRO_m;
                if((d_tofBeta0 != 0.0) && (d_charge0 != 0.0)) h2_dEdx_PRO_pq -> Fill(d_mc0,1.0/d_tofBeta0);
              }//if(b_PRO0)
            else if(b_PI0)
              {
                d_M0 = d_PI_m;
                if((d_tofBeta0 != 0.0) && (d_charge0 != 0.0)) h2_dEdx_PI_pq -> Fill(d_mc0,1.0/d_tofBeta0);
              }//else if(b_PI0)
            else if(b_K0)
              {
                d_M0 = d_K_m;
                if((d_tofBeta0 != 0.0) && (d_charge0 != 0.0)) h2_dEdx_K_pq -> Fill(d_mc0,1.0/d_tofBeta0);
              }//else if(b_K0)

            for(int j = 0; j < v_tracks1.size(); j++)
                //      for(int j = i+1; j < v_tracks1.size(); j++)
              {
                StPicoTrack * picoTrack1 = v_tracks1[j];

                if(!picoTrack1 || (picoTrack0->id() == picoTrack1->id())) continue;
                double d_TPCnSigmaPion1   = fabs(picoTrack1->nSigmaPion());
                double d_TPCnSigmaProton1 = fabs(picoTrack1->nSigmaProton());
                double d_TPCnSigmaKaon1   = fabs(picoTrack1->nSigmaKaon());

                // PID Cuts
                bool b_PI1  = fabs(d_TPCnSigmaPion1) < d_SigmaCutLevel;
                bool b_PRO1 = fabs(d_TPCnSigmaProton1) < d_SigmaCutLevel;
                bool b_K1   = fabs(d_TPCnSigmaKaon1) < d_SigmaCutLevel;

                //pick the most likely particle from the TPC nsigmas
                if( b_PI1
                   && (d_TPCnSigmaPion1 < d_TPCnSigmaProton1)
                   && (d_TPCnSigmaPion1 < d_TPCnSigmaKaon1) )
                  { b_PI1 = true; b_PRO1 = false;}// b_K1 = false;}

                if( b_PRO1
                   && (d_TPCnSigmaProton1 < d_TPCnSigmaPion1)
                   && (d_TPCnSigmaProton1 < d_TPCnSigmaKaon1) )
                  { b_PRO1 = true; b_PI1 = false;}// b_K1 = false;}

                if( b_K1
                   && (d_TPCnSigmaKaon1 < d_TPCnSigmaProton1)
                   && (d_TPCnSigmaKaon1 < d_TPCnSigmaPion1) )
                  { b_K1 = true; b_PRO1 = false; b_PI1 = false;}

                //double d_TOFnSigmaKaon1   = (v_tracks1[j]->btofPidTraits()).sigmaKaon();
                //bool b_K_tof1   = fabs(d_TOFnSigmaKaon1) < d_SigmaCutLevel;
                //tof cuts
                // if(!b_K_tof1) b_K1 = false;

                double d_M1 = -9999.0;
                if(b_PRO1) d_M1 = d_PRO_m;
                else if(b_PI1)  d_M1 = d_PI_m;
                else if(b_K1)   d_M1 = d_K_m;

                double d_px1      = picoTrack1->gMom().x();
                double d_py1      = picoTrack1->gMom().y();
                double d_pz1      = picoTrack1->gMom().z();
                double d_pT1      = picoTrack1->gPt();
                double d_mom1     = sqrt(d_pT1*d_pT1 + d_pz1*d_pz1);
                StPicoBTofPidTraits *trait = NULL;
                if(picoTrack1->isTofTrack()) trait = dst->btofPidTraits(picoTrack1->bTofPidTraitsIndex());
                double d_tofBeta1 = -999;
                if(trait) d_tofBeta1 = trait->btofBeta();
                //double d_dEdxFit1 = picoTrack1->probPidTraits().dEdxFit();
                //double d_dEdxTru1 = picoTrack1->probPidTraits().dEdxTruncated();
                double d_charge1  = picoTrack1->charge();

                double d_mc1      = d_mom1/d_charge1;
                double d_eta1     = picoTrack1->gMom().Eta();
                double d_phi1     = picoTrack1->gMom().Phi();

                //final cuts
                if(d_charge0 == d_charge1) continue;
                bool b_PHI    = b_K0 && b_K1;
                bool b_RHO    = b_PI0 && b_PI1;
                bool b_K0S    = b_RHO;

                bool b_LAMBDA = (b_PRO0 && b_PI1)||(b_PRO1 && b_PI0);
                bool b_V0 = b_PHI || b_RHO || b_LAMBDA;
                if(!b_V0) continue;

                // computational intensive, do all cuts before
                StPicoPhysicalHelix    trackhelix1 = picoTrack1->helix(B);
                // Get Mom at decayvtx:
                pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);

                TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, d_MagField*kilogauss);
                TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, d_MagField*kilogauss);

                TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
                TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);

                double d_dca_products =(v3D_x_daughter0-v3D_x_daughter1).Mag();
                //      if(d_dca_products > d_cut_dca_daughters) continue;
                if(d_dca_products > d_cut_dca_daughters_lam) b_LAMBDA = false;
                if(d_dca_products > d_cut_dca_daughters_k0s) b_RHO    = false;

                TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)*0.5;
                TVector3 v3D_xvec_decayl = v3D_x_mother - v3D_vtx;
                TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;

                double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);

                //DCA btw mother and vtx
                double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );

                h_DCA_mother_r -> Fill(d_dca_mother);
                //if(d_dca_mother > d_cut_dca_mother) continue;
                if(d_dca_mother > d_cut_dca_mother_lam) b_LAMBDA = false;
                if(d_dca_mother > d_cut_dca_mother_k0s) b_RHO    = false;

                double d_mother_decay_length =  v3D_xvec_decayl.Mag();
                h_decay_length -> Fill(d_mother_decay_length);
                //if(d_mother_decay_length < d_cut_mother_decay_length) continue;
                if(d_mother_decay_length < d_cut_mother_decay_length_lam) b_LAMBDA = false;
                if(d_mother_decay_length < d_cut_mother_decay_length_k0s) b_K0S    = false;
                if(d_mother_decay_length > d_cut_mother_decay_length_RHO) b_RHO    = false;
                if(b_PHI) h_PHI_decay_length -> Fill(d_mother_decay_length);

                double d_E0 = sqrt(v3D_p_daughter0.Mag2()+d_M0*d_M0);
                double d_E1 = sqrt(v3D_p_daughter1.Mag2()+d_M1*d_M1);
                double d_inv_m = sqrt(d_M0*d_M0
                                      +d_M1*d_M1
                                      +2.0*d_E0*d_E1
                                      -2.0*(v3D_p_daughter0.Dot(v3D_p_daughter1)) );

                if(b_PHI) h_nodecaylength_cut_inv_m_PHI -> Fill(d_inv_m);
                if(d_mother_decay_length > d_cut_mother_decay_length_PHI) b_PHI    = false;
                if(b_PHI) h_nodipangle_cut_inv_m_PHI -> Fill(d_inv_m);

                // Cut on dip angle
                double d_dip_angle = TMath::ACos((d_pT0*d_pT1+d_pz0*d_pz1) / (d_mom0*d_mom1) );
                if(b_PHI) h_dip_angle -> Fill(d_dip_angle);
                if(d_dip_angle < 0.04) b_PHI = false;

                //mass of hypothosized mother particle (GeV)
                double d_mother_m = -9999;
                if(b_PHI)
                  {
                    h_inv_m_PHI    -> Fill(d_inv_m);
                  }//if(b_PHI)
                if(b_RHO)
                  {
                    h_inv_m_RHO    -> Fill(d_inv_m);
                  }//if(b_RHO)
                if(b_K0S)
                  {
                    h_inv_m_K0S    -> Fill(d_inv_m);
                  }//if(b_K0S)
                if(b_LAMBDA)
                  {
                    h_inv_m_LAMBDA -> Fill(d_inv_m);
                  }//if(b_LAMBDA)
              }//for(int j = i+1; j < v_tracks.size(); j++)
          }//for(int i = 0; i < v_tracks.size()-1;i++))

        // testing purposes, repeat for primary tracks
        // if(v_pri_tracks_pl.size() <= 0) return;
        // if(v_pri_tracks_mi.size() <= 0) return;
        vector<StPicoTrack *> v_pri_tracks0 = v_pri_tracks_pl;
        vector<StPicoTrack *> v_pri_tracks1 = v_pri_tracks_mi;


        for(int i = 0; i < v_pri_tracks0.size();i++)

            //  for(int i = 0; i < v_pri_tracks0.size()-1;i++)
          {
            //      cout<<" . "<<endl;
            StPicoTrack * picoTrack0 = v_pri_tracks0[i];
            if(!picoTrack0) continue;

            double d_TPCnSigmaPion0   = fabs(picoTrack0->nSigmaPion());
            double d_TPCnSigmaProton0 = fabs(picoTrack0->nSigmaProton());
            double d_TPCnSigmaKaon0   = fabs(picoTrack0->nSigmaKaon());

            // PID Cuts
            bool b_PI0  = fabs(d_TPCnSigmaPion0) < d_SigmaCutLevel;
            bool b_PRO0 = fabs(d_TPCnSigmaProton0) < d_SigmaCutLevel;
            bool b_K0   = fabs(d_TPCnSigmaKaon0) < d_SigmaCutLevel;

            //pick the most likely particle from the TPC nsigmas
            if( b_PI0
               && (d_TPCnSigmaPion0 < d_TPCnSigmaProton0)
               && (d_TPCnSigmaPion0 < d_TPCnSigmaKaon0) )
              { b_PI0 = true; b_PRO0 = false;}// b_K0 = false;}

            if( b_PRO0
               && (d_TPCnSigmaProton0 < d_TPCnSigmaPion0)
               && (d_TPCnSigmaProton0 < d_TPCnSigmaKaon0) )
              { b_PRO0 = true; b_PI0 = false;}// b_K0 = false;}

            if( b_K0
               && (d_TPCnSigmaKaon0 < d_TPCnSigmaProton0)
               && (d_TPCnSigmaKaon0 < d_TPCnSigmaPion0) )
              { b_K0 = true; b_PRO0 = false; b_PI0 = false;}

            //double d_TOFnSigmaKaon0   = (picoTrack0->btofPidTraits()).sigmaKaon();
            //bool b_K_tof0   = fabs(d_TOFnSigmaKaon0) < d_SigmaCutLevel;
            //tof cuts
            //      if(!b_K_tof0) b_K0 = false;

            double d_M0 = -9999.0;
            double d_px0      = picoTrack0->pMom().x();
            double d_py0      = picoTrack0->pMom().y();
            double d_pz0      = picoTrack0->pMom().z();
            double d_pT0      = picoTrack0->pPt();
            double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
            StPicoBTofPidTraits *trait = NULL;
            if(picoTrack0->isTofTrack()) trait = dst->btofPidTraits(picoTrack0->bTofPidTraitsIndex());
            double d_tofBeta0 = -999;
            if(trait) d_tofBeta0 = trait->btofBeta();
            //double d_dEdxFit0 = picoTrack0->probPidTraits().dEdxFit();
            //double d_dEdxTru0 = picoTrack0->probPidTraits().dEdxTruncated();
            double d_charge0  = picoTrack0->charge();
            double d_mc0      = d_mom0/d_charge0;
            double d_eta0     = picoTrack0->pMom().Eta();
            double d_phi0     = picoTrack0->pMom().Phi();

            TVector3 v3D_obj_p0 = picoTrack0->pMom();

            StPicoPhysicalHelix    trackhelix0 = picoTrack0->helix(B);

            if(b_PRO0)
              {
                d_M0 = d_PRO_m;
              }//if(b_PRO0)
            else if(b_PI0)
              {
                d_M0 = d_PI_m;
              }//else if(b_PI0)
            else if(b_K0)
              {
                d_M0 = d_K_m;
              }//else if(b_K0)

            for(int j = 0; j < v_pri_tracks1.size(); j++)
                //      for(int j = i+1; j < v_pri_tracks1.size(); j++)
              {
                //      cout<<"    .."<<endl;
                StPicoTrack * picoTrack1 = v_pri_tracks1[j];

                if(!picoTrack1 || (picoTrack0->id() == picoTrack1->id())) continue;
                double d_TPCnSigmaPion1   = fabs(picoTrack1->nSigmaPion());
                double d_TPCnSigmaProton1 = fabs(picoTrack1->nSigmaProton());
                double d_TPCnSigmaKaon1   = fabs(picoTrack1->nSigmaKaon());

                // PID Cuts
                bool b_PI1  = fabs(d_TPCnSigmaPion1) < d_SigmaCutLevel;
                bool b_PRO1 = fabs(d_TPCnSigmaProton1) < d_SigmaCutLevel;
                bool b_K1   = fabs(d_TPCnSigmaKaon1) < d_SigmaCutLevel;

                //pick the most likely particle from the TPC nsigmas
                if( b_PI1
                   && (d_TPCnSigmaPion1 < d_TPCnSigmaProton1)
                   && (d_TPCnSigmaPion1 < d_TPCnSigmaKaon1) )
                  { b_PI1 = true; b_PRO1 = false;}// b_K1 = false;}

                if( b_PRO1
                   && (d_TPCnSigmaProton1 < d_TPCnSigmaPion1)
                   && (d_TPCnSigmaProton1 < d_TPCnSigmaKaon1) )
                  { b_PRO1 = true; b_PI1 = false;}// b_K1 = false;}

                if( b_K1
                   && (d_TPCnSigmaKaon1 < d_TPCnSigmaProton1)
                   && (d_TPCnSigmaKaon1 < d_TPCnSigmaPion1) )
                  { b_K1 = true; b_PRO1 = false; b_PI1 = false;}

                //double d_TOFnSigmaKaon1   = (picoTrack1->btofPidTraits()).sigmaKaon();
                //bool b_K_tof1   = fabs(d_TOFnSigmaKaon1) < d_SigmaCutLevel;
                //tof cuts
                // if(!b_K_tof1) b_K1 = false;

                double d_M1 = -9999.0;
                if(b_PRO1) d_M1 = d_PRO_m;
                else if(b_PI1)  d_M1 = d_PI_m;
                else if(b_K1)   d_M1 = d_K_m;

                double d_px1      = picoTrack1->pMom().x();
                double d_py1      = picoTrack1->pMom().y();
                double d_pz1      = picoTrack1->pMom().z();
                double d_pT1      = picoTrack1->pPt();
                double d_mom1     = sqrt(d_pT1*d_pT1 + d_pz1*d_pz1);
                StPicoBTofPidTraits *trait = NULL;
                if(picoTrack1->isTofTrack()) trait = dst->btofPidTraits(picoTrack1->bTofPidTraitsIndex());
                double d_tofBeta1 = -999;
                if(trait) d_tofBeta1 = trait->btofBeta();
                //double d_dEdxFit1 = picoTrack1->probPidTraits().dEdxFit();
                //double d_dEdxTru1 = picoTrack1->probPidTraits().dEdxTruncated();
                double d_charge1  = picoTrack1->charge();

                double d_mc1      = d_mom1/d_charge1;
                double d_eta1     = picoTrack1->pMom().Eta();
                double d_phi1     = picoTrack1->pMom().Phi();

                TVector3 v3D_obj_p1 = picoTrack1->pMom();

                //final cuts
                if(d_charge0 == d_charge1) continue;
                bool b_PHI    = b_K0 && b_K1;
                bool b_RHO    = b_PI0 && b_PI1;
                bool b_K0S    = b_RHO;

                bool b_LAMBDA = (b_PRO0 && b_PI1)||(b_PRO1 && b_PI0);
                bool b_V0 = b_PHI || b_RHO || b_LAMBDA;
                if(!b_V0) continue;

                // computational intensive, do all cuts before
                StPicoPhysicalHelix    trackhelix1 = picoTrack1->helix(B);
                // Get Mom at decayvtx:
                pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);

                TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, d_MagField*kilogauss);
                TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, d_MagField*kilogauss);

                TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
                TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);

                double d_dca_products =(v3D_x_daughter0-v3D_x_daughter1).Mag();
                //      if(d_dca_products > d_cut_dca_daughters) continue;
                if(d_dca_products > d_cut_dca_daughters_lam) b_LAMBDA = false;
                if(d_dca_products > d_cut_dca_daughters_k0s) b_RHO    = false;

                TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)*0.5;
                TVector3 v3D_xvec_decayl = v3D_x_mother - v3D_vtx;
                TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;

                double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);

                //DCA btw mother and vtx
                double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );

                //if(d_dca_mother > d_cut_dca_mother) continue;
                if(d_dca_mother > d_cut_dca_mother_lam) b_LAMBDA = false;
                if(d_dca_mother > d_cut_dca_mother_k0s) b_RHO    = false;

                double d_mother_decay_length =  v3D_xvec_decayl.Mag();
                //if(d_mother_decay_length < d_cut_mother_decay_length) continue;
                if(d_mother_decay_length < d_cut_mother_decay_length_lam) b_LAMBDA = false;
                if(d_mother_decay_length < d_cut_mother_decay_length_k0s) b_K0S    = false;
                if(d_mother_decay_length > d_cut_mother_decay_length_RHO) b_RHO    = false;

                // double d_E0 = sqrt(v3D_p_daughter0.Mag2()+d_M0*d_M0);
                // double d_E1 = sqrt(v3D_p_daughter1.Mag2()+d_M1*d_M1);
                // double d_inv_m = sqrt(d_M0*d_M0
                //             +d_M1*d_M1
                //             +2.0*d_E0*d_E1
                //             -2.0*(v3D_p_daughter0.Dot(v3D_p_daughter1)) );

                double d_E0 = sqrt(v3D_obj_p0.Mag2()+d_M0*d_M0);
                double d_E1 = sqrt(v3D_obj_p1.Mag2()+d_M1*d_M1);
                double d_inv_m = sqrt(d_M0*d_M0
                                      +d_M1*d_M1
                                      +2.0*d_E0*d_E1
                                      -2.0*(v3D_obj_p0.Dot(v3D_obj_p1)) );




                if(d_mother_decay_length > d_cut_mother_decay_length_PHI) b_PHI    = false;

                // Cut on dip angle
                double d_dip_angle = TMath::ACos((d_pT0*d_pT1+d_pz0*d_pz1) / (d_mom0*d_mom1) );
                if(d_dip_angle < 0.04) b_PHI = false;

                //mass of hypothosized mother particle (GeV)
                double d_mother_m = -9999;
                if(b_PHI)
                  {
                    h_prim_inv_m_PHI    -> Fill(d_inv_m);
                  }//if(b_PHI)
              }//for(int j = i+1; j < v_pri_tracks.size(); j++)
          }//for(int i = 0; i < v_pri_tracks.size()-1;i++))

        // END testing purposes, repeat for primary tracks
        //=========================================================
        //  END Invariant Mass Nested Track Loop
        //=========================================================



        //   //      trackinfo.reset();
        //   // trackinfo.i_event = i_event;
        //   // trackinfo.d_zvtx  = d_zvtx;

        //   // trackinfo.runNumber   = runNumber;
        //   // trackinfo.eventNumber = eventNumber;

        //   // if(d_charge < 0) trackinfo.pT      = (-1.0)*d_pT;
        //   // else             trackinfo.pT      = d_pT;
        //   trackinfo.b_PI  = b_PI;
        //   trackinfo.b_PRO = b_PRO;
        //   trackinfo.b_K   = b_K;

        //   trackinfo.px      = d_px;
        //   trackinfo.py      = d_py;
        //   trackinfo.pz      = d_pz;
        //   // // //      trackinfo.tofBeta = d_tofBeta;
        //   // trackinfo.dEdxFit = d_dEdxFit;
        //   // // //trackinfo.d_dEdxTru = d_dEdxTru;
        //   trackinfo.DCA_r   = d_DCA_r;
        //   // trackinfo.eta     = d_eta;
        //   // trackinfo.phi     = d_phi;

        //   h_pT       -> Fill(d_pT);
        //   if((d_tofBeta != 0.0) && (d_charge != 0.0)) h2_dEdx_pq -> Fill(d_mc,1.0/d_tofBeta);

        //   //==================================================
        //   //========== === Fill TTree ========================
        //   trackinfo.nGoodTracks = nGoodTracks;
        //   //      t_tracks         -> Fill();
        //   //========== END Fill TTree ========================

        //   // Loop over tracks again


        //   //pi minus
        //   bool b_pi =
        //     (
        //      ((d_tofBeta < 1.3) && (d_mc < 0.54))
        //      ||((d_tofBeta < 1.1) && (d_mc < 0.9))
        //      ||((d_tofBeta < 1.04) && (d_mc < 1.4))
        //      );
        //   double d_pi_mass = 0.13957018;

        //   if(!b_pi) continue;
        //   if((d_pT > 1) || (d_pT < 0.1)) continue;
        //   double d_E = sqrt(d_mom*d_mom + d_pi_mass*d_pi_mass);
        //   if(d_E - d_pz == 0.0) continue;
        //   double d_y = 0.5*TMath::Log((d_E + d_pz)/(d_E - d_pz));
        //   double d_mT = sqrt(d_pT*d_pT + d_pi_mass*d_pi_mass);
        //   h_rapidity_pm -> Fill(d_y+1.52);

        //   //for now only fill ultra central
        //   if(b_cent_05)  h2_mT_rapidity_pm -> Fill(d_mT - d_pi_mass,d_y);
        // } // while ( tracks

        //End of Make
        i_event++;
        h_evt -> Fill(0.5);

    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    h_evt_vs_cuts->SetBinContent(1,ievtcut[0]);
    h_evt_vs_cuts->SetBinContent(2,ievtcut[1]);
    h_evt_vs_cuts->SetBinContent(3,ievtcut[2]);
    h_evt_vs_cuts->SetBinContent(4,ievtcut[3]);
    h_evt_vs_cuts->SetBinContent(5,ievtcut[4]);
    h_evt_vs_cuts->SetBinContent(6,ievtcut[5]);
    h_trk_vs_cuts->SetBinContent(1,itrkcut[0]);
    h_trk_vs_cuts->SetBinContent(2,itrkcut[1]);
    h_trk_vs_cuts->SetBinContent(3,itrkcut[2]);
    h_trk_vs_cuts->SetBinContent(4,itrkcut[3]);
    h_trk_vs_cuts->SetBinContent(5,itrkcut[4]);
    h_trk_vs_cuts->SetBinContent(6,itrkcut[5]);

    outputFile->cd();

    bool b_write_tree = false;
    if(b_write_tree)  t_tracks   -> Write();
    h_evt_vs_cuts->Write();
    h_trk_vs_cuts->Write();
    t_K ->Write();
    h_evt      -> Write();
    h_vtx      -> Write();
    h_pT       -> Write();
    h2_dEdx_PI_pq -> Write();
    h2_dEdx_PRO_pq -> Write();
    h2_dEdx_K_pq -> Write();
    h2_dEdx_All_pq->Write();
    h2_m2_QA_pq->Write();
    h2_m2_QA_pT->Write();
    h_mult     -> Write();
    h_DCA_r    -> Write();
    h_DCA_PRO    -> Write();
    h_DCA_PIM    -> Write();
    h_DCA_mother_r    -> Write();
    h_rapidity_pm -> Write();
    h2_mT_rapidity_pm -> Write();

    h_nodecaylength_cut_inv_m_PHI -> Write();
    h_nodipangle_cut_inv_m_PHI -> Write();
    h_inv_m_PHI    -> Write();
    h_prim_inv_m_PHI -> Write();
    h_inv_m_RHO    -> Write();
    h_inv_m_K0S    -> Write();
    h_inv_m_LAMBDA -> Write();

    h_decay_length -> Write();

    h_K_DCA_r       -> Write();
    h_K_obj_DCA_r   -> Write();
    h_K_diff_DCA_r  -> Write();


    h2_pidprob_vs_nsigma_PI  -> Write();
    h2_pidprob_vs_nsigma_PRO -> Write();
    h2_pidprob_vs_nsigma_K   -> Write();

    h2_TOF_nsigma_vs_TPC_nsigma_PI  -> Write();
    h2_TOF_nsigma_vs_TPC_nsigma_PRO    -> Write();
    h2_TOF_nsigma_vs_TPC_nsigma_K      -> Write();

    h2_TOF_nsigma_K_vs_PI  -> Write();
    h2_TPC_nsigma_K_vs_PI  -> Write();
    h2_TOF_nsigma_K_vs_PRO -> Write();
    h2_TPC_nsigma_K_vs_PRO -> Write();

    h_PHI_decay_length     -> Write();
    h_dip_angle   -> Write();

    h_K_minus_pT   -> Write();
    h_PI_minus_pT  -> Write();
    h_PRO_minus_pT -> Write();

    h_K_plus_pT   -> Write();
    h_PI_plus_pT  -> Write();
    h_PRO_plus_pT -> Write();

    h_K_plus_y -> Write();
    h_PI_plus_y -> Write();
    h_PRO_plus_y -> Write();

    h_K_minus_y -> Write();
    h_PI_minus_y -> Write();
    h_PRO_minus_y -> Write();

    h2_K_plus_pT_vs_y -> Write();
    h2_PI_plus_pT_vs_y -> Write();
    h2_PRO_plus_pT_vs_y -> Write();

    h2_K_minus_pT_vs_y -> Write();
    h2_PI_minus_pT_vs_y -> Write();
    h2_PRO_minus_pT_vs_y -> Write();


    gROOT->GetListOfFiles()->Remove(outputFile);
    outputFile->Close();

    /*
    //---
    TH1D * h_inv_m_PHI    = new TH1D("h_inv_m_PHI","h_inv_m_PHI",10000,0.0,10.0);
    TH1D * h_inv_m_K0S    = new TH1D("h_inv_m_K0S","h_inv_m_K0S",1000,0.0,1.0);
    TH1D * h_inv_m_RHO    = new TH1D("h_inv_m_RHO","h_inv_m_RHO",2000,0.45,1.2);
    TH1D * h_inv_m_LAMBDA = new TH1D("h_inv_m_LAMBDA","h_inv_m_LAMBDA",1000,1.0,1.18);
    //---

    TH1D * h_diff_zvtx    = new TH1D("h_diff_zvtx","h_diff_zvtx",100,-0.6,0.6);

    //  TFile * tf_in= new TFile("out_01.root","READ");
    //  TFile * tf_in = new TFile("/star/data01/pwg/jbryslawskyj/out_20K.root","READ");
    TFile * tf_in = new TFile("/star/data01/pwg/jbryslawskyj/SUM_out.root","READ");
    // TH1D * h_inv_m_LAMBDA = (TH1D*) tf_in -> Get("h_inv_m_LAMBDA");
    // h_inv_m_LAMBDA -> Draw();

    TTree * t_tracks = (TTree*) tf_in -> Get("t_tracks");
    int N_entries = t_tracks->GetEntries();
    cout<<" N: "<<N_entries<<" "<<t_tracks<<endl;
    if(N_entries <= 0) { cout<<" ! NO Track ENTRIES ! STOPPING"<<endl; return;}

    vector<StPhysicalHelixD> v_helix_pl;
    vector<StPhysicalHelixD> v_helix_mi;

    st_track *m_trackinfo = new st_track();
    StPhysicalHelixD trackhelix0;
    StPhysicalHelixD trackhelix1;
    // TBranch *branch  = t_tracks->GetBranch("trackinfo");
    // branch->SetAddress(&m_trackinfo);
    // branch -> GetEntry(0);
    // StPhysicalHelixD tmp_helix  =  m_trackinfo->trackhelix;
    // //  cout<<"i: "<<tmp_helix.period()<<endl;

    t_tracks -> SetBranchAddress("trackinfo",&m_trackinfo);


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
    double d_cut_mother_decay_length_RHO = 0.; // must be LESS    than this
    double d_cut_mother_decay_length_PHI = 1.1;//0.1; // must be LESS    than this

    double d_PI_m    = 0.13957018;
    double d_PRO_m   = 0.9382720813;
    double d_K_m     = 0.493677;

    int j_start = 0;
    unsigned long long pre_event_pl_runnumber = 9999;
    unsigned long long this_mixed_event_pl_runnumber = 9999;

    // regular events
    // for(int i = 0; i < N_entries; i++)
    //   {
    //     t_tracks -> GetEntry(i);
    //     bool b_pos_charge0 = t_tracks -> GetLeaf("b_pos_charge")->GetValue(0);
    //     bool b_PI0  = t_tracks -> GetLeaf("b_PI")->GetValue(0);
    //     bool b_PRO0 = t_tracks -> GetLeaf("b_PRO")->GetValue(0);
    //     bool b_K0   = t_tracks -> GetLeaf("b_K")->GetValue(0);

    //     double d_M0 = -9999.0;
    //     if(b_PI0)  d_M0 = d_PI_m;
    //     if(b_PRO0) d_M0 = d_PRO_m;
    //     if(b_K0)   d_M0 = d_K_m;

    //     unsigned int i_runNumber0   = t_tracks -> GetLeaf("runNumber")->GetValue(0);
    //     unsigned int i_eventNumber0 = t_tracks -> GetLeaf("eventNumber")->GetValue(0);
    //     unsigned long long event_pl_runnumber0 = (i_runNumber0-16140000)*100000+i_eventNumber0;

    //     if(event_pl_runnumber0 != pre_event_pl_runnumber) {j_start = i; pre_event_pl_runnumber = event_pl_runnumber0;}

    //     //    branch -> GetEntry(i);
    //     trackhelix0  =  m_trackinfo->trackhelix;
    //     //      cout<<"i: "<<i<<trackhelix0.phase()<<endl;
    //     //      cout<<"  "<<i<<" "<<i_runNumber<<" "<<i_eventNumber<<" "<<event_pl_runnumber<< endl;

    //     double d_MagField0 = t_tracks -> GetLeaf("MagField")->GetValue(0);

    //     float d_xvtx = m_trackinfo -> x_vtx;
    //     float d_yvtx = m_trackinfo -> y_vtx;
    //     float d_zvtx = m_trackinfo -> z_vtx;

    //     // float d_xvtx = 0.0;
    //     // float d_yvtx = 0.0;
    //     // float d_zvtx = 0.0;

    //     StThreeVectorD v3D_vtx0;
    //     v3D_vtx0.setX(d_xvtx);
    //     v3D_vtx0.setY(d_yvtx);
    //     v3D_vtx0.setZ(d_zvtx);
    //     if(!b_pos_charge0) continue;

    //     for(int j = j_start; j < N_entries; j++)
    //         {
    //       t_tracks -> GetEntry(j);
    //           unsigned int i_runNumber1   = t_tracks -> GetLeaf("runNumber")->GetValue(0);
    //           unsigned int i_eventNumber1 = t_tracks -> GetLeaf("eventNumber")->GetValue(0);
    //           unsigned long long event_pl_runnumber1 = (i_runNumber1-16140000)*100000+i_eventNumber1;

    //           if(event_pl_runnumber0 != event_pl_runnumber1) break;//continue;

    //           bool b_pos_charge1 = t_tracks -> GetLeaf("b_pos_charge")->GetValue(0);
    //           bool b_PI1  = t_tracks -> GetLeaf("b_PI")->GetValue(0);
    //           bool b_PRO1 = t_tracks -> GetLeaf("b_PRO")->GetValue(0);
    //           bool b_K1   = t_tracks -> GetLeaf("b_K")->GetValue(0);
    //           //final cuts
    //           if(b_pos_charge0 == b_pos_charge1) continue;

    //       if(b_pos_charge1) continue;

    //           bool b_PHI    = b_K0 && b_K1;
    //           bool b_RHO    = b_PI0 && b_PI1;
    //           bool b_K0S    = b_RHO;

    //           bool b_LAMBDA = (b_PRO0 && b_PI1)||(b_PRO1 && b_PI0);
    //           bool b_V0 = b_PHI || b_RHO || b_LAMBDA;
    //           if(!b_V0) continue;


    //           double d_M1 = -9999.0;
    //           if(b_PI1)  d_M1 = d_PI_m;
    //           if(b_PRO1) d_M1 = d_PRO_m;
    //           if(b_K1)   d_M1 = d_K_m;

    //           trackhelix1  =  m_trackinfo->trackhelix;

    //           //--- --- --- ---
    //           // Get Mom at decayvtx:
    //           pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);

    //           TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, d_MagField0*kilogauss);
    //           TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, d_MagField0*kilogauss);

    //           TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
    //           TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);

    //           double d_dca_products =(v3D_x_daughter0-v3D_x_daughter1).Mag();
    //           //      if(d_dca_products > d_cut_dca_daughters) continue;
    //           if(d_dca_products > d_cut_dca_daughters_lam) b_LAMBDA = false;
    //           if(d_dca_products > d_cut_dca_daughters_k0s) b_RHO    = false;

    //           TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)/2.0;
    //           TVector3 v3D_xvec_decayl = v3D_x_mother - v3D_vtx0;
    //           TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;

    //           double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);

    //           //DCA btw mother and vtx
    //           double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );

    //           //if(d_dca_mother > d_cut_dca_mother) continue;
    //           if(d_dca_mother > d_cut_dca_mother_lam) b_LAMBDA = false;
    //           if(d_dca_mother > d_cut_dca_mother_k0s) b_RHO    = false;

    //           double d_mother_decay_length =  v3D_xvec_decayl.Mag();
    //           // ////h_decay_length -> Fill(d_mother_decay_length);
    //           //if(d_mother_decay_length < d_cut_mother_decay_length) continue;
    //           if(d_mother_decay_length < d_cut_mother_decay_length_lam) b_LAMBDA = false;
    //           if(d_mother_decay_length < d_cut_mother_decay_length_k0s) b_K0S    = false;
    //           if(d_mother_decay_length > d_cut_mother_decay_length_RHO) b_RHO    = false;
    //           if(d_mother_decay_length > d_cut_mother_decay_length_PHI) b_PHI    = false;

    //           double d_E0 = sqrt(v3D_p_daughter0.Mag2()+d_M0*d_M0);
    //           double d_E1 = sqrt(v3D_p_daughter1.Mag2()+d_M1*d_M1);
    //           double d_inv_m = sqrt(d_M0*d_M0
    //                       +d_M1*d_M1
    //                       +2.0*d_E0*d_E1
    //                       -2.0*(v3D_p_daughter0.Dot(v3D_p_daughter1)) );

    //           if(b_PHI)
    //             {
    //               h_inv_m_PHI    -> Fill(d_inv_m);
    //             }//if(b_PHI)
    //           if(b_RHO)
    //             {
    //               h_inv_m_RHO    -> Fill(d_inv_m);
    //             }//if(b_RHO)
    //           if(b_K0S)
    //             {
    //               h_inv_m_K0S    -> Fill(d_inv_m);
    //             }//if(b_K0S)
    //           if(b_LAMBDA)
    //             {
    //               h_inv_m_LAMBDA -> Fill(d_inv_m);
    //             }//if(b_LAMBDA)

    //         }//for(int j = 0; j < N_entries; j++)
    //   }//for(int i = 0; i < N_entries; i++)

    //mixed events
    // for(int i = 0; i < N_entries; i++)
    //   {
    //     t_tracks -> GetEntry(i);

    //   }//for(int i = 0; i < N_entries; i++)
    //  N_entries = 50000;
    //  N_entries = 10 000 000;
    //make a list of events used in mixed events
    bool b_found_mixed_evt = false;
    bool b_half  = false;
    bool b_quart = false;
    std::set<unsigned long long> s_used_mx_events;
    for(int i = 0; i < N_entries; i++)
      {
        if(( i > ((double) N_entries/4.0))&&(!b_quart)&&(!b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
        if(( i > ((double) (3.0*N_entries)/4.0))&&(!b_quart)&&(b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
        if(( i > ((double) N_entries/2.0))&&(!b_half)) {cout<<" halfway"<<endl; b_half = true; b_quart = false;}
        t_tracks -> GetEntry(i);
        bool b_pos_charge0 = t_tracks -> GetLeaf("b_pos_charge")->GetValue(0);
        bool b_PI0  = t_tracks -> GetLeaf("b_PI")->GetValue(0);
        bool b_PRO0 = t_tracks -> GetLeaf("b_PRO")->GetValue(0);
        bool b_K0   = t_tracks -> GetLeaf("b_K")->GetValue(0);

        double d_M0 = -9999.0;
        if(b_PI0)  d_M0 = d_PI_m;
        if(b_PRO0) d_M0 = d_PRO_m;
        if(b_K0)   d_M0 = d_K_m;

        unsigned int i_runNumber0   = t_tracks -> GetLeaf("runNumber")->GetValue(0);
        unsigned int i_eventNumber0 = t_tracks -> GetLeaf("eventNumber")->GetValue(0);
        unsigned long long event_pl_runnumber0 = (i_runNumber0-16140000)*100000+i_eventNumber0;

        //cout<<endl<<" 0: "<<event_pl_runnumber0<<" pre: "<<pre_event_pl_runnumber<<" found: "<<b_found_mixed_evt<<endl;

        if(event_pl_runnumber0 != pre_event_pl_runnumber) {j_start = i; pre_event_pl_runnumber = event_pl_runnumber0; b_found_mixed_evt = false; }

        //    branch -> GetEntry(i);
        trackhelix0  =  m_trackinfo->trackhelix;
        //      cout<<"i: "<<i<<trackhelix0.phase()<<endl;
        //      cout<<"  "<<i<<" "<<i_runNumber<<" "<<i_eventNumber<<" "<<event_pl_runnumber<< endl;

        double d_MagField0 = t_tracks -> GetLeaf("MagField")->GetValue(0);

        float d_xvtx0 = m_trackinfo -> x_vtx;
        float d_yvtx0 = m_trackinfo -> y_vtx;
        float d_zvtx0 = m_trackinfo -> z_vtx;

        // float d_xvtx = 0.0;
        // float d_yvtx = 0.0;
        // float d_zvtx = 0.0;

        StThreeVectorD v3D_vtx0;
        v3D_vtx0.setX(d_xvtx0);
        v3D_vtx0.setY(d_yvtx0);
        v3D_vtx0.setZ(d_zvtx0);
        if(!b_pos_charge0) continue;

        for(int j = j_start; j < N_entries; j++)
          {
            t_tracks -> GetEntry(j);
            unsigned int i_runNumber1   = t_tracks -> GetLeaf("runNumber")->GetValue(0);
            unsigned int i_eventNumber1 = t_tracks -> GetLeaf("eventNumber")->GetValue(0);
            unsigned long long event_pl_runnumber1 = (i_runNumber1-16140000)*100000+i_eventNumber1;

            float d_xvtx1 = m_trackinfo -> x_vtx;
            float d_yvtx1 = m_trackinfo -> y_vtx;
            float d_zvtx1 = m_trackinfo -> z_vtx;

            StThreeVectorD v3D_vtx1;
            v3D_vtx1.setX(d_xvtx1);
            v3D_vtx1.setY(d_yvtx1);
            v3D_vtx1.setZ(d_zvtx1);
            //cout<<"         1precheck: "<<event_pl_runnumber1<<" this: "<<this_mixed_event_pl_runnumber<<" found: "<<b_found_mixed_evt<<endl;
            //if(event_pl_runnumber0 != event_pl_runnumber1) {h_diff_zvtx -> Fill( d_zvtx1 - d_zvtx0); break;}
            if(event_pl_runnumber0 == event_pl_runnumber1)
              {
                //cout<<"                            same event continueing"<<endl;
                continue;
              }
            if(b_found_mixed_evt && (event_pl_runnumber1 != this_mixed_event_pl_runnumber))
              {
                //cout<<" hit next mixed event breaking "<<event_pl_runnumber1 <<endl;
                break;
              }
            if(fabs(d_zvtx1 - d_zvtx0)>0.5)
              {
                // cout<<"                            bad zvtx continuing "<<event_pl_runnumber1<<endl;
                continue;
              }
            // if( (!b_found_mixed_evt) && (s_used_mx_events.find(event_pl_runnumber1)!=s_used_mx_events.end() ) )
            //   {
            //     //          cout<<"                            event used continuning "<<endl;
            //     continue;
            //   }

            //-- found good event
            //if(!b_found_mixed_evt) cout<<"      found 1: "<<event_pl_runnumber1<<" this: "<<this_mixed_event_pl_runnumber<<" 0: "<<event_pl_runnumber0<<endl;
            if(!b_found_mixed_evt) s_used_mx_events.insert(event_pl_runnumber1);
            b_found_mixed_evt = true;
            this_mixed_event_pl_runnumber = event_pl_runnumber1;
            //      cout<<"           running mixed event ana with 1: "<<event_pl_runnumber1<<" 0: "<<event_pl_runnumber0<<endl;
            bool b_pos_charge1 = t_tracks -> GetLeaf("b_pos_charge")->GetValue(0);
            bool b_PI1  = t_tracks -> GetLeaf("b_PI")->GetValue(0);
            bool b_PRO1 = t_tracks -> GetLeaf("b_PRO")->GetValue(0);
            bool b_K1   = t_tracks -> GetLeaf("b_K")->GetValue(0);
            //final cuts
            if(b_pos_charge0 == b_pos_charge1) continue;

            if(b_pos_charge1) continue;

            bool b_PHI    = b_K0 && b_K1;
            bool b_RHO    = b_PI0 && b_PI1;
            bool b_K0S    = b_RHO;

            bool b_LAMBDA = (b_PRO0 && b_PI1)||(b_PRO1 && b_PI0);
            bool b_V0 = b_PHI || b_RHO || b_LAMBDA;
            if(!b_V0) continue;


            double d_M1 = -9999.0;
            if(b_PI1)  d_M1 = d_PI_m;
            if(b_PRO1) d_M1 = d_PRO_m;
            if(b_K1)   d_M1 = d_K_m;

            trackhelix1  =  m_trackinfo->trackhelix;

            //--- --- --- ---
            // Get Mom at decayvtx:
            pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);

            TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, d_MagField0*kilogauss);
            TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, d_MagField0*kilogauss);

            TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
            TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);

            double d_dca_products =(v3D_x_daughter0-v3D_x_daughter1).Mag();
            //      if(d_dca_products > d_cut_dca_daughters) continue;
            if(d_dca_products > d_cut_dca_daughters_lam) b_LAMBDA = false;
            if(d_dca_products > d_cut_dca_daughters_k0s) b_RHO    = false;

            TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)/2.0;
            TVector3 v3D_xvec_decayl = v3D_x_mother - v3D_vtx0;
            TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;

            double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);

            //DCA btw mother and vtx
            double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );

            //if(d_dca_mother > d_cut_dca_mother) continue;
            if(d_dca_mother > d_cut_dca_mother_lam) b_LAMBDA = false;
            if(d_dca_mother > d_cut_dca_mother_k0s) b_RHO    = false;

            double d_mother_decay_length =  v3D_xvec_decayl.Mag();
            // ////h_decay_length -> Fill(d_mother_decay_length);
            //if(d_mother_decay_length < d_cut_mother_decay_length) continue;
            if(d_mother_decay_length < d_cut_mother_decay_length_lam) b_LAMBDA = false;
            if(d_mother_decay_length < d_cut_mother_decay_length_k0s) b_K0S    = false;
            if(d_mother_decay_length > d_cut_mother_decay_length_RHO) b_RHO    = false;
            if(d_mother_decay_length > d_cut_mother_decay_length_PHI) b_PHI    = false;

            double d_E0 = sqrt(v3D_p_daughter0.Mag2()+d_M0*d_M0);
            double d_E1 = sqrt(v3D_p_daughter1.Mag2()+d_M1*d_M1);
            double d_inv_m = sqrt(d_M0*d_M0
                                  +d_M1*d_M1
                                  +2.0*d_E0*d_E1
                                  -2.0*(v3D_p_daughter0.Dot(v3D_p_daughter1)) );

            if(b_PHI)
              {
                h_inv_m_PHI    -> Fill(d_inv_m);
              }//if(b_PHI)
            if(b_RHO)
              {
                h_inv_m_RHO    -> Fill(d_inv_m);
              }//if(b_RHO)
            if(b_K0S)
              {
                h_inv_m_K0S    -> Fill(d_inv_m);
              }//if(b_K0S)
            if(b_LAMBDA)
              {
                h_inv_m_LAMBDA -> Fill(d_inv_m);
              }//if(b_LAMBDA)

          }//for(int j = 0; j < N_entries; j++)
      }//for(int i = 0; i < N_entries; i++)

    TFile * tf_inv_out = new TFile("inv_out.root","RECREATE");
    h_inv_m_PHI    -> Write();
    h_inv_m_K0S    -> Write();
    h_inv_m_RHO    -> Write();
    h_inv_m_LAMBDA -> Write();
    h_diff_zvtx    -> Write();
    */

  picoReader->Finish();

    //outputFile->Write();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    << std::endl;
}
