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

bool b_K_min_check = false;//true; //turn off pid for k minus to get more statistics //  What is k minus?
bool b_K_pl_check  = false;//Tue Sep  4 21:19:14 EDT 2018 //true; //turn off pid for k minus to get more statistics
// save tracks for inv mass calc // Tracks: does the track here mean channel or trajectory?
struct st_track         //Build a structure named st_strack to combine data of different kinds.
{
  bool    b_pos_charge; //What is b_pos_charge?
  bool    b_bad_TOF;    //What is b_bad_TOF?
  int     runNumber;    //What is run?
  int     eventNumber;
  int     nGoodTracks;  //What is nGoodTracks
  double  px;           //momentum-x
  double  py;           //momentum-y
  double  pz;
  double  x_vtx;        //What is x_vtx, velocity?
  double  y_vtx;
  double  z_vtx;
  double  DCA_r;        //What is DCA_r
  double  d_TPCnSigmaKaon; //What is d_TPCnSigmaKaon?
  double  d_TOFnSigmaKaon; //What is d_TOFnSigmaKaon?
};//struct st_track

struct st_event              //A structure named st_event
{
  vector<st_track> v_trk_pl; //What is v_trk_pl
  vector<st_track> v_trk_mi; /*What is v_trk_mi*/
};//struct st_event

void read_K_tree()           /*Main body*/
{
  //K +- mass
  const  double d_K_m     = 0.493677;   /*Kaon mass in GeV*/
  //old invariant mass binning: 1000,0.9,1.1    -> binwidth = 0.0002 GeV
  //binning from STAR paper: 5*9 = 45,0.98,1.07 -> binwidth = 0.002  GeV

  //Mon Aug 28 10:36:25 EDT 2017
  TH1D * h_prim_inv_m_PHI    = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",100,0.9,1.1);     /*histogram_primary?_invariant_mass_$\phi$*/
  TH1D * h_mx_prim_inv_m_PHI = new TH1D("h_mx_prim_inv_m_PHI","h_prim_inv_m_PHI",100,0.9,1.1);  //BG mixed

  //  TH1D * h_dip_angle = new TH1D("h_dip_angle","h_dip_angle",200,0.0,2.0);
  TH1D * h_dip_angle = new TH1D("h_dip_angle","h_dip_angle",1000,-1,1.0);                       //Dip angle of What? Magnetic field?

  TH1D * h_Kmin_TPC2_TOF3_pT = new TH1D("h_Kmin_TPC2_TOF3_pT","h_Kmin_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_Kpl_TPC2_TOF3_pT  = new TH1D("h_Kpl_TPC2_TOF3_pT","h_Kpl_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_phi_TPC2_TOF3_pT  = new TH1D("h_phi_TPC2_TOF3_pT","h_phi_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_phi_TPC2_TOF3_mT  = new TH1D("h_phi_TPC2_TOF3_mT","h_phi_TPC2_TOF3_mT",200,0.0,10);

  TH1D * h_phi_TPC2_TOF3_y  = new TH1D("h_phi_TPC2_TOF3_y","h_phi_TPC2_TOF3_y",200,-10.,10);
  TH1D * h_phi_TPC2_TOF3_tight_y  = new TH1D("h_phi_TPC2_TOF3_tight_y","h_phi_TPC2_TOF3_tight_y",200,-10.0,10.0);

  //  TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_02_08_2017.root","READ");
  //  TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_11_08_2017.root","READ");

  //// TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_20_08_2017.root","READ");

  /////TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_10_02_2018.root","READ");
  //k  TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_09_08_2018.root","READ");
  //  TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_10_08_2018.root","READ");

  //  TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_31_08_2018.root","READ");
  //    TFile * tf_in  = new TFile("/home/jbryslaw/LES_Ana/SUM_03_09_2018.root","READ");
  TFile * tf_in  = new TFile("Phi_test.picoDst.result.root","READ"); //same as SUM_03_09_2018.root but on i5 //output

  TTree * t_K = (TTree*) tf_in -> Get("t_K");                                               //Get the tree_K from events
  int N_entries = t_K -> GetEntries();                                                      //How many dots

  int pre_runNumber   = -9999;                                                              //Why negative?
  int pre_eventNumber = -9999;

  vector<st_track> v_trk_pl;
  vector<st_track> v_trk_mi;

  vector<st_event> v_evt;

  TH2D * h2_TOF_beta_pq                  = new TH2D("h2_TOF_beta_pq","h2_TOF_beta_pq",500,-2,2,500,0,3);
  TH2D * h2_TOF_beta_pq1                 = new TH2D("h2_TOF_beta_pq1","h2_TOF_beta_pq1",500,-2,2,500,0,3);
  TH2D * h2_TOF_beta_pq2                 = new TH2D("h2_TOF_beta_pq2","h2_TOF_beta_pq2",500,-2,2,500,0,3);

  TH2D * h2_dEdx_pq                      = new TH2D("h2_dEdx_pq","h2_dEdx_pq",500,-2,2,500,2e-6,10e-6);

  TH3D * h3_TOF_beta_pq_vs_nsig          = new TH3D("h3_TOF_beta_pq_vs_nsig","h3_TOF_beta_pq_vs_nsig",500,-2,2,500,0,3,100,0.0,10.0);
  TH3D * h3_dEdx_pq_vs_nsig              = new TH3D("h3_dEdx_pq_vs_nsig","h3_dEdx_pq_vs_nsig",500,-2,2,500,2e-6,10e-6,100,0.0,10.0);

  TH2D * h2_TOF_beta_TPC2sig_TOF_3sig_pq = new TH2D("h2_TOF_beta_TPC2sig_TOF_3sig_pq","h2_TOF_beta_TPC2sig_TOF_3sig_pq",500,-2,2,500,0,3);
  TH2D * h2_dEdx_TPC2sig_TOF_3sig_pq     = new TH2D("h2_dEdx_TPC2sig_TOF_3sig_pq","h2_dEdx_TPC2sig_TOF_3sig_pq",500,-2,2,500,2e-6,10e-6);

  TH1D * h_K_pl_pT  = new TH1D("h_K_pl_pT","h_K_pl_pT",200,0.0,10);
  TH1D * h_K_min_pT = new TH1D("h_K_min_pT","h_K_min_pT",200,0.0,10);


  TH1D * h_K_pl_eta  = new TH1D("h_K_pl_eta","h_K_pl_eta",200,-10,10);
  TH1D * h_K_min_eta = new TH1D("h_K_min_eta","h_K_min_eta",200,-10,10);

  int N_events = 0;
  long int N_max_events =1e16;//1000000;//1e16;//1000000;//1e16;
  for(long int i = 0; i< N_entries; i++)
    {
      //      if(N_events >= 3) break;
      //if(i > 50) break;
      if(i > N_max_events) break;
      //      if(i< N_entries-  50) continue;
      t_K -> GetEntry(i);
      bool    b_pos_charge    = t_K-> GetLeaf("b_pos_charge")->GetValue(0);
      bool    b_bad_TOF       = t_K-> GetLeaf("b_bad_TOF")->GetValue(0);
      int     runNumber       = t_K-> GetLeaf("runNumber")->GetValue(0);
      int     eventNumber     = t_K-> GetLeaf("eventNumber")->GetValue(0);
      int     nGoodTracks     = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
      double  px              = t_K-> GetLeaf("px")->GetValue(0);
      double  py              = t_K-> GetLeaf("py")->GetValue(0);
      double  pz              = t_K-> GetLeaf("pz")->GetValue(0);
      double  x_vtx           = t_K-> GetLeaf("x_vtx")->GetValue(0);
      double  y_vtx           = t_K-> GetLeaf("y_vtx")->GetValue(0);
      double  z_vtx           = t_K-> GetLeaf("z_vtx")->GetValue(0);
      double  DCA_r           = t_K-> GetLeaf("DCA_r")->GetValue(0);
      double  d_TPCnSigmaKaon = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
      double  d_TOFnSigmaKaon = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);

      double pT  = sqrt(px*px+py*py);
      double mom = sqrt(px*px+py*py+pz*pz);
      double E   = sqrt((px*px+py*py+pz*pz)+d_K_m*d_K_m);
      double y   = ((E-pz) != 0.0) ? 0.5*TMath::Log( (E + pz) / (E - pz) ) : -9999;
      double eta = ((mom - pz) != 0.0) ? 0.5*TMath::Log( (mom + pz) / (mom - pz) ) : -9999.0;


      //after 10_02_2018
      double d_tofBeta        = t_K -> GetLeaf("tofBeta")->GetValue(0);
      double dEdxFit          = t_K -> GetLeaf("dEdxFit")->GetValue(0);
      double d_dEdxTru        = t_K -> GetLeaf("d_dEdxTru")->GetValue(0);
      double d_inv_tofBeta    = -9999.0;
      if(d_tofBeta != 0.0) d_inv_tofBeta = 1.0 / d_tofBeta;

      //      cout<<" "<<i<<" "<<runNumber<<" "<<eventNumber;
      bool b_new_event = false;
      if((pre_runNumber != runNumber)||(pre_eventNumber != eventNumber))
	{
	  N_events++;
	  b_new_event     = true;
	  pre_runNumber   = runNumber;
	  pre_eventNumber = eventNumber;

	  st_event evt;
	  evt.v_trk_pl = v_trk_pl;
	  evt.v_trk_mi = v_trk_mi;

	  v_evt.push_back(evt);
	  v_trk_pl.clear();
	  v_trk_mi.clear();
	}//if((pre_runNumber != runNumber)&&(pre_eventNumber != eventNumber))

      // if(b_new_event)  cout<<" new event N: "<<N_events<< endl;
      // else cout<<endl;

      // PID QA histos
      double d_q = 1.;
      if(!b_pos_charge) d_q = -1.;
      double d_p = sqrt(px*px+py*py+pz*pz);
      double d_pq = fabs(d_p) * d_q;

      h2_TOF_beta_pq -> Fill(d_pq,d_inv_tofBeta);
      h2_dEdx_pq      -> Fill(d_pq,d_dEdxTru);

      h3_TOF_beta_pq_vs_nsig  -> Fill(d_pq,d_inv_tofBeta,d_TOFnSigmaKaon);
      h3_dEdx_pq_vs_nsig      -> Fill(d_pq,d_dEdxTru,d_TPCnSigmaKaon);

      //if((fabs(d_TPCnSigmaKaon) < 3.0) && (fabs(d_TOFnSigmaKaon) < 2.0))
      if((fabs(d_TPCnSigmaKaon) < 3.0) && (fabs(d_TOFnSigmaKaon) < 2.0)) //test
	{
	  h2_TOF_beta_TPC2sig_TOF_3sig_pq -> Fill(d_pq,d_inv_tofBeta);
	  h2_dEdx_TPC2sig_TOF_3sig_pq     -> Fill(d_pq,d_dEdxTru);
	}//if((d_TPCnSigmaKaon < 3.0) && (d_TOFnSigmaKaon < 2.0))

      st_track trk;
      trk.b_pos_charge    = b_pos_charge;
      trk.b_bad_TOF       = b_bad_TOF;
      trk.runNumber       = runNumber;
      trk.eventNumber     = eventNumber;
      trk.nGoodTracks     = nGoodTracks;
      trk.px              = px;
      trk.py              = py;
      trk.pz              = pz;
      trk.x_vtx           = x_vtx;
      trk.y_vtx           = y_vtx;
      trk.z_vtx           = z_vtx;
      trk.DCA_r           = DCA_r;
      trk.d_TPCnSigmaKaon = d_TPCnSigmaKaon;
      trk.d_TOFnSigmaKaon = d_TOFnSigmaKaon;

      if(b_pos_charge) v_trk_pl.push_back(trk);
      else v_trk_mi.push_back(trk);
      //if(i%1000 == 0) cout<< "v_trk_pl size "<<v_trk_pl.size()<<" "
                         // << "v_trk_mi size "<<v_trk_mi.size()<<" "
                          //<< "i "<<i<<endl;



      //fill Ks to their own histos if they pass the ana cuts
      //eta cuts
      //cout << eta << endl;// for test
     // h_K_min_pT -> Fill(pT); // test
      bool b_K_eta   = (eta >= -1.47) && (eta <= 0.0);//cout
      if(!b_K_eta) continue;

      //y rapidity cuts
      if((y<-2.0)||(y>0.0)) continue;
      //cout /*<< "y rapidity = " */<< y << endl; //test
      //cout << "pT = " << pT << endl;// for test

      //pT cut
      if(pT > 3.0) continue;
      //cout << pT << endl;// for test
      //cout << "d_TPCnSigmaKaon = "<<d_TPCnSigmaKaon << endl;// for test
      //cout << "d_TOFnSigmaKaon = "<<d_TOFnSigmaKaon << endl;// for test
      if((fabs(d_TPCnSigmaKaon) < 3.0) && (fabs(d_TOFnSigmaKaon) < 2.0)) //test
	{
	  if(b_pos_charge) h_K_pl_pT -> Fill(pT);
	  else  h_K_min_pT -> Fill(pT);

	  if(b_pos_charge) h_K_pl_eta -> Fill(eta);
	  else  h_K_min_eta -> Fill(eta);
    //cout << "b_pos_charge =" << b_pos_charge << endl; // test

	}//if((fabs(d_TPCnSigmaKaon) < 3.0) && (fabs(d_TOFnSigmaKaon) < 2.0))
    }//for(long int i = 0; i< N_entries; i++)

  // new TCanvas();
  // h_K_pl_pT  -> Draw();
  // new TCanvas();
  // h_K_min_pT -> Draw();

  //   new TCanvas();
  // h_K_pl_eta  -> Draw();
  // new TCanvas();
  // h_K_min_eta -> Draw();

  // MAIN EVENT LOOP
  // TH1D * h_prim_inv_dummy;
  // unsigned int nGoodTracks_cut;
  // double d_DCA_r_cut;
  // double d_TPCnSigmaKaon_cut;
  // double d_TOFnSigmaKaon_cut;

  // TTree * t_prim_inv = new TTree("t_prim_inv","t_prim_inv");
  // t_prim -> Branch("h_prim_inv",      &h_prim_inv_dummy);
  // t_prim -> Branch("nGoodTracks_cut",     &nGoodTracks_cut,     "nGoodTracks_cut/I");
  // t_prim -> Branch("d_DCA_r_cut",           &d_DCA_r_cut          , "d_DCA_r_cut/D");
  // t_prim -> Branch("d_TPCnSigmaKaon_cut", &d_TPCnSigmaKaon_cut, "d_TPCnSigmaKaon_cut/D");
  // TOF nsigma not working (Sat Aug 12 17:57:34 EDT 2017)

  //  t_prim -> Branch("d_TOFnSigmaKaon_cut", &d_TOFnSigmaKaon_cut, "d_TOFnSigmaKaon_cut/D");

  //later do this for each ngooodtracks bins
  // TH3D * h3_prim_inv_multi_nsig = new TH3D("h3_prim_inv_multi_nsig","Kaon inv m vs Multiplicity vs nsig",
  // 					   100,0.9,1.1,
  // 					   250,0.,250.,
  // 					   100,0,4);
  // TH3D * h3_prim_inv_nsig = new TH3D("h3_prim_inv_nsig","Kaon inv m vs TOF nsig vs TPC nsig",
  // 				     100,0.9,1.1,
  // 				     100,0,4,
  // 				     100,0,4);

  TH3D * a_h3_prim_inv_nsig[8];
  TString TS_cent[8] =
    {
      "Min Bias",
      "0-10% Central",
      "10-30 Central",
      ">30% Central",
      "0.5<pT<1",
      "1<pT<2",
      "2<pT<3",
      "0.5<pT"
    };
    TString TS_ctr[8] =
    {
      "MB",
      "010",
      "1030",
      "30100",
      "pT05_1",
      "pT1_2",
      "pT2_3",
      "pT3_4"
    };

  for(int i = 0; i < 8; i++)
    {
      // a_h3_prim_inv_nsig[i] = new TH3D(Form("h3_prim_inv_nsig_%s",TS_ctr[i].Data()),
      // 				       Form("Kaon inv m vs TOF nsig vs TPC nsig %s",TS_cent[i].Data()),
      // 				       100,0.9,1.1,
      // 				       100,0,4,
      // 				       100,0,4);
      a_h3_prim_inv_nsig[i] = new TH3D(Form("h3_prim_inv_nsig_%s",TS_ctr[i].Data()),
				       Form("Kaon inv m vs TOF nsig vs TPC nsig %s",TS_cent[i].Data()),
				       100,0.9,1.1,
				       100,0,4,
				       100,0,4);


    }//for(int i = 0; i < 8; i++)

  TH1D * h_TPC2_TOF3_invM = new TH1D("h_TPC2_TOF3_invM","M_{K+K-} TPC < 2#sigma TOF < 3#sigma",100,0.9,1.1);


  cout<<endl<<"  LOOPING OVER EVENTS "<< v_evt.size()<<endl<<endl;
  for(int i = 0; i < v_evt.size(); i++)
    {
      //      cout<<" new event"<<endl;
      st_event evt = v_evt[i];
      vector<st_track> v_trk_pl = v_evt[i].v_trk_pl;
      vector<st_track> v_trk_mi = v_evt[i].v_trk_mi;
      for(int j = 0; j < v_trk_pl.size(); j++)
	{
	  st_track trk0 = v_trk_pl[j];

	  bool    b_pos_charge0    = trk0.b_pos_charge;
	  bool    b_bad_TOF0       = trk0.b_bad_TOF;
	  int     runNumber0       = trk0.runNumber;
	  int     eventNumber0     = trk0.eventNumber;
	  int     nGoodTracks0     = trk0.nGoodTracks;
	  double  px0              = trk0.px;
	  double  py0              = trk0.py;
	  double  pz0              = trk0.pz;
	  double  x_vtx0           = trk0.x_vtx;
	  double  y_vtx0           = trk0.y_vtx;
	  double  z_vtx0           = trk0.z_vtx;
	  double  DCA_r0           = trk0.DCA_r;
	  double  d_TPCnSigmaKaon0 = trk0.d_TPCnSigmaKaon;
	  double  d_TOFnSigmaKaon0 = trk0.d_TOFnSigmaKaon;

	  // if(b_K_pl_check)
	  //   {
	  //     d_TPCnSigmaKaon0 = 0.0;
	  //     d_TOFnSigmaKaon0 = 0.0;
	  //   }//if(b_K_pl_check)


	  //cout<<" "<<j<<" "<<runNumber0<<" "<<eventNumber0<<" px0: "<<px0<<endl;
	  double d_pT0 = sqrt(px0*px0+py0*py0);
	  double d_mom0 = sqrt(px0*px0+py0*py0+pz0*pz0);

	  double eta0 = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : 1.0;

	  for(int k = 0; k < v_trk_mi.size(); k++)
	    {
	      st_track trk1 = v_trk_mi[k];

	      bool    b_pos_charge1    = trk1.b_pos_charge;
	      bool    b_bad_TOF1       = trk1.b_bad_TOF;
	      int     runNumber1       = trk1.runNumber;
	      int     eventNumber1     = trk1.eventNumber;
	      int     nGoodTracks1     = trk1.nGoodTracks;
	      double  px1              = trk1.px;
	      double  py1              = trk1.py;
	      double  pz1              = trk1.pz;
	      double  x_vtx1           = trk1.x_vtx;
	      double  y_vtx1           = trk1.y_vtx;
	      double  z_vtx1           = trk1.z_vtx;
	      double  DCA_r1           = trk1.DCA_r;
	      double  d_TPCnSigmaKaon1 = trk1.d_TPCnSigmaKaon;
	      double  d_TOFnSigmaKaon1 = trk1.d_TOFnSigmaKaon;
	      //cout<<" "<<k<<" "<<runNumber1<<" "<<eventNumber1<<" px1: "<<px1<<endl;
	      //turn off pid for k minus to get more statistics
 	      if(b_K_min_check)
		{
		  d_TPCnSigmaKaon1 = 0.0;
		  d_TOFnSigmaKaon1 = 0.0;
		}//

	      double d_pT1 = sqrt(px1*px1+py1*py1);
	      double d_mom1 = sqrt(px1*px1+py1*py1+pz1*pz1);

	      double eta1 = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;
	      //double check events
	      if(eventNumber0 != eventNumber1) continue;
	      if(eventNumber1 != eventNumber1) continue;

	      bool b_PHI = true;

	      // Cut on dip angle
	      double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
	      h_dip_angle -> Fill(d_dip_angle);
	      //	      if(d_dip_angle < 0.04) b_PHI = false;

	      double d_M1 = d_K_m;
	      double d_M0 = d_K_m;
	      double d_E0 = sqrt((px0*px0+py0*py0+pz0*pz0)+d_M0*d_M0);
	      double d_E1 = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);

	      double d_inv_m = sqrt(d_M0*d_M0
				    +d_M1*d_M1
				    +2.0*d_E0*d_E1
				    -2.0*(px0*px1+py0*py1+pz0*pz1) );

	      if(b_PHI) h_prim_inv_m_PHI    -> Fill(d_inv_m);

	      double d_min_DCA     = (DCA_r0 < DCA_r1) ? DCA_r0 : DCA_r1;
	      double d_max_TPCnsig = (fabs(d_TPCnSigmaKaon0) > fabs(d_TPCnSigmaKaon1)) ? d_TPCnSigmaKaon0 : d_TPCnSigmaKaon1;
	      double d_max_TOFnsig = (fabs(d_TOFnSigmaKaon0) > fabs(d_TOFnSigmaKaon1)) ? d_TOFnSigmaKaon0 : d_TOFnSigmaKaon1;
	      int i_max_multi      = (nGoodTracks0 > nGoodTracks1) ? nGoodTracks0 : nGoodTracks1;
	      double d_max_multi   = (double) i_max_multi;

	      // "Min Bias"      -> 0-240
	      // "0-10% Central  -> 121-240
	      // "10-30 Central  -> 48-121
	      // ">30% Central   -> 0-48
	      bool a_b_cent0[4];
	      a_b_cent0[0] = ( nGoodTracks0 > 0)   && ( nGoodTracks0 < 240);
	      a_b_cent0[1] = ( nGoodTracks0 > 121) && ( nGoodTracks0 < 240);
	      a_b_cent0[2] = ( nGoodTracks0 > 48)  && ( nGoodTracks0 <= 121);
	      a_b_cent0[3] = ( nGoodTracks0 > 0)   && ( nGoodTracks0 <= 48);

	      bool a_b_cent1[4];
	      a_b_cent1[0] = ( nGoodTracks1 > 0)   && ( nGoodTracks1 < 240);
	      a_b_cent1[1] = ( nGoodTracks1 > 121) && ( nGoodTracks1 < 240);
	      a_b_cent1[2] = ( nGoodTracks1 > 48)  && ( nGoodTracks1 <= 121);
	      a_b_cent1[3] = ( nGoodTracks1 > 0)   && ( nGoodTracks1 <= 48);

	      bool a_b_pT0[4];
	      a_b_pT0[0]    = (d_pT0 > 0.5) && (d_pT0 <= 1.0);
	      a_b_pT0[1]    = (d_pT0 > 1.0) && (d_pT0 <= 2.0);
	      a_b_pT0[2]    = (d_pT0 > 2.0) && (d_pT0 <= 3.0);
	      a_b_pT0[3]    = (d_pT0 > 0.5);

	      bool a_b_pT1[4];
	      a_b_pT1[0]    = (d_pT1 > 0.5) && (d_pT1 <= 1.0);
	      a_b_pT1[1]    = (d_pT1 > 1.0) && (d_pT1 <= 2.0);
	      a_b_pT1[2]    = (d_pT1 > 2.0) && (d_pT1 <= 3.0);
	      a_b_pT1[3]    = (d_pT1 > 0.5);

	      //eta fiducial cut
	      bool b_K0_eta   = (eta0 >= -1.47) && (eta0 <= 0.0);
	      bool b_K1_eta   = (eta1 >= -1.47) && (eta1 <= 0.0);

	      // pT cuts
	      // K cuts:
	      //	      if((d_pT1 > 2.0) || (d_pT0 > 2.0) ) continue;
	      // phi cuts:
	      double d_pT_phi = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);// px0*py0*px0*py0 + px1*py1*px1*py1 + 2.*px0*px1  + 2.*py0*py1;
	      double m_phi = 1.019455;
	      double d_mT_phi = sqrt(d_pT_phi*d_pT_phi + m_phi*m_phi );

	      double d_phi_pz = pz0+pz1;
	      double d_phi_E  = d_E0+d_E1;
	      double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;

	      h_phi_TPC2_TOF3_y -> Fill(d_phi_y);

	      //if( (!b_K0_eta) || (!b_K1_eta)) continue;
	      //if( (!b_K0_eta )) continue;

	      //if((d_pT_phi < 1.0) || (d_pT_phi > 3.0) ) continue;
	      //if(d_pT_phi > 3.0) continue;
	      if(d_mT_phi > 2.0) continue;
	      //y rapidity cuts
	      //if((d_phi_y<-2.0)||(d_phi_y>0.0)) continue;

	      //only fill if inv m is in range
	      if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;

	      if(nGoodTracks0 != nGoodTracks1) continue;
	      for( int l = 0; l < 4; l++) if(b_PHI && a_b_cent0[l] && a_b_cent1[l]) a_h3_prim_inv_nsig[l] -> Fill(d_inv_m,fabs(d_max_TOFnsig),fabs(d_max_TPCnsig));
	      //make the pT inclusive (one particle in the pair must pass the cuts
	      for( int l = 4; l < 8; l++) if(b_PHI && (a_b_pT0[l-4] || a_b_pT1[l-4]) ) a_h3_prim_inv_nsig[l] -> Fill(d_inv_m,fabs(d_max_TOFnsig),fabs(d_max_TPCnsig));


	      if(b_K_pl_check)
		{
		  d_max_TPCnsig = d_TPCnSigmaKaon1;
		  d_max_TOFnsig = d_TOFnSigmaKaon1;
		}//if(b_K_pl_check)


	      if( (fabs(d_max_TPCnsig) <= 2.0) && (fabs(d_max_TOFnsig) <= 3.0) )//test
		{
		  if(j==0) h_Kmin_TPC2_TOF3_pT -> Fill( d_pT1);//Kmin_pT);
		  if(k==0) h_Kpl_TPC2_TOF3_pT  -> Fill( d_pT0);// Kpl_pT);

		  //Px  = p0x+p1x, Py = p0y+p1y
		  //PT2 = Px2+Py2 = p0xp1x2+2p0xp1x +  p0yp1y2+2p0yp1y
		  if(b_PHI&&a_b_cent0[0]&&a_b_cent1[0])
		    {
		      //4 sigma for TPC2 TOF3 1.00138 : 1.03956
		      if((d_inv_m > 1.00138) && (d_inv_m < 1.03956)) h_phi_TPC2_TOF3_tight_y -> Fill(d_phi_y);
		      h_phi_TPC2_TOF3_pT -> Fill(d_pT_phi);
		      h_phi_TPC2_TOF3_mT -> Fill(d_mT_phi-m_phi);

		      h_TPC2_TOF3_invM -> Fill(d_inv_m);
		    }//if(b_PHI)
		}//if( (d_max_TPCnsig < 2.0) && (d_max_TOFnsig < 3.0) )

	      //	      for( int k = 0; k < 4; k++) if(b_PHI && a_b_cent[k]) a_h3_prim_inv_nsig[k] -> Fill(d_inv_m,d_max_TOFnsig,d_max_TPCnsig);
	      //if(b_PHI) h3_prim_inv_multi_nsig -> Fill(d_inv_m,d_min_DCA,d_min_TPCnsig);
	      // //if(b_PHI) h3_prim_inv_nsig -> Fill(d_inv_m,d_max_TOFnsig,d_max_TPCnsig);
	      //	      if(b_PHI) h3_prim_inv_nsig -> Fill(d_inv_m,d_max_multi,d_max_TPCnsig);
	    }//for(int k = 1; k < v_trk_mi.size(); k++)
	}//for(int j = 0; j < v_trk_pl.size(); j++)
    }//for(int i = 0; i < v_evt.size(); i++)

  h_prim_inv_m_PHI -> SaveAs("h_prim_inv_m_PHI.root");
  for( int k = 0; k < 8; k++) a_h3_prim_inv_nsig[k] -> SaveAs(Form("h3_prim_inv_nsig_%s.root",TS_ctr[k].Data()));
    //h3_prim_inv_nsig -> SaveAs("h3_prim_inv_nsig.root");
  // new TCanvas();
  // h_prim_inv_m_PHI -> Draw();

    h_TPC2_TOF3_invM    -> SaveAs("h_TPC2_TOF3_invM_out.root");
    ///    return;

  // MIXED EVENT LOOP
  // TH3D * h3_mx_prim_inv_multi_nsig = new TH3D("h3_mx_prim_inv_multi_nsig","Mixed Events Kaon inv m vs multi vs nsig",
  // 					      100,0.9,1.1,
  // 					    250,0.,250.,
  // 					      100,0,4);
  // TH3D * h3_mx_prim_inv_nsig = new TH3D("h3_mx_prim_inv_nsig","Kaon inv m vs TOF nsig vs TPC nsig",
  // 				     100,0.9,1.1,
  // 				     100,0,4,
  // 				     100,0,4);
  TH3D * a_h3_mx_prim_inv_nsig[8];
  for(int i = 0; i < 8; i++)
    {
      a_h3_mx_prim_inv_nsig[i] = new TH3D(Form("h3_mx_prim_inv_nsig_%s",TS_ctr[i].Data()),
					  Form("Kaon inv m vs TOF nsig vs TPC nsig %s",TS_cent[i].Data()),
					  100,0.9,1.1,
					  100,0,4,
					  100,0,4);

    }//for(int i = 0; i < 4; i++)

  TH1D * h_mx_TPC2_TOF3_invM = new TH1D("h_mx_TPC2_TOF3_invM","Mixed Events M_{K+K-} TPC < 2#sigma TOF < 3#sigma",100,0.9,1.1);

  int j_start = 0;
  unsigned long long pre_event_pl_runnumber = 9999;
  unsigned long long this_mixed_event_pl_runnumber = 9999;

  //make a list of events used in mixed events
  bool b_found_mixed_evt = false;
  bool b_half  = false;
  bool b_quart = false;
  std::set<unsigned long long> s_used_mx_events;
  for(int i = 0; i < N_entries; i++)
    {
      if(i > N_max_events) break;
      if(( i > ((double) N_entries/4.0))&&(!b_quart)&&(!b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
      if(( i > ((double) (3.0*N_entries)/4.0))&&(!b_quart)&&(b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
      if(( i > ((double) N_entries/2.0))&&(!b_half)) {cout<<" halfway"<<endl; b_half = true; b_quart = false;}
      t_K -> GetEntry(i);
      bool b_pos_charge0 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);

      unsigned int i_runNumber0   = t_K -> GetLeaf("runNumber")->GetValue(0);
      unsigned int i_eventNumber0 = t_K -> GetLeaf("eventNumber")->GetValue(0);
      unsigned long long event_pl_runnumber0 = (i_runNumber0-16140000)*100000+i_eventNumber0;

      if(event_pl_runnumber0 != pre_event_pl_runnumber) {j_start = i; pre_event_pl_runnumber = event_pl_runnumber0; b_found_mixed_evt = false; }

      int     nGoodTracks0     = t_K-> GetLeaf("nGoodTracks")->GetValue(0);

      double d_xvtx0 = t_K -> GetLeaf("x_vtx")->GetValue(0);
      double d_yvtx0 = t_K -> GetLeaf("y_vtx")->GetValue(0);
      double d_zvtx0 = t_K -> GetLeaf("z_vtx")->GetValue(0);

      double px0 = t_K -> GetLeaf("px")->GetValue(0);
      double py0 = t_K -> GetLeaf("py")->GetValue(0);
      double pz0 = t_K -> GetLeaf("pz")->GetValue(0);
      double d_pT0 = sqrt(px0*px0+py0*py0);
      double d_mom0 = sqrt(px0*px0+py0*py0+pz0*pz0);
      double eta0 = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : 1.0;

      double  DCA_r0           = t_K-> GetLeaf("DCA_r")->GetValue(0);
      double  d_TPCnSigmaKaon0 = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
      double  d_TOFnSigmaKaon0 = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);


      // Use only positive chargges in 1st loop?
      //            if(!b_pos_charge0) continue;

      for(int j = j_start; j < N_entries; j++)
      	{
	  t_K -> GetEntry(j);
      	  unsigned int i_runNumber1   = t_K -> GetLeaf("runNumber")->GetValue(0);
      	  unsigned int i_eventNumber1 = t_K -> GetLeaf("eventNumber")->GetValue(0);
      	  unsigned long long event_pl_runnumber1 = (i_runNumber1-16140000)*100000+i_eventNumber1;

	  int     nGoodTracks1     = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
	  double d_xvtx1 = t_K -> GetLeaf("x_vtx")->GetValue(0);
	  double d_yvtx1 = t_K -> GetLeaf("y_vtx")->GetValue(0);
	  double d_zvtx1 = t_K -> GetLeaf("z_vtx")->GetValue(0);

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

	  if(!b_found_mixed_evt) s_used_mx_events.insert(event_pl_runnumber1);
	  b_found_mixed_evt = true;
	  this_mixed_event_pl_runnumber = event_pl_runnumber1;

	  bool b_pos_charge1 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);

      	  //final cuts
      	  if(b_pos_charge0 == b_pos_charge1) continue;

	  double px1 = t_K -> GetLeaf("px")->GetValue(0);
	  double py1 = t_K -> GetLeaf("py")->GetValue(0);
	  double pz1 = t_K -> GetLeaf("pz")->GetValue(0);


	  double d_M1 = d_K_m;
	  double d_M0 = d_K_m;
	  double d_E0 = sqrt((px0*px0+py0*py0+pz0*pz0)+d_M0*d_M0);
	  double d_E1 = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);

	  double d_inv_m = sqrt(d_M0*d_M0
				+d_M1*d_M1
				+2.0*d_E0*d_E1
				-2.0*(px0*px1+py0*py1+pz0*pz1) );

	  double d_pT1 = sqrt(px1*px1+py1*py1);
	  double d_mom1 = sqrt(px1*px1+py1*py1+pz1*pz1);
	  double eta1 = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;

	  double  DCA_r1           = t_K-> GetLeaf("DCA_r")->GetValue(0);
	  double  d_TPCnSigmaKaon1 = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
	  double  d_TOFnSigmaKaon1 = t_K-> GetLeaf("d_TOFnSigmaKaon")->GetValue(0);


	  bool b_PHI = true;
	  double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
	  //	  if(d_dip_angle < 0.04) b_PHI = false;

	  double d_min_DCA     = (DCA_r0 < DCA_r1) ? DCA_r0 : DCA_r1;
	  double d_max_TPCnsig = (fabs(d_TPCnSigmaKaon0) > fabs(d_TPCnSigmaKaon1)) ? d_TPCnSigmaKaon0 : d_TPCnSigmaKaon1;
	  double d_max_TOFnsig = (fabs(d_TOFnSigmaKaon0) > fabs(d_TOFnSigmaKaon1)) ? d_TOFnSigmaKaon0 : d_TOFnSigmaKaon1;
	  int i_max_multi      = (nGoodTracks0 > nGoodTracks1) ? nGoodTracks0 : nGoodTracks1;
	  double d_max_multi   = (double) i_max_multi;

	  // if(
	  //    (fabs(d_TPCnSigmaKaon0) > 1)
	  //    ||(fabs(d_TPCnSigmaKaon1) > 1)
	  //    ) b_PHI = false;

	  if(b_PHI) h_mx_prim_inv_m_PHI    -> Fill(d_inv_m);


	  // "Min Bias"      -> 0-240
	  // "0-10% Central  -> 121-240
	  // "10-30 Central  -> 48-121
	  // ">30% Central   -> 0-48
	  bool a_b_cent0[4];
	  a_b_cent0[0] = ( nGoodTracks0 > 0)   && ( nGoodTracks0 < 240);
	  a_b_cent0[1] = ( nGoodTracks0 > 121) && ( nGoodTracks0 < 240);
	  a_b_cent0[2] = ( nGoodTracks0 > 48)  && ( nGoodTracks0 < 121);
	  a_b_cent0[3] = ( nGoodTracks0 > 0)   && ( nGoodTracks0 < 48);

	  bool a_b_cent1[4];
	  a_b_cent1[0] = ( nGoodTracks1 > 0)   && ( nGoodTracks1 < 240);
	  a_b_cent1[1] = ( nGoodTracks1 > 121) && ( nGoodTracks1 < 240);
	  a_b_cent1[2] = ( nGoodTracks1 > 48)  && ( nGoodTracks1 < 121);
	  a_b_cent1[3] = ( nGoodTracks1 > 0)   && ( nGoodTracks1 < 48);

	  bool a_b_pT0[4];
	  a_b_pT0[0]    = (d_pT0 > 0.5) && (d_pT0 <= 1.0);
	  a_b_pT0[1]    = (d_pT0 > 1.0) && (d_pT0 <= 2.0);
	  a_b_pT0[2]    = (d_pT0 > 2.0) && (d_pT0 <= 3.0);
	  a_b_pT0[3]    = (d_pT0 > 0.5);

	  bool a_b_pT1[4];
	  a_b_pT1[0]    = (d_pT1 > 0.5) && (d_pT1 <= 1.0);
	  a_b_pT1[1]    = (d_pT1 > 1.0) && (d_pT1 <= 2.0);
	  a_b_pT1[2]    = (d_pT1 > 2.0) && (d_pT1 <= 3.0);
	  a_b_pT1[3]    = (d_pT1 > 0.5);



	  //eta fiducial cut
	  bool b_K0_eta   = (eta0 >= -1.47) && (eta0 <= 0.0);
	  bool b_K1_eta   = (eta1 >= -1.47) && (eta1 <= 0.0);

	  // pT cuts
	  // K cuts:
	  //	      if((d_pT1 > 2.0) || (d_pT0 > 2.0) ) continue;
	  // phi cuts:
	  double d_pT_phi = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);//px0*py0*px0*py0 + px1*py1*px1*py1 + 2.*px0*px1  + 2.*py0*py1;
	  double m_phi = 1.019455;
	  double d_mT_phi = sqrt(d_pT_phi*d_pT_phi + m_phi*m_phi );

	  double d_phi_pz = pz0+pz1;
	  double d_phi_E  = d_E0+d_E1;
	  double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;

	  //	  if( (!b_K0_eta) || (!b_K1_eta)) continue;
	  //if((d_pT_phi < 1.0) || (d_pT_phi > 3.0) ) continue;
	  if(d_pT_phi > 3.0) continue;
	  //	      if(d_mT_phi > 2.0) continue;
	  //y rapidity cuts
	  if((d_phi_y<-2.0)||(d_phi_y>0.0)) continue;

	  //only fill if inv m is in range
	  if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;

	  //Ana PID cuts
	  if(b_K_pl_check)
	    {
	       d_max_TPCnsig = d_TPCnSigmaKaon1;
	      d_max_TOFnsig = d_TOFnSigmaKaon1;
	    }//if(b_K_pl_check)


	  for( int k = 0; k < 4; k++) if(b_PHI && a_b_cent0[k] && a_b_cent1[k]) a_h3_mx_prim_inv_nsig[k] -> Fill(d_inv_m,d_max_TOFnsig,d_max_TPCnsig);
	  //make the pT inclusive (one particle in the pair must pass the cuts
	  for( int k = 4; k < 8; k++) if(b_PHI && (a_b_pT0[k-4] || a_b_pT1[k-4])) a_h3_mx_prim_inv_nsig[k] -> Fill(d_inv_m,d_max_TOFnsig,d_max_TPCnsig);


	  if(b_PHI && a_b_cent0[0] && a_b_cent1[0] && (fabs(d_max_TPCnsig) < 2.0) && (fabs(d_max_TOFnsig) < 3.0) ) h_mx_TPC2_TOF3_invM -> Fill(d_inv_m); //test

	  //  if(b_PHI) h3_mx_prim_inv_nsig -> Fill(d_inv_m,d_max_TOFnsig,d_max_TPCnsig);
	  //	  if(b_PHI) h3_mx_prim_inv_nsig -> Fill(d_inv_m,d_max_multi,d_max_TPCnsig);
	  //if(b_PHI) cout<<" "<<event_pl_runnumber0<<"  "<<event_pl_runnumber1<<"  "<<nGoodTracks0<<"  "<<nGoodTracks1<<" "<< i_max_multi <<" "<<d_max_multi <<endl;
	  //	  if((i>10000)||(j>10000)) break;
	}//for(int j = 0; j < N_entries; j++)

    }//  for(int i = 0; i < N_entries; i++)
  h_mx_prim_inv_m_PHI    -> SaveAs("h_mx_prim_inv_m_PHI.root");

  //  h3_mx_prim_inv_nsig -> SaveAs("h3_mx_prim_inv_nsig.root");
  for( int k = 0; k < 8; k++)
    {
      if(k >= 8) break;
      cout<<"saving: "<<k<<endl;
      a_h3_mx_prim_inv_nsig[k] -> SaveAs(Form("h3_mx_prim_inv_nsig_%s.root",TS_ctr[k].Data()));
    }//for( int k = 0; k < 8; k++)

  TFile * tf_out = new TFile("h_K_PID_QA.root","RECREATE");
  h2_TOF_beta_pq                  -> Write();
  h2_dEdx_pq                      -> Write();

  h3_TOF_beta_pq_vs_nsig          -> Write();
  h3_dEdx_pq_vs_nsig              -> Write();

  h2_TOF_beta_TPC2sig_TOF_3sig_pq -> Write();
  h2_dEdx_TPC2sig_TOF_3sig_pq     -> Write();


  h_Kmin_TPC2_TOF3_pT  -> Write();
  h_Kpl_TPC2_TOF3_pT   -> Write();
  h_phi_TPC2_TOF3_pT   -> Write();
  h_phi_TPC2_TOF3_mT   -> Write();

  h_phi_TPC2_TOF3_y       -> Write();

  h_phi_TPC2_TOF3_tight_y -> SetTitle("#phi Uncorrected dN/dy");
  h_phi_TPC2_TOF3_tight_y -> SetXTitle("y");
  h_phi_TPC2_TOF3_tight_y -> SetYTitle("Counts");
  h_phi_TPC2_TOF3_tight_y -> Write();

  h_K_pl_pT  -> Write();
  h_K_min_pT -> Write();

  h_K_pl_eta  -> Write();
  h_K_min_eta -> Write();

  h_dip_angle                     -> Write();

  h_mx_TPC2_TOF3_invM -> SaveAs("h_mx_TPC2_TOF3_invM_out.root");

  new TCanvas();
  h_TPC2_TOF3_invM    -> Draw();
  h_mx_TPC2_TOF3_invM -> SetFillColor(kYellow);
  //  h_mx_TPC2_TOF3_invM -> Draw("sames");

  return;
}//void read_K_tree()

//is only read by compiler
#ifndef __CINT__
int main(int argc, char * argv[]){
  gSystem->Load("libTree");
  TApplication App("Analysis", &argc, argv);

  read_K_tree();
  pthread_exit(NULL);

  return 0;
}
#endif
