#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

// Define Global contstant
const Int_t Ncentrailities = 7;
const Int_t Nrapbins       = 15;

// void AverageFlow(const char* inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/bbcEventPlane/phi_v1_vs_invM_first_run.root")
void AverageFlow(const char* inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/resultSysAna/vertexZcut_9_3BDBD46209544A6CEF9A1AEEB1659D51.picoDst.result.combined.root")
{
  TFile* F = TFile::Open(inFile);
  if (!F || F->IsZombie()){
    std::cout<<"Cannot open the file"<<inFile<<std::endl;
    return;
  }

  TProfile3D *tp3d_proton = (TProfile3D*)F->Get("profile3D_proton_v1");
  TProfile3D *tp3d_pionPlus = (TProfile3D*)F->Get("profile3D_pionPlus_v1");
  TProfile3D *tp3d_pionMinus = (TProfile3D*)F->Get("profile3D_pionMinus_v1");
  TProfile3D *tp3d_kaonPlus = (TProfile3D*)F->Get("profile3D_kaonPlus_v1");
  TProfile3D *tp3d_kaonMinus = (TProfile3D*)F->Get("profile3D_kaonMinus_v1");

  TProfile *profile_correlation_tpc_east_tpc_west_input = (TProfile*)F->Get("profile_correlation_tpc_east_tpc_west");
  TProfile *profile_correlation_tpc_east_bbc_east_input = (TProfile*)F->Get("profile_correlation_tpc_east_bbc_east");
  TProfile *profile_correlation_tpc_west_bbc_east_input = (TProfile*)F->Get("profile_correlation_tpc_west_bbc_east");

  TH1D *h_correlation_tpc_east_tpc_west_input = new TH1D("h_correlation_tpc_east_tpc_west_input","<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>",7,0.5,7+0.5);
  TH1D *h_correlation_tpc_east_bbc_east_input = new TH1D("h_correlation_tpc_east_bbc_east_input","<cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>",7,0.5,7+0.5);
  TH1D *h_correlation_tpc_west_bbc_east_input = new TH1D("h_correlation_tpc_west_bbc_east_input","<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})>",7,0.5,7+0.5);

  for(int i = 1;i<=7;i++)
  {
    h_correlation_tpc_east_tpc_west_input->SetBinContent(i,profile_correlation_tpc_east_tpc_west_input->GetBinContent(i));
    h_correlation_tpc_east_tpc_west_input->SetBinError(i,profile_correlation_tpc_east_tpc_west_input->GetBinError(i));

    h_correlation_tpc_east_bbc_east_input->SetBinContent(i,profile_correlation_tpc_east_bbc_east_input->GetBinContent(i));
    h_correlation_tpc_east_bbc_east_input->SetBinError(i,profile_correlation_tpc_east_bbc_east_input->GetBinError(i));

    h_correlation_tpc_west_bbc_east_input->SetBinContent(i,profile_correlation_tpc_west_bbc_east_input->GetBinContent(i));
    h_correlation_tpc_west_bbc_east_input->SetBinError(i,profile_correlation_tpc_west_bbc_east_input->GetBinError(i));


  }

  TProfile2D *tp2d_proton = tp3d_proton->Project3DProfile("xz");
  TProfile2D *tp2d_pionPlus = tp3d_pionPlus->Project3DProfile("xz");
  TProfile2D *tp2d_pionMinus = tp3d_pionMinus->Project3DProfile("xz");
  TProfile2D *tp2d_kaonPlus = tp3d_kaonPlus->Project3DProfile("xz");
  TProfile2D *tp2d_kaonMinus = tp3d_kaonMinus->Project3DProfile("xz");

  // TProfile *profile_resolution2 = new TProfile("profile_resolution2","<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})>*<cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>/<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>",7,0.5,7+0.5,-1.0,1.0,"");
  // profile_resolution2->GetXaxis()->SetTitle("Centrality (%)");
  // profile_resolution2->GetYaxis()->SetTitle("Resolution Squared");
  //
  // profile_resolution2->Multiply(profile_correlation_tpc_east_bbc_east_input,profile_correlation_tpc_west_bbc_east_input,1,1,"");
  // profile_resolution2->Divide(profile_correlation_tpc_east_tpc_west_input);

  TH1D *h_resolution2 = new TH1D("h_resolution2","#frac{<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})> #times <cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>}{<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>}",7,0.5,7+0.5);
  h_resolution2->GetXaxis()->SetTitle("Centrality (%)");
  h_resolution2->GetYaxis()->SetTitle("Resolution Squared");

  TH1D *h_resolution = new TH1D("h_resolution","#sqrt{#frac{<cos(#psi^{TPC west}_{1} #minus #psi^{BBC east}_{1})> #times <cos(#psi^{TPC east}_{1} #minus #psi^{BBC east}_{1})>}{<cos(#psi^{TPC east}_{1} #minus #psi^{TPC west}_{1})>}}",7,0.5,7+0.5);
  h_resolution->GetXaxis()->SetTitle("Centrality (%)");
  h_resolution->GetYaxis()->SetTitle("Resolution");

  h_resolution2->Multiply(h_correlation_tpc_east_bbc_east_input,h_correlation_tpc_west_bbc_east_input,1,1,"");
  h_resolution2->Divide(h_correlation_tpc_east_tpc_west_input);



  h_resolution2->Draw();
  TH1D * h_PRO_v1_vs_y_10_30 = new TH1D("h_PRO_v1_vs_y_10_30","h_PRO_v1_vs_y_10_30",15,-2.9,0.1);
  TH1D * h_pionPlus_v1_vs_y_10_30 = new TH1D("h_pionPlus_v1_vs_y_10_30","h_pionPlus_v1_vs_y_10_30",15,-2.9,0.1);
  TH1D * h_pionMinus_v1_vs_y_10_30 = new TH1D("h_pionMinus_v1_vs_y_10_30","h_pionMinus_v1_vs_y_10_30",15,-2.9,0.1);
  TH1D * h_kaonPlus_v1_vs_y_10_30 = new TH1D("h_kaonPlus_v1_vs_y_10_30","h_kaonPlus_v1_vs_y_10_30",15,-2.9,0.1);
  TH1D * h_kaonMinus_v1_vs_y_10_30 = new TH1D("h_kaonMinus_v1_vs_y_10_30","h_kaonMinus_v1_vs_y_10_30",15,-2.9,0.1);


  double  value1[Ncentrailities][Nrapbins] = {0};
  double  error1[Ncentrailities][Nrapbins] = {0};
  double  value1_avg[Nrapbins]  = {0};
  double  weight1[Nrapbins]     = {0};
  double  value1_fnal[Nrapbins] = {0};
  double  errorbar1[Nrapbins] = {0};

  double  value2[Ncentrailities][Nrapbins] = {0};
  double  error2[Ncentrailities][Nrapbins] = {0};
  double  value2_avg[Nrapbins]  = {0};
  double  weight2[Nrapbins]     = {0};
  double  value2_fnal[Nrapbins] = {0};
  double  errorbar2[Nrapbins] = {0};

  double  value3[Ncentrailities][Nrapbins] = {0};
  double  error3[Ncentrailities][Nrapbins] = {0};
  double  value3_avg[Nrapbins]  = {0};
  double  weight3[Nrapbins]     = {0};
  double  value3_fnal[Nrapbins] = {0};
  double  errorbar3[Nrapbins] = {0};

  double  value4[Ncentrailities][Nrapbins] = {0};
  double  error4[Ncentrailities][Nrapbins] = {0};
  double  value4_avg[Nrapbins]  = {0};
  double  weight4[Nrapbins]     = {0};
  double  value4_fnal[Nrapbins] = {0};
  double  errorbar4[Nrapbins] = {0};

  double  value5[Ncentrailities][Nrapbins] = {0};
  double  error5[Ncentrailities][Nrapbins] = {0};
  double  value5_avg[Nrapbins]  = {0};
  double  weight5[Nrapbins]     = {0};
  double  value5_fnal[Nrapbins] = {0};
  double  errorbar5[Nrapbins] = {0};

  double  d_resolution2[Ncentrailities] ={0};
  double  d_resolution2_err[Ncentrailities] ={0};

  double  d_resolution[Ncentrailities] ={0};
  double  d_resolution_err[Ncentrailities] ={0};

  for(int i=0;i<7;i++)
  {
    d_resolution2[i] = h_resolution2->GetBinContent(i+1);
    d_resolution2_err[i] = h_resolution2->GetBinError(i+1);
    if(d_resolution2[i]>0)
    {
      d_resolution[i] = TMath::Sqrt(d_resolution2[i]);
      d_resolution_err[i] = 0.5 * d_resolution2_err[i]/TMath::Sqrt(d_resolution2[i]);
      h_resolution->SetBinContent(i+1,d_resolution[i]);
      h_resolution->SetBinError(i+1,d_resolution_err[i]);
    }
    std::cout<< Form("resolution %d = ",i+1) << d_resolution[i]<<std::endl;

  }
  h_resolution->Draw();
  // d_resolution[0] =0.201086 ;
  // d_resolution[1] =0.312114 ; //0.312115
  // d_resolution[2] =0.398419 ;
  // d_resolution[3] =0.454645 ;
  // d_resolution[4] =0.498949 ;
  // d_resolution[5] =0.596083 ; //0.596082
  // d_resolution[6] =0 ;



  // for(int i=0;i<Ncentrailities;i++)
  // {
  //   for(int j=0;j<Nrapbins;j++)
  //   {
  //     if(d_resolution[i+2])
  //     {
  //       if(d_resolution[i+2]>0)
  //     {
  //       double d_v1_raw_pro           = tp2d_proton->GetBinContent(j+1,i+3);
  //       double d_v1_raw_pro_err       = tp2d_proton->GetBinError(j+1,i+3);
  //       double d_v1_raw_pionPlus      = tp2d_pionPlus->GetBinContent(j+1,i+3);
  //       double d_v1_raw_pionPlus_err  = tp2d_pionPlus->GetBinError(j+1,i+3);
  //       double d_v1_raw_pionMinus     = tp2d_pionMinus->GetBinContent(j+1,i+3);
  //       double d_v1_raw_pionMinus_err = tp2d_pionMinus->GetBinError(j+1,i+3);
  //       double d_v1_raw_kaonPlus      = tp2d_kaonPlus->GetBinContent(j+1,i+3);
  //       double d_v1_raw_kaonPlus_err  = tp2d_kaonPlus->GetBinError(j+1,i+3);
  //       double d_v1_raw_kaonMinus     = tp2d_kaonMinus->GetBinContent(j+1,i+3);
  //       double d_v1_raw_kaonMinus_err = tp2d_kaonMinus->GetBinError(j+1,i+3);
  //
  //       value1[i][j] = d_v1_raw_pro/d_resolution[i+2];
  //       error1[i][j] = TMath::Sqrt((d_v1_raw_pro_err/d_resolution[i+2])*(d_v1_raw_pro_err/d_resolution[i+2])
  //                    + (0.5*d_v1_raw_pro*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_pro*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
  //
  //       value2[i][j] = d_v1_raw_pionPlus/d_resolution[i+2];
  //       error2[i][j] = TMath::Sqrt((d_v1_raw_pionPlus_err/d_resolution[i+2])*(d_v1_raw_pionPlus_err/d_resolution[i+2])
  //                    + (0.5*d_v1_raw_pionPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_pionPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
  //
  //       value3[i][j] = d_v1_raw_pionMinus/d_resolution[i+2];
  //       error3[i][j] = TMath::Sqrt((d_v1_raw_pionMinus_err/d_resolution[i+2])*(d_v1_raw_pionMinus_err/d_resolution[i+2])
  //                    + (0.5*d_v1_raw_pionMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_pionMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
  //
  //       value4[i][j] = d_v1_raw_kaonPlus/d_resolution[i+2];
  //       error4[i][j] = TMath::Sqrt((d_v1_raw_kaonPlus_err/d_resolution[i+2])*(d_v1_raw_kaonPlus_err/d_resolution[i+2])
  //                    + (0.5*d_v1_raw_kaonPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_kaonPlus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
  //
  //       value5[i][j] = d_v1_raw_kaonMinus/d_resolution[i+2];
  //       error5[i][j] = TMath::Sqrt((d_v1_raw_kaonMinus_err/d_resolution[i+2])*(d_v1_raw_kaonMinus_err/d_resolution[i+2])
  //                    + (0.5*d_v1_raw_kaonMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3)))*(0.5*d_v1_raw_kaonMinus*d_resolution2_err[i+2]/(pow(d_resolution[i+2],3))));
  //
  //     }
  //    }
  //   }
  // }

  // for(int j=0;j<Nrapbins;j++)
  // {
  //   for(int i=0;i<Ncentrailities;i++)
  //   {
  //     if(error1[i][j] != 0)
  //     {
  //       value1_avg[j] += (value1[i][j]) / ((error1[i][j]) * (error1[i][j]));
  //       weight1[j]    += 1/((error1[i][j]) * (error1[i][j]));
  //     }
  //
  //     if(error2[i][j] != 0)
  //     {
  //       value2_avg[j] += (value2[i][j]) / ((error2[i][j]) * (error2[i][j]));
  //       weight2[j]    += 1/((error2[i][j]) * (error2[i][j]));
  //     }
  //
  //     if(error3[i][j] != 0)
  //     {
  //       value3_avg[j] += (value3[i][j]) / ((error3[i][j]) * (error3[i][j]));
  //       weight3[j]    += 1/((error3[i][j]) * (error3[i][j]));
  //     }
  //
  //     if(error4[i][j] != 0)
  //     {
  //       value4_avg[j] += (value4[i][j]) / ((error4[i][j]) * (error4[i][j]));
  //       weight4[j]    += 1/((error4[i][j]) * (error4[i][j]));
  //     }
  //
  //     if(error5[i][j] != 0)
  //     {
  //       value5_avg[j] += (value5[i][j]) / ((error5[i][j]) * (error5[i][j]));
  //       weight5[j]    += 1/((error5[i][j]) * (error5[i][j]));
  //     }
  //
  //   }
  //   if(weight1[j]!=0)
  //   {
  //     value1_fnal[j] = value1_avg[j]/weight1[j];
  //     errorbar1[j] = TMath::Sqrt(1/weight1[j]);
  //   }
  //
  //   if(weight2[j]!=0)
  //   {
  //     value2_fnal[j] = value2_avg[j]/weight2[j];
  //     errorbar2[j] = TMath::Sqrt(1/weight2[j]);
  //   }
  //
  //   if(weight3[j]!=0)
  //   {
  //     value3_fnal[j] = value3_avg[j]/weight3[j];
  //     errorbar3[j] = TMath::Sqrt(1/weight3[j]);
  //   }
  //
  //   if(weight4[j]!=0)
  //   {
  //     value4_fnal[j] = value4_avg[j]/weight4[j];
  //     errorbar4[j] = TMath::Sqrt(1/weight4[j]);
  //   }
  //
  //   if(weight5[j]!=0)
  //   {
  //     value5_fnal[j] = value5_avg[j]/weight5[j];
  //     errorbar5[j] = TMath::Sqrt(1/weight5[j]);
  //   }
  //
  // }

  // for(int j=0;j<Nrapbins;j++)
  // {
  //   h_PRO_v1_vs_y_10_30->SetBinContent((j+1),value1_fnal[j]);
  //   h_PRO_v1_vs_y_10_30->SetBinError((j+1),errorbar1[j]);
  //
  //   h_pionPlus_v1_vs_y_10_30->SetBinContent((j+1),value2_fnal[j]);
  //   h_pionPlus_v1_vs_y_10_30->SetBinError((j+1),errorbar2[j]);
  //
  //   h_pionMinus_v1_vs_y_10_30->SetBinContent((j+1),value3_fnal[j]);
  //   h_pionMinus_v1_vs_y_10_30->SetBinError((j+1),errorbar3[j]);
  //
  //   h_kaonPlus_v1_vs_y_10_30->SetBinContent((j+1),value4_fnal[j]);
  //   h_kaonPlus_v1_vs_y_10_30->SetBinError((j+1),errorbar4[j]);
  //
  //   h_kaonMinus_v1_vs_y_10_30->SetBinContent((j+1),value5_fnal[j]);
  //   h_kaonMinus_v1_vs_y_10_30->SetBinError((j+1),errorbar5[j]);
  //
  //
  //
  // }
  TFile* o1 = new TFile("averageflow_vertexZcut_9.root","RECREATE");
  // h_PRO_v1_vs_y_10_30->Write();
  // h_pionPlus_v1_vs_y_10_30->Write();
  // h_pionMinus_v1_vs_y_10_30->Write();
  // h_kaonPlus_v1_vs_y_10_30->Write();
  // h_kaonMinus_v1_vs_y_10_30->Write();
  // profile_resolution2->Write();
  profile_correlation_tpc_east_tpc_west_input->Write();
  profile_correlation_tpc_east_bbc_east_input->Write();
  profile_correlation_tpc_west_bbc_east_input->Write();
  h_correlation_tpc_east_tpc_west_input->Write();
  h_correlation_tpc_east_bbc_east_input->Write();
  h_correlation_tpc_west_bbc_east_input->Write();
  h_resolution2->Write();
  h_resolution->Write();


  return;

}
