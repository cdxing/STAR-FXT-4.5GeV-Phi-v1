#include <TStyle.h>


void sysErrBar6DataPlot()
{
  TFile *f1 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/resultCutVariation/Phi_flow_fitting_out_primary.root","READ");
  std::cout<< "test 5" << std::endl;

  TFile *fsyserr1 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_1.root","READ");
  TFile *fsyserr2 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_2.root","READ");
  TFile *fsyserr3 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_3.root","READ");
  TFile *fsyserr4 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_4.root","READ");
  TFile *fsyserr5 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_5.root","READ");
  TFile *fsyserr6 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_6.root","READ");
  TFile *fsyserr7 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_7.root","READ");
  TFile *fsyserr8 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_8.root","READ");
  TFile *fsyserr9 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_9.root","READ");
  TFile *fsyserr10 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_10.root","READ");
  TFile *fsyserr11 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_11.root","READ");
  TFile *fsyserr12 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_12.root","READ");
  TFile *fsyserr13 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_13.root","READ");
  TFile *fsyserr14 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_14.root","READ");
  TFile *fsyserr15 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_15.root","READ");
  TFile *fsyserr16 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_16.root","READ");
  TFile *fsyserr17 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_17.root","READ");
  TFile *fsyserr18 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_18.root","READ");
  TFile *fsyserr19 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_19.root","READ");
  TFile *fsyserr20 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_20.root","READ");
  TFile *fsyserr21 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_21.root","READ");
  TFile *fsyserr22 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_plus_bg_polN.root","READ");
  TFile *fsyserr23 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_eff_corr_EP.root","READ");
  TFile *fsyserr24 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_eff_corr_invM.root","READ");

  TH1D  *h1 = (TH1D*) f1->Get("h_phi_v1_vs_y_flip");

  double bin4err[24]={};
  double bin5err[24]={};
  double bin6err[24]={};
  std::cout<< "test 4" << std::endl;

  bin4err[0] = ((TH1D*) fsyserr1 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[1] = ((TH1D*) fsyserr2 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[2] = ((TH1D*) fsyserr3 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[3] = ((TH1D*) fsyserr4 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[4] = ((TH1D*) fsyserr5 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[5] = ((TH1D*) fsyserr6 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[6] = ((TH1D*) fsyserr7 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[7] = ((TH1D*) fsyserr8 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[8] = ((TH1D*) fsyserr9 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[9] = ((TH1D*) fsyserr10 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[10] = ((TH1D*) fsyserr11 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[11] = ((TH1D*) fsyserr12 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[12] = ((TH1D*) fsyserr13 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[13] = ((TH1D*) fsyserr14 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[14] = ((TH1D*) fsyserr15 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[15] = ((TH1D*) fsyserr16 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[16] = ((TH1D*) fsyserr17 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[17] = ((TH1D*) fsyserr18 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[18] = ((TH1D*) fsyserr19 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[19] = ((TH1D*) fsyserr20 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[20] = ((TH1D*) fsyserr21 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[21] = ((TH1D*) fsyserr22 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[22] = ((TH1D*) fsyserr23 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  bin4err[23] = ((TH1D*) fsyserr24 ->Get("h_bin4"))->GetStdDev(1); delete h_bin4;
  std::cout<< "test 3.2" << std::endl;

  bin5err[0] = ((TH1D*) fsyserr1 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[1] = ((TH1D*) fsyserr2 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[2] = ((TH1D*) fsyserr3 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[3] = ((TH1D*) fsyserr4 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[4] = ((TH1D*) fsyserr5 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[5] = ((TH1D*) fsyserr6 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[6] = ((TH1D*) fsyserr7 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[7] = ((TH1D*) fsyserr8 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[8] = ((TH1D*) fsyserr9 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[9] = ((TH1D*) fsyserr10 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[10] = ((TH1D*) fsyserr11 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[11] = ((TH1D*) fsyserr12 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[12] = ((TH1D*) fsyserr13 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[13] = ((TH1D*) fsyserr14 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[14] = ((TH1D*) fsyserr15 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[15] = ((TH1D*) fsyserr16 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[16] = ((TH1D*) fsyserr17 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[17] = ((TH1D*) fsyserr18 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[18] = ((TH1D*) fsyserr19 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[19] = ((TH1D*) fsyserr20 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[20] = ((TH1D*) fsyserr21 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[21] = ((TH1D*) fsyserr22 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[22] = ((TH1D*) fsyserr23 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  bin5err[23] = ((TH1D*) fsyserr24 ->Get("h_bin5"))->GetStdDev(1); delete h_bin5;
  std::cout<< "test 3.1" << std::endl;

  bin6err[0] = ((TH1D*) fsyserr1 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "test 3.04" << std::endl;

  bin6err[1] = ((TH1D*) fsyserr2 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "bin6err[1] = "<<  bin6err[1]<< std::endl;

  std::cout<< "test 3.031" << std::endl;

  bin6err[2] = ((TH1D*) fsyserr3 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "test 3.03" << std::endl;

  bin6err[3] = ((TH1D*) fsyserr4 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[4] = ((TH1D*) fsyserr5 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[5] = ((TH1D*) fsyserr6 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "test 3.02" << std::endl;

  bin6err[6] = ((TH1D*) fsyserr7 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[7] = ((TH1D*) fsyserr8 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[8] = ((TH1D*) fsyserr9 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[9] = ((TH1D*) fsyserr10 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[10] = ((TH1D*) fsyserr11 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[11] = ((TH1D*) fsyserr12 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[12] = ((TH1D*) fsyserr13 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "test 3.01" << std::endl;

  bin6err[13] = ((TH1D*) fsyserr14 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[14] = ((TH1D*) fsyserr15 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[15] = ((TH1D*) fsyserr16 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[16] = ((TH1D*) fsyserr17 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[17] = ((TH1D*) fsyserr18 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[18] = ((TH1D*) fsyserr19 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[19] = ((TH1D*) fsyserr20 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[20] = ((TH1D*) fsyserr21 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[21] = ((TH1D*) fsyserr22 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[22] = ((TH1D*) fsyserr23 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  bin6err[23] = ((TH1D*) fsyserr24 ->Get("h_bin6"))->GetStdDev(1); delete h_bin6;
  std::cout<< "test 3" << std::endl;

  TF1   *tf1 = new TF1("tf1","[0] + [1]*x ",-1.53,1.53);
  h1->Fit(tf1,"E","R",-1.53,1.53);
    TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
    c1->DrawFrame(-1.53, -0.3, 1.53, 0.3);
    double x[6]    = {1, 2, 3, 4, 5, 6};
    double zero[6] = {0, 0, 0, 0, 0, 0};
    x[0] = h1->GetBinCenter(26);
    x[1] = h1->GetBinCenter(76);
    x[2] = h1->GetBinCenter(126);
    x[3] = h1->GetBinCenter(180);
    x[4] = h1->GetBinCenter(230);
    x[5] = h1->GetBinCenter(280);

    // data set (1) with stat and sys errors
    double py1[6]      = {1.2, 1.15, 1.19, 0.9, 1.4, 0};
    py1[0] = h1->GetBinContent(26);
    py1[1] = h1->GetBinContent(76);
    py1[2] = h1->GetBinContent(126);
    py1[3] = h1->GetBinContent(180);
    py1[4] = h1->GetBinContent(230);
    py1[5] = h1->GetBinContent(280);

    double ey_stat1[6] = {0.2, 0.18, 0.17, 0.2, 0.4, 0};
    ey_stat1[0] = h1->GetBinError(26);
    ey_stat1[1] = h1->GetBinError(76);
    ey_stat1[2] = h1->GetBinError(126);
    ey_stat1[3] = h1->GetBinError(180);
    ey_stat1[4] = h1->GetBinError(230);
    ey_stat1[5] = h1->GetBinError(280);
    std::cout<< "test 2" << std::endl;

    double ey_sys1[6]  = {0.5, 0.71, 0.76, 0.5, 0.45, 0};
    ey_sys1[0]=sqrt(
      bin6err[0]**2+
      bin6err[1]**2+
      bin6err[2]**2+
      bin6err[3]**2+
      bin6err[4]**2+
      bin6err[5]**2+
      bin6err[6]**2+
      bin6err[7]**2+
      bin6err[8]**2+
      bin6err[9]**2+
      bin6err[10]**2+
      bin6err[11]**2+
      bin6err[12]**2+
      bin6err[13]**2+
      bin6err[14]**2+
      bin6err[15]**2+
      bin6err[16]**2+
      bin6err[17]**2+
      bin6err[18]**2+
      bin6err[19]**2+
      bin6err[20]**2+
      bin6err[21]**2+
      bin6err[22]**2+
      bin6err[23]**2
    );
    // ey_sys1[0]=sqrt(
    //   0.01878**2+
    //   0.02005**2+
    //   0.00417**2+
    //   0.00009302**2+
    //   0.002504**2+
    //   0.000000001317**2+
    //   6]**2+
    //   7]**2+
    //   8]**2+
    //   9]**2+
    //   10]**2+
    //   11]**2+
    //   12]**2+
    //   13]**2+
    //   14]**2+
    //   15]**2+
    //   16]**2+
    //   17]**2+
    //   18]**2+
    //   19]**2+
    //   20]**2+
    //   21]**2+
    //   22]**2+
    //   bin6err[23]**2
    // );

    ey_sys1[1]=sqrt(
      bin5err[0]**2+
      bin5err[1]**2+
      bin5err[2]**2+
      bin5err[3]**2+
      bin5err[4]**2+
      bin5err[5]**2+
      bin5err[6]**2+
      bin5err[7]**2+
      bin5err[8]**2+
      bin5err[9]**2+
      bin5err[10]**2+
      bin5err[11]**2+
      bin5err[12]**2+
      bin5err[13]**2+
      bin5err[14]**2+
      bin5err[15]**2+
      bin5err[16]**2+
      bin5err[17]**2+
      bin5err[18]**2+
      bin5err[19]**2+
      bin5err[20]**2+
      bin5err[21]**2+
      bin5err[22]**2+
      bin5err[23]**2
    );
    ey_sys1[2]=sqrt(
      bin4err[0]**2+
      bin4err[1]**2+
      bin4err[2]**2+
      bin4err[3]**2+
      bin4err[4]**2+
      bin4err[5]**2+
      bin4err[6]**2+
      bin4err[7]**2+
      bin4err[8]**2+
      bin4err[9]**2+
      bin4err[10]**2+
      bin4err[11]**2+
      bin4err[12]**2+
      bin4err[13]**2+
      bin4err[14]**2+
      bin4err[15]**2+
      bin4err[16]**2+
      bin4err[17]**2+
      bin4err[18]**2+
      bin4err[19]**2+
      bin4err[20]**2+
      bin4err[21]**2+
      bin4err[22]**2+
      bin4err[23]**2
    );
    ey_sys1[3]=ey_sys1[2];
    ey_sys1[4]=ey_sys1[1];
    ey_sys1[5]=ey_sys1[0];

    std::cout<< "test 1" << std::endl;

    // Now draw data set (1)
    // We first have to draw it only with the stat errors
    gStyle->SetOptDate(0);

    TGraphErrors *graph1 = new TGraphErrors(6, x, py1, zero, ey_stat1);
    graph1->SetMarkerStyle(20);
    graph1->Draw("PZ");
    tf1->Draw("same");
    // Now we have to somehow depict the sys errors
    TGraphErrors *graph1_sys = new TGraphErrors(6, x, py1, zero, ey_sys1);
    graph1_sys->Draw("[]");
}
