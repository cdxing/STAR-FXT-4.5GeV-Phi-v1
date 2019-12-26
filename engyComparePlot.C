#include <TStyle.h>


void engyComparePlot()
{
  // TFile *fsyserr1 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_1.root","READ");
  // TFile *fsyserr2 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_2.root","READ");
  // TFile *fsyserr3 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_3.root","READ");
  // TFile *fsyserr4 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_4.root","READ");
  // TFile *fsyserr5 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_5.root","READ");
  // TFile *fsyserr6 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_6.root","READ");
  // TFile *fsyserr7 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_7.root","READ");
  // TFile *fsyserr8 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_8.root","READ");
  // TFile *fsyserr9 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_9.root","READ");
  // TFile *fsyserr10 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_10.root","READ");
  // TFile *fsyserr11 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_11.root","READ");
  // TFile *fsyserr12 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_12.root","READ");
  // TFile *fsyserr13 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_13.root","READ");
  // TFile *fsyserr14 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_merged_14.root","READ");
  // TFile *fsyserr15 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_15.root","READ");
  // TFile *fsyserr16 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_16.root","READ");
  // TFile *fsyserr17 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_17.root","READ");
  // TFile *fsyserr18 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_18.root","READ");
  // TFile *fsyserr19 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_19.root","READ");
  // TFile *fsyserr20 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_20.root","READ");
  // TFile *fsyserr21 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_hadd_21.root","READ");
  // TFile *fsyserr22 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_plus_bg_polN.root","READ");
  // TFile *fsyserr23 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_eff_corr_EP.root","READ");
  // TFile *fsyserr24 = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result21cuts_v2/sys_err_21cuts_eff_corr_invM.root","READ");
  //
  //   double slopeErr[24];
  //   slopeErr[0] = ((TH1D*) fsyserr1 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[1] = ((TH1D*) fsyserr2 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[2] = ((TH1D*) fsyserr3 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[3] = ((TH1D*) fsyserr4 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[4] = ((TH1D*) fsyserr5 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[5] = ((TH1D*) fsyserr6 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[6] = ((TH1D*) fsyserr7 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[7] = ((TH1D*) fsyserr8 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[8] = ((TH1D*) fsyserr9 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[9] = ((TH1D*) fsyserr10 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[10] = ((TH1D*) fsyserr11 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[11] = ((TH1D*) fsyserr12 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[12] = ((TH1D*) fsyserr13 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[13] = ((TH1D*) fsyserr14 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[14] = ((TH1D*) fsyserr15 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[15] = ((TH1D*) fsyserr16 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[16] = ((TH1D*) fsyserr17 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[17] = ((TH1D*) fsyserr18 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[18] = ((TH1D*) fsyserr19 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[19] = ((TH1D*) fsyserr20 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[20] = ((TH1D*) fsyserr21 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[21] = ((TH1D*) fsyserr22 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[22] = ((TH1D*) fsyserr23 ->Get("h_slope"))->GetStdDev(1);
  //   slopeErr[23] = ((TH1D*) fsyserr24 ->Get("h_slope"))->GetStdDev(1);

    TCanvas *c1 = new TCanvas("c1","c1",200,10,1024,768);
    c1->DrawFrame(3., -0.05, 250., 0.05);

    //4.5 GeV phi
    double x0[1]    = {4.5};
    double zero0[1] = {0};
    double py0[1]      = {-0.01301};
    double ey_stat0[1] = {0.01566};
    ey_stat0[0] *= sqrt(2);
    double ey_sys0[1]  = {0.013012};
    // ey_sys0[0]=sqrt(
    //   slopeErr[0]**2+
    //   slopeErr[1]**2+
    //   slopeErr[2]**2+
    //   slopeErr[3]**2+
    //   slopeErr[4]**2+
    //   slopeErr[5]**2+
    //   slopeErr[6]**2+
    //   slopeErr[7]**2+
    //   slopeErr[8]**2+
    //   slopeErr[9]**2+
    //   slopeErr[10]**2+
    //   slopeErr[11]**2+
    //   slopeErr[12]**2+
    //   slopeErr[13]**2+
    //   slopeErr[14]**2+
    //   slopeErr[15]**2+
    //   slopeErr[16]**2+
    //   slopeErr[17]**2+
    //   slopeErr[18]**2+
    //   slopeErr[19]**2+
    //   slopeErr[20]**2+
    //   slopeErr[21]**2+
    //   slopeErr[22]**2+
    //   slopeErr[21]**2
    // );
    ey_sys0[0] *= sqrt(2);
    std::cout<< "ey_sys0[0] = "<< ey_sys0[0] << std::endl;
    //0.013012
    // BES-I phi
    double x1[6]    = { 11.5, 14.5, 19.6, 27, 39, 200};
    double zero1[6] = {0};
    double py1[6]      = {0.0199, -0.0311, -0.0304, -0.02,   -0.013,   -0.00366944};
    double ey_stat1[6] = {0.0183, 0.00938, 0.00763, 0.00644, 0.006142, 0.00203012};
    double ey_sys1[6]  = {0.03,   0.021,   0.014,   0.015,   0.012,    0.000533};

    // BES-I Proton
    double x2[8]    = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero2[8] = {0};
    double py2[8]      = {0.022534   , -0.000486069, -0.00684062, -0.00690107, -0.00561899, -0.00415636, -0.00127671, 0.000601299};
    double ey_stat2[8] = {0.000863638, 0.000618966 , 0.000843574, 0.000853584, 0.000836661, 0.000565664, 0.00168219 , 0.000738392};
    double ey_sys2[8]  = {0.00447883 , 0.00225644  , 0.0021     , 0.000933001, 0.0008     , 0.0001     , 0.0031039  , 0.00011};

    // BES-I KPlus
    double x3[8]       = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero3[8]    = {0};
    double py3[8]      = {-0.0137338, -0.0135527,  -0.0127973, -0.0124038,  -0.0104937,  -0.0082035,  -0.00771108, -0.00270561};
    double ey_stat3[8] = {0.00166579, 0.00100918,  0.001172,   0.00103033,  0.000927541, 0.000581183, 0.00185789,  0.000290576};
    double ey_sys3[8]  = {0.00076095, 0.000472041, 0.00104204, 0.000836145, 0.00021,     0.0002,      0.0029149,   0.000137};

    // BES-I KMinus
    double x4[8]       = {7.7, 11.5, 14.5, 19.6, 27, 39, 62.4, 200};
    double zero4[8]    = {0};
    double py4[8]      = {-0.00756657, -0.0163133,  -0.0160707,  -0.0174652,  -0.0147908, -0.011649,   -0.010862,  -0.00336346};
    double ey_stat4[8] = {0.00290359,  0.00150192,  0.00162578,  0.00132712,  0.00112134, 0.000672674, 0.00205018, 0.000302512};
    double ey_sys4[8]  = {0.000950781, 0.000974057, 0.000974057, 0.000458308, 0.0002,     0.0001,      0.0030239,  0.000194};

    // data set (1) with stat and sys errors


    // Now draw data set (1)
    // We first have to draw it only with the stat errors
    gStyle->SetOptDate(0);
    gStyle->SetEndErrorSize(6);
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogx();

    //4.5 GeV phi
    TGraphErrors *graph0 = new TGraphErrors(1, x0, py0, zero0, ey_stat0);
    graph0->SetMarkerStyle(28);
    graph0->SetMarkerColor(kRed);
    graph0->SetLineColor(kRed);
    graph0->SetMarkerSize(2);
    graph0->Draw("P");

    // BES-I phi
    TGraphErrors *graph1 = new TGraphErrors(6, x1, py1, zero1, ey_stat1);
    graph1->SetMarkerStyle(34);
    graph1->SetMarkerColor(kGreen+3);
    graph1->SetLineColor(kGreen+3);
    graph1->SetMarkerSize(2);
    graph1->Draw("P");

    // BES-I Proton
    TGraphErrors *graph2 = new TGraphErrors(8, x2, py2, zero2, ey_stat2);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerColor(kBlack);
    graph2->SetLineColor(kBlack);
    graph2->SetMarkerSize(2);
    graph2->Draw("P");

    // BES-I KPlus
    TGraphErrors *graph3 = new TGraphErrors(8, x3, py3, zero3, ey_stat3);
    graph3->SetMarkerStyle(22);
    graph3->SetMarkerColor(kBlue);
    graph3->SetLineColor(kBlue);
    graph3->SetMarkerSize(2);
    graph3->Draw("P");

    // BES-I KMinus
    TGraphErrors *graph4 = new TGraphErrors(8, x4, py4, zero4, ey_stat4);
    graph4->SetMarkerStyle(23);
    graph4->SetMarkerColor(6);
    graph4->SetLineColor(6);
    graph4->SetMarkerSize(2);
    graph4->Draw("P");

    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(graph0,"#phi @ 4.5 GeV","p");
    legend->AddEntry(graph1,"#phi @ BES-I","p");
    legend->AddEntry(graph2,"p @ BES-I","p");
    legend->AddEntry(graph3,"K^{+} @ BES-I","p");
    legend->AddEntry(graph4,"K^{-} @ BES-I","p");

    legend->Draw();
    // Now we have to somehow depict the sys errors

    //4.5 GeV phi
    TGraphErrors *graph0_sys = new TGraphErrors(1, x0, py0, zero0, ey_sys0);
    graph0_sys->SetLineColor(kRed);
    graph0_sys->Draw("[]");

    cout << "py0 = "     << py0[0] << endl;
    cout << "ey_stat0 = "<< ey_stat0[0] << endl;
    cout << "ey_sys0 = " << ey_sys0[0] << endl;

    // BES-I phi
    TGraphErrors *graph1_sys = new TGraphErrors(6, x1, py1, zero1, ey_sys1);
    graph1_sys->SetLineColor(kGreen+3);
    graph1_sys->Draw("[]");

    // BES-I Proton
    TGraphErrors *graph2_sys = new TGraphErrors(8, x2, py2, zero2, ey_sys2);
    graph2_sys->SetLineColor(kBlack);
    graph2_sys->Draw("[]");

    // BES-I KPlus
    TGraphErrors *graph3_sys = new TGraphErrors(8, x3, py3, zero3, ey_sys3);
    graph3_sys->SetLineColor(kBlue);
    graph3_sys->Draw("[]");

    // BES-I KMinus
    TGraphErrors *graph3_sys = new TGraphErrors(8, x4, py4, zero4, ey_sys4);
    graph3_sys->SetLineColor(6);
    graph3_sys->Draw("[]");

    TLine *line = new TLine(3, 0, 250, 0);
    line->SetLineStyle(7);
    line->Draw();
}
