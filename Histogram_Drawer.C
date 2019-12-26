#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TProfile.h"

void Histogram_Drawer(const char* inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/Event_Plane/v1_BBCEP_TOF_cut.root")
{
  TFile* F = TFile::Open("/star/data01/pwg/dchen/Ana/fxtPicoAna/files/Event_Plane/v1_BBCEP_TOF_cut.root");
  if (!F || F->IsZombie()){
    std::cout<<"Cannot open the file"<<inFile<<std::endl;
    return;
  }
  TH2D* H2_1 = F->Get("h2_PRO_plus_v1_vs_y_sub1");
  TH2D* H2_2 = F->Get("h2_PRO_plus_v1_vs_y_sub2");
  TH2D* H2_3 = F->Get("h2_PRO_plus_v1_vs_y");

  TH2D* H2_4 = F->Get("h2_PRO_plus_v1_vs_pT_sub1");
  TH2D* H2_5 = F->Get("h2_PRO_plus_v1_vs_pT_sub2");
  TH2D* H2_6 = F->Get("h2_PRO_plus_v1_vs_pT");

  H2_1->RebinX(10);
  H2_2->RebinX(10);
  H2_3->RebinX(10);
  H2_4->RebinX(10);
  H2_5->RebinX(10);
  H2_6->RebinX(10);


  TProfile* H2_1_pfx = H2_1->ProfileX("H2_1_pfx", 1, -1, "s");
  TProfile* H2_2_pfx = H2_1->ProfileX("H2_2_pfx", 1, -1, "s");
  TProfile* H2_3_pfx = H2_1->ProfileX("H2_3_pfx", 1, -1, "s");

  TProfile* H2_4_pfx = H2_1->ProfileX("H2_4_pfx", 1, -1, "s");
  TProfile* H2_5_pfx = H2_1->ProfileX("H2_5_pfx", 1, -1, "s");
  TProfile* H2_6_pfx = H2_1->ProfileX("H2_6_pfx", 1, -1, "s");

  H2_1_pfx->SetTitle("h2_PRO_plus_v1_vs_y_sub1");
  H2_1_pfx->GetXaxis()->SetTitle("v1");
  H2_1_pfx->GetYaxis()->SetTitle("y");

  H2_2_pfx->SetTitle("h2_PRO_plus_v1_vs_y_sub2");
  H2_2_pfx->GetXaxis()->SetTitle("v1");
  H2_2_pfx->GetYaxis()->SetTitle("y");

  H2_3_pfx->SetTitle("h2_PRO_plus_v1_vs_y");
  H2_3_pfx->GetXaxis()->SetTitle("v1");
  H2_3_pfx->GetYaxis()->SetTitle("y");

  H2_4_pfx->SetTitle("h2_PRO_plus_v1_vs_pT_sub1");
  H2_4_pfx->GetXaxis()->SetTitle("v1");
  H2_4_pfx->GetYaxis()->SetTitle("pT");

  H2_5_pfx->SetTitle("h2_PRO_plus_v1_vs_pT_sub2");
  H2_5_pfx->GetXaxis()->SetTitle("v1");
  H2_5_pfx->GetYaxis()->SetTitle("pT");

  H2_6_pfx->SetTitle("h2_PRO_plus_v1_vs_pT");
  H2_6_pfx->GetXaxis()->SetTitle("v1");
  H2_6_pfx->GetYaxis()->SetTitle("pT");

  TString outFile = "flowXProject";
  TFile *outputFile = new TFile(outFile,"recreate");
  outFile.Append(".picoDst.result.root");
  outputFile->cd();

  H2_1_pfx->Write();
  H2_2_pfx->Write();
  H2_3_pfx->Write();

  H2_4_pfx->Write();
  H2_5_pfx->Write();
  H2_6_pfx->Write();

  H2_1->Write();
  H2_2->Write();
  H2_3->Write();

  H2_4->Write();
  H2_5->Write();
  H2_6->Write();

  // TH2D* H2_1 = F->Get("h2_m2_QA_pq");
  // TH2D* H2_2 = F->Get("h2_m2_QA_pT");
  //
  // TH2D* H2_3 = F->Get("h2_E_plus_pT_vs_y");
  // TH2D* H2_4 = (TH2D*)F->Get("h2_K_plus_pT_vs_y");
  // TH2D* H2_5 = (TH2D*)F->Get("h2_PI_plus_pT_vs_y");
  // TH2D* H2_6 = (TH2D*)F->Get("h2_PRO_plus_pT_vs_y");
  //
  // TH2D* H2_7 = (TH2D*)F->Get("h2_E_minus_pT_vs_y");
  // TH2D* H2_8 = (TH2D*)F->Get("h2_K_minus_pT_vs_y");
  // TH2D* H2_9 = (TH2D*)F->Get("h2_PI_minus_pT_vs_y");
  // TH2D* H2_10 = (TH2D*)F->Get("h2_PRO_minus_pT_vs_y");
  //
  // TH1D* H1_1 = (TH1D*)F->Get("h_E_plus_mT_Diff");
  // TH1D* H1_2 = (TH1D*)F->Get("h_K_plus_mT_Diff");
  // TH1D* H1_3 = (TH1D*)F->Get("h_PI_plus_mT_Diff");
  // TH1D* H1_4 = (TH1D*)F->Get("h_PRO_plus_mT_Diff");
  //
  // TH1D* H1_5 = (TH1D*)F->Get("h_E_minus_mT_Diff");
  // TH1D* H1_6 = (TH1D*)F->Get("h_K_minus_mT_Diff");
  // TH1D* H1_7 = (TH1D*)F->Get("h_PI_minus_mT_Diff");
  // TH1D* H1_8 = (TH1D*)F->Get("h_PRO_minus_mT_Diff");
  //
  // TH1D* H1_9 = (TH1D*)F->Get("h_E_minus_pT");
  // TH1D* H1_10 = (TH1D*)F->Get("h_K_minus_pT");
  // TH1D* H1_11 = (TH1D*)F->Get("h_PI_minus_pT");
  // TH1D* H1_12 = (TH1D*)F->Get("h_PRO_minus_pT");
  //
  // TH1D* H1_13 = (TH1D*)F->Get("h_E_plus_pT");
  // TH1D* H1_14 = (TH1D*)F->Get("h_K_plus_pT");
  // TH1D* H1_15 = (TH1D*)F->Get("h_PI_plus_pT");
  // TH1D* H1_16 = (TH1D*)F->Get("h_PRO_plus_pT");
  //
  // TH1D* H1_17 = (TH1D*)F->Get("h_E_plus_y");
  // TH1D* H1_18 = (TH1D*)F->Get("h_K_plus_y");
  // TH1D* H1_19 = (TH1D*)F->Get("h_PI_plus_y");
  // TH1D* H1_20 = (TH1D*)F->Get("h_PRO_plus_y");
  //
  // TH1D* H1_21 = (TH1D*)F->Get("h_E_minus_y");
  // TH1D* H1_22 = (TH1D*)F->Get("h_K_minus_y");
  // TH1D* H1_23 = (TH1D*)F->Get("h_PI_minus_y");
  // TH1D* H1_24 = (TH1D*)F->Get("h_PRO_minus_y");
  //
  // TCanvas* c1 = new TCanvas("c1", "Canvas_1", 1600, 800);
  // TCanvas* c2 = new TCanvas("c2", "Canvas_2", 1600, 800);
  // TCanvas* c3 = new TCanvas("c3", "Canvas_3", 1600, 800);
  // TCanvas* c4 = new TCanvas("c4", "Canvas_4", 1600, 800);
  //
  // TCanvas* c5 = new TCanvas("c5", "Canvas_5", 1600, 800);
  // TCanvas* c6 = new TCanvas("c6", "Canvas_6", 1600, 800);
  // TCanvas* c7 = new TCanvas("c7", "Canvas_7", 1600, 800);
  // TCanvas* c8 = new TCanvas("c8", "Canvas_8", 1600, 800);
  //
  // TCanvas* c9 = new TCanvas("c9", "Canvas_9", 1600, 800);
  // TCanvas* c10 = new TCanvas("c10", "Canvas_10", 1600, 800);
  // TCanvas* c11 = new TCanvas("c11", "Canvas_11", 1600, 800);
  // TCanvas* c12 = new TCanvas("c12", "Canvas_12", 1600, 800);
  //
  // TCanvas* c13 = new TCanvas("c13", "Canvas_13", 1600, 800);
  // TCanvas* c14 = new TCanvas("c14", "Canvas_14", 1600, 800);
  // TCanvas* c15 = new TCanvas("c15", "Canvas_15", 1600, 800);
  // TCanvas* c16 = new TCanvas("c16", "Canvas_16", 1600, 800);
  //
  // TCanvas* c17 = new TCanvas("c17", "Canvas_17", 1600, 800);
  // TCanvas* c18 = new TCanvas("c18", "Canvas_18", 1600, 800);
  // TCanvas* c19 = new TCanvas("c19", "Canvas_19", 1600, 800);
  // TCanvas* c20 = new TCanvas("c20", "Canvas_20", 1600, 800);
  //
  // TCanvas* c21 = new TCanvas("c21", "Canvas_21", 1600, 800);
  // TCanvas* c22 = new TCanvas("c22", "Canvas_22", 1600, 800);
  // TCanvas* c23 = new TCanvas("c23", "Canvas_23", 1600, 800);
  // TCanvas* c24 = new TCanvas("c24", "Canvas_24", 1600, 800);
  //
  // TCanvas* c25 = new TCanvas("c25", "Canvas_25", 1600, 800);
  // TCanvas* c26 = new TCanvas("c26", "Canvas_26", 1600, 800);
  // TCanvas* c27 = new TCanvas("c27", "Canvas_27", 1600, 800);
  // TCanvas* c28 = new TCanvas("c28", "Canvas_28", 1600, 800);
  //
  // TCanvas* c29 = new TCanvas("c29", "Canvas_29", 1600, 800);
  // TCanvas* c30 = new TCanvas("c30", "Canvas_30", 1600, 800);
  // TCanvas* c31 = new TCanvas("c31", "Canvas_31", 1600, 800);
  // TCanvas* c32 = new TCanvas("c32", "Canvas_32", 1600, 800);
  //
  // TCanvas* c33 = new TCanvas("c33", "Canvas_33", 1600, 800);
  // TCanvas* c34 = new TCanvas("c34", "Canvas_34", 1600, 800);
  //
  // c1->SetLogz();
  // // c2->SetLogz();
  // c3->SetLogz();
  // c4->SetLogz();
  // c5->SetLogz();
  // c6->SetLogz();
  // c7->SetLogz();
  // c8->SetLogz();
  // c9->SetLogz();
  // c10->SetLogz();
  //
  // c11->SetLogy();
  // c12->SetLogy();
  // c13->SetLogy();
  // c14->SetLogy();
  //
  // c15->SetLogy();
  // c16->SetLogy();
  // c17->SetLogy();
  // c18->SetLogy();
  //
  // c19->SetLogy();
  // c20->SetLogy();
  // c21->SetLogy();
  // c22->SetLogy();
  //
  // c23->SetLogy();
  // c24->SetLogy();
  // c25->SetLogy();
  // c26->SetLogy();
  //
  // c27->SetLogy();
  // c28->SetLogy();
  // c29->SetLogy();
  // c30->SetLogy();
  //
  // c31->SetLogy();
  // c32->SetLogy();
  // c33->SetLogy();
  // c34->SetLogy();
  //
  //
  // c1->cd();
  // H2_1->SetTitle("M^{2} VS p/q");
  // H2_1->GetXaxis()->SetTitle("p/q (GeV/c)");
  // H2_1->GetYaxis()->SetTitle("M^{2} (GeV/c^{2})^{2}");
  // H2_1->Draw("COLz");
  //
  // c2->cd();
  // H2_2->SetTitle("M^{2} VS pT/q");
  // H2_2->GetXaxis()->SetTitle("pT/q (GeV/c)");
  // H2_2->GetYaxis()->SetTitle("M^{2} (GeV/c^{2})^{2}");
  // H2_2->Draw("COLz");
  //
  // c3->cd();
  // H2_3->SetTitle("e^{+} pT VS y");
  // H2_3->GetXaxis()->SetTitle("y");
  // H2_3->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_3->Draw("COLz");
  //
  // c4->cd();
  // H2_4->SetTitle("K^{+} pT VS y");
  // H2_4->GetXaxis()->SetTitle("y");
  // H2_4->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_4->Draw("COLz");
  //
  // c5->cd();
  // H2_5->SetTitle("Pion^{+} pT VS y");
  // H2_5->GetXaxis()->SetTitle("y");
  // H2_5->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_5->Draw("COLz");
  //
  // c6->cd();
  // H2_6->SetTitle("Proton pT VS y");
  // H2_6->GetXaxis()->SetTitle("y");
  // H2_6->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_6->Draw("COLz");
  //
  // c7->cd();
  // H2_7->SetTitle("e^{-} pT VS y");
  // H2_7->GetXaxis()->SetTitle("y");
  // H2_7->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_7->Draw("COLz");
  //
  // c8->cd();
  // H2_8->SetTitle("K^{-} pT VS y");
  // H2_8->GetXaxis()->SetTitle("y");
  // H2_8->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_8->Draw("COLz");
  //
  // c9->cd();
  // H2_9->SetTitle("Pion^{-} pT VS y");
  // H2_9->GetXaxis()->SetTitle("y");
  // H2_9->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_9->Draw("COLz");
  //
  // c10->cd();
  // H2_10->SetTitle("Anti-proton pT VS y");
  // H2_10->GetXaxis()->SetTitle("y");
  // H2_10->GetYaxis()->SetTitle("pT/q (GeV/c)");
  // H2_10->Draw("COLz");
  //
  // c11->cd();
  // H1_1->SetTitle("e^{+} (mT - m0) Distribution");
  // H1_1->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_1->GetYaxis()->SetTitle("Counts");
  // H1_1->Draw();
  //
  // c12->cd();
  // H1_2->SetTitle("K^{+} (mT - m0) Distribution");
  // H1_2->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_2->GetYaxis()->SetTitle("Counts");
  // H1_2->Draw();
  //
  // c13->cd();
  // H1_3->SetTitle("Pion^{+} (mT - m0) Distribution");
  // H1_3->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_3->GetYaxis()->SetTitle("Counts");
  // H1_3->Draw();
  //
  // c14->cd();
  // H1_4->SetTitle("Proton^{+} (mT - m0) Distribution");
  // H1_4->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_4->GetYaxis()->SetTitle("Counts");
  // H1_4->Draw();
  //
  // c15->cd();
  // H1_5->SetTitle("e^{-} (mT - m0) Distribution");
  // H1_5->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_5->GetYaxis()->SetTitle("Counts");
  // H1_5->Draw();
  //
  // c16->cd();
  // H1_6->SetTitle("K^{-} (mT - m0) Distribution");
  // H1_6->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_6->GetYaxis()->SetTitle("Counts");
  // H1_6->Draw();
  //
  // c17->cd();
  // H1_7->SetTitle("Pion^{-} (mT - m0) Distribution");
  // H1_7->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_7->GetYaxis()->SetTitle("Counts");
  // H1_7->Draw();
  //
  // c18->cd();
  // H1_8->SetTitle("Proton^{-} (mT - m0) Distribution");
  // H1_8->GetXaxis()->SetTitle("mT-m0 (GeV/c)");
  // H1_8->GetYaxis()->SetTitle("Counts");
  // H1_8->Draw();
  //
  // c19->cd();
  // H1_9->SetTitle("e^{-} pT Distribution");
  // H1_9->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_9->GetYaxis()->SetTitle("Counts");
  // H1_9->Draw();
  //
  // c20->cd();
  // H1_10->SetTitle("K^{-} pT Distribution");
  // H1_10->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_10->GetYaxis()->SetTitle("Counts");
  // H1_10->Draw();
  //
  // c21->cd();
  // H1_11->SetTitle("Pion^{-} pT Distribution");
  // H1_11->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_11->GetYaxis()->SetTitle("Counts");
  // H1_11->Draw();
  //
  // c22->cd();
  // H1_12->SetTitle("Proton^{-} pT Distribution");
  // H1_12->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_12->GetYaxis()->SetTitle("Counts");
  // H1_12->Draw();
  //
  // c23->cd();
  // H1_13->SetTitle("e^{+} pT Distribution");
  // H1_13->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_13->GetYaxis()->SetTitle("Counts");
  // H1_13->Draw();
  //
  // c24->cd();
  // H1_14->SetTitle("K^{+} pT Distribution");
  // H1_14->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_14->GetYaxis()->SetTitle("Counts");
  // H1_14->Draw();
  //
  // c25->cd();
  // H1_15->SetTitle("Pion^{+} pT Distribution");
  // H1_15->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_15->GetYaxis()->SetTitle("Counts");
  // H1_15->Draw();
  //
  // c26->cd();
  // H1_16->SetTitle("Proton^{+} pT Distribution");
  // H1_16->GetXaxis()->SetTitle("pT (GeV/c)");
  // H1_16->GetYaxis()->SetTitle("Counts");
  // H1_16->Draw();
  //
  // c27->cd();
  // H1_17->SetTitle("e^{+} y Distribution");
  // H1_17->GetXaxis()->SetTitle("y");
  // H1_17->GetYaxis()->SetTitle("Counts");
  // H1_17->Draw();
  //
  // c28->cd();
  // H1_18->SetTitle("K^{+} y Distribution");
  // H1_18->GetXaxis()->SetTitle("y");
  // H1_18->GetYaxis()->SetTitle("Counts");
  // H1_18->Draw();
  //
  // c29->cd();
  // H1_19->SetTitle("Pion^{+} y Distribution");
  // H1_19->GetXaxis()->SetTitle("y");
  // H1_19->GetYaxis()->SetTitle("Counts");
  // H1_19->Draw();
  //
  // c30->cd();
  // H1_20->SetTitle("Proton^{+} y Distribution");
  // H1_20->GetXaxis()->SetTitle("y");
  // H1_20->GetYaxis()->SetTitle("Counts");
  // H1_20->Draw();
  //
  // c31->cd();
  // H1_21->SetTitle("e^{-} y Distribution");
  // H1_21->GetXaxis()->SetTitle("y");
  // H1_21->GetYaxis()->SetTitle("Counts");
  // H1_21->Draw();
  //
  // c32->cd();
  // H1_22->SetTitle("K^{-} y Distribution");
  // H1_22->GetXaxis()->SetTitle("y");
  // H1_22->GetYaxis()->SetTitle("Counts");
  // H1_22->Draw();
  //
  // c33->cd();
  // H1_23->SetTitle("Pion^{-} y Distribution");
  // H1_23->GetXaxis()->SetTitle("y");
  // H1_23->GetYaxis()->SetTitle("Counts");
  // H1_23->Draw();
  //
  // c34->cd();
  // H1_24->SetTitle("Proton^{-} y Distribution");
  // H1_24->GetXaxis()->SetTitle("y");
  // H1_24->GetYaxis()->SetTitle("Counts");
  // H1_24->Draw();

  // TFile * output =  new TFile("output.root","recreate");
  //
  // H2_1->Write();
  // H2_2->Write();
  //
  // H2_3->Write();
  // H2_4->Write();
  // H2_5->Write();
  // H2_6->Write();
  // H2_7->Write();
  // H2_8->Write();
  // H2_9->Write();
  // H2_10->Write();
  //
  // H1_1->Write();
  // H1_3->Write();
  // H1_2->Write();
  // H1_4->Write();
  // H1_5->Write();
  // H1_6->Write();
  // H1_7->Write();
  // H1_8->Write();
  // H1_9->Write();
  // H1_10->Write();
  // H1_11->Write();
  // H1_12->Write();
  // H1_13->Write();
  // H1_14->Write();
  // H1_15->Write();
  // H1_16->Write();
  // H1_17->Write();
  // H1_18->Write();
  // H1_19->Write();
  // H1_20->Write();
  // H1_21->Write();
  // H1_22->Write();
  // H1_24->Write();
  // H1_23->Write();

  return;
}
