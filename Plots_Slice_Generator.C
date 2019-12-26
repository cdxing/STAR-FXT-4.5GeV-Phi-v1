#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"

void Plots_Slice_Generator(const char* inFile = "./result/0FCC1C2E3DC45954DCA54F56103F1937.picoDst.result.combined.root")
{
  TFile* F = TFile::Open("./result/0FCC1C2E3DC45954DCA54F56103F1937.picoDst.result.combined.root");
  if (!F || F->IsZombie()){
    std::cout<<"Cannot open the file"<<inFile<<std::endl;
    return;
  }

  TH2D* H2 = (TH2D*)F->Get("h2_m2_QA_pT");
  TH1D* H1_slice[100];
  for(Int_t ibin = 0; ibin < 100; ibin++){
      H1_slice[ibin] = H2->ProjectionY(Form("Slice_%d",ibin),ibin*50+1,ibin*50+50);//->Draw();
  }

  return;
}
