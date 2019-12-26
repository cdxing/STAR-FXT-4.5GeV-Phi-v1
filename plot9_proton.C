//2019/7/1 Got from Yang Wu to analyze the fit on m2

// set how many Gaussian distributions to be seen in one distribution?
Int_t fitOrder = 3;

// pre-defined Gaussian function
Double_t GaussianFunction(Double_t *x, Double_t *par) {
    return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2.0));
}

// pre-defined multiple Gaussian function
Double_t twoGaussianFunction(Double_t *x, Double_t *par) {
    Double_t myFunc;
    for(Int_t i=0;i<2;i++) {
        if(i == 0) myFunc += GaussianFunction(x,par);
        myFunc += GaussianFunction(x,&par[3*i]);
    }
    return myFunc;
}

// pre-defined multiple Gaussian function
Double_t multipleGaussianFunction(Double_t *x, Double_t *par) {
    Double_t myFunc;
    for(Int_t i=0;i<fitOrder;i++) {
        if(i == 0) myFunc += GaussianFunction(x,par);
        myFunc += GaussianFunction(x,&par[3*i]);
    }
    return myFunc;
}

Double_t protonMass2 = 0.938272 * 0.938272;
Double_t pionMass2 = 0.139568 * 0.139568;
Double_t kaonMass2 = 0.493646 * 0.493646;
Double_t electronMass2 = 0.000510998 * 0.000510998;

void plot9_proton(Int_t momentumBin = 33) {

    TFile *inputFile1 = new TFile("all_mass2_total.root","read");

    TH1D *input_mass2_pT = (TH1D*)inputFile1->Get("extra_mass2_pT_2sigmaProton");
    Char_t name[100], description[200];
    TH1D *projection_mass2_in_pT_bin[60];
    for(Int_t pTbin=1;pTbin<=60;pTbin++) {
        sprintf(name,"projection_mass2_in_pT_bin%d",pTbin);
        sprintf(description,"m^{2} in %1.1f < p_{T} < %1.1f [GeV/c]",-3.1 + 0.1*pTbin,-3.0 + 0.1*pTbin);
        projection_mass2_in_pT_bin[pTbin-1] = new TH1D(name,description,500,-0.6,5.0);
        projection_mass2_in_pT_bin[pTbin-1]->GetXaxis()->SetTitle("m^{2} [GeV/c^{2}]^{2}");
        projection_mass2_in_pT_bin[pTbin-1]->GetYaxis()->SetTitle("# of tracks");
        for(Int_t mbin=1;mbin<=500;mbin++) {
            projection_mass2_in_pT_bin[pTbin-1]->SetBinContent(mbin,input_mass2_pT->GetBinContent(pTbin,mbin));
            projection_mass2_in_pT_bin[pTbin-1]->SetBinError(mbin,input_mass2_pT->GetBinError(pTbin,mbin));
        }
    }

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    gStyle->SetLegendBorderSize(0);

    TCanvas *c1 = new TCanvas("c1","c1",1024,768);
    //c1->Divide(3,2);
    // proton flow
    c1->cd(1);
    gPad->SetLeftMargin(0.1);gPad->SetRightMargin(0.03);gPad->SetBottomMargin(0.09);
    gPad->SetLogy(1);

    Int_t maxBin = projection_mass2_in_pT_bin[momentumBin-1]->GetMaximumBin();
    Double_t maxBinContent = projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(maxBin);
    projection_mass2_in_pT_bin[momentumBin-1]->GetXaxis()->SetRangeUser(-0.6,1.2);
    projection_mass2_in_pT_bin[momentumBin-1]->GetYaxis()->SetRangeUser(1.0,1e7);
    projection_mass2_in_pT_bin[momentumBin-1]->Draw();

    //gPad->Update();

    //projection_mass2_in_pT_bin[momentumBin-1]->Draw("same");

    TH1D *subRegion[5];
    TF1 *subFunc[10];
    TFitResultPtr r[5];
    Int_t usableDataPoints[5] = {0,0,0,0,0};
    Double_t width1 = 0.1;

    Double_t subRegion1_low = -0.05;
    Double_t subRegion1_high = 0.05;
    Double_t subRegion2_low =  0.17;
    Double_t subRegion2_high = 0.28;
    Double_t subRegion3_low =  0.8;
    Double_t subRegion3_high = 1.0;

    subRegion[1] = new TH1D("subRegion1","subRegion1",500,-0.6,5.0);
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion1_low);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion1_high);
    for(Int_t subBin=lowBin;subBin<=highBin;subBin++) {
        subRegion[1]->SetBinContent(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(subBin));
        subRegion[1]->SetBinError(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(subBin));
        usableDataPoints[1] += 1;
    }
  }
    subFunc[1] = new TF1("subFitFunc1",GaussianFunction,subRegion1_low,subRegion1_high,3);
    /*subFunc[1]->SetParameter(0,1e4);
    subFunc[1]->SetParLimits(0,1.0,subRegion[1]->GetMaximum());
    subFunc[1]->SetParameter(1,electronMass2);
    subFunc[1]->SetParLimits(1,electronMass2-0.002,electronMass2+0.002);
    subFunc[1]->SetParameter(2,subRegion[1]->GetRMS() * 0.5);
    subFunc[1]->SetParLimits(2,0.01,0.02);*/
    subFunc[1]->SetParameter(0,subRegion[1]->GetMaximum());
    subFunc[1]->SetParLimits(0,1.0,subRegion[1]->GetMaximum());
    subFunc[1]->SetParameter(1,pionMass2);
    subFunc[1]->SetParLimits(1,pionMass2-0.02,pionMass2+0.02);
    subFunc[1]->SetParameter(2,subRegion[1]->GetRMS() * 0.5);
    subFunc[1]->SetParLimits(2,0.01,0.1);
    r[1] = subRegion[1]->Fit("subFitFunc1","QS","",subRegion1_low,subRegion1_high);
    std::cout<<"Sub-section "<<1<<" fitting chi2 / dnf("<<usableDataPoints[1]-5<<") is "<<r[1]->Chi2()/(Double_t)(usableDataPoints[1] - 5)<<std::endl;

    subRegion[2] = new TH1D("subRegion2","subRegion1",500,-0.6,5.0);
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion2_low);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion2_high);
    for(Int_t subBin=lowBin;subBin<=highBin;subBin++) {
        subRegion[2]->SetBinContent(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(subBin));
        subRegion[2]->SetBinError(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(subBin));
        usableDataPoints[2] += 1;
    }
  }
    subFunc[2] = new TF1("subFitFunc2",GaussianFunction,subRegion2_low,subRegion2_high,3);
    subFunc[2]->SetParameter(0,subRegion[2]->GetMaximum());
    subFunc[2]->SetParLimits(0,1.0,subRegion[2]->GetMaximum());
    subFunc[2]->SetParameter(1,kaonMass2);
    subFunc[2]->SetParLimits(1,kaonMass2-0.025,kaonMass2+0.025);
    subFunc[2]->SetParameter(2,subRegion[2]->GetRMS());
    subFunc[2]->SetParLimits(2,0.01,0.1);
    r[2] = subRegion[2]->Fit("subFitFunc2","QS","",subRegion2_low,subRegion2_high);
    std::cout<<"Sub-section "<<2<<" fitting chi2 / dnf("<<usableDataPoints[2]-2<<") is "<<r[2]->Chi2()/(Double_t)(usableDataPoints[2] - 2)<<std::endl;

    subRegion[3] = new TH1D("subRegion3","subRegion1",500,-0.6,5.0);
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion3_low);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(subRegion3_high);
    for(Int_t subBin=lowBin;subBin<=highBin;subBin++) {
        subRegion[3]->SetBinContent(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(subBin));
        subRegion[3]->SetBinError(subBin,projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(subBin));
        usableDataPoints[3] += 1;
    }
  }
    subFunc[3] = new TF1("subFitFunc3",GaussianFunction,subRegion3_low,subRegion3_high,3);
    subFunc[3]->SetParameter(0,subRegion[3]->GetMaximum());
    subFunc[3]->SetParLimits(0,1.0,1e6);
    subFunc[3]->SetParameter(1,protonMass2);
    subFunc[3]->SetParLimits(1,protonMass2-0.05,protonMass2+0.05);
    subFunc[3]->SetParameter(2,subRegion[3]->GetRMS());
    subFunc[3]->SetParLimits(2,0.01,0.3);
    r[3] = subRegion[3]->Fit("subFitFunc3","QS","",subRegion3_low,subRegion3_high);
    std::cout<<"Sub-section "<<3<<" fitting chi2 / dnf("<<usableDataPoints[3]-2<<") is "<<r[3]->Chi2()/(Double_t)(usableDataPoints[3] - 2)<<std::endl;

    TF1 *f1 = new TF1("fitFunc",multipleGaussianFunction,-0.06,1.2,3*fitOrder);

    Double_t width2 = 0.1;
    f1->SetParameter(0,r[1]->Parameter(0)); f1->SetParLimits(0,1.0,subRegion[1]->GetMaximum());
    f1->SetParameter(1,pionMass2); f1->SetParLimits(1,pionMass2-width2,pionMass2+width2);
    f1->SetParameter(2,r[1]->Parameter(2)); f1->SetParLimits(2,0.01,0.5);
    // parameters of second gaussian
    f1->SetParameter(3,r[2]->Parameter(0)); f1->SetParLimits(3,1.0,1e6);
    f1->SetParameter(4,kaonMass2); f1->SetParLimits(4,kaonMass2-width2,kaonMass2+width2);
    f1->SetParameter(5,r[2]->Parameter(2)); f1->SetParLimits(5,0.01,0.5);
    // parameters of third gaussian
    f1->SetParameter(6,r[3]->Parameter(0)); f1->SetParLimits(6,1.0,1e6);
    f1->SetParameter(7,protonMass2); f1->SetParLimits(7,protonMass2-width2,protonMass2+width2);
    f1->SetParameter(8,r[3]->Parameter(2)); f1->SetParLimits(8,0.01,0.5);

    //TFitResultPtr r0 = projection_mass2_in_pT_bin[momentumBin-1]->Fit(f1,"S","",-10.0,10.0);

    projection_mass2_in_pT_bin[momentumBin-1]->GetXaxis()->SetRangeUser(-0.6,2.0);
    projection_mass2_in_pT_bin[momentumBin-1]->GetYaxis()->SetRangeUser(1.0,1e6);
    projection_mass2_in_pT_bin[momentumBin-1]->Draw();

    TBox *protonBox = new TBox(0.6,1,1.1,1e6);
    protonBox->SetFillStyle(3002);
    protonBox->SetLineWidth(1);
    protonBox->SetFillColor(8);
    protonBox->Draw("same");

    Double_t protonTotal = 0.0, protonTotalError = 0.0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.6);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(1.1);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        protonTotal += projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(bin);
        protonTotalError += projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(bin);
    }
  }

    TBox *kaonBox = new TBox(0.15,1,0.34,1e5);
    kaonBox->SetFillStyle(3002);
    kaonBox->SetLineWidth(1);
    kaonBox->SetFillColor(6);
    //kaonBox->Draw("same");

    Double_t kaonTotal = 0.0, kaonTotalError = 0.0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.15);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.34);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        kaonTotal += projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(bin);
        kaonTotalError += projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(bin);
    }
  }

    TBox *pionBox = new TBox(-0.1,1,0.12,1e5);
    pionBox->SetFillStyle(3002);
    pionBox->SetLineWidth(1);
    pionBox->SetFillColor(4);
    //pionBox->Draw("same");

    Double_t pionTotal = 0.0, pionTotalError = 0.0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(-0.1);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.12);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        pionTotal += projection_mass2_in_pT_bin[momentumBin-1]->GetBinContent(bin);
        pionTotalError += projection_mass2_in_pT_bin[momentumBin-1]->GetBinError(bin);
    }
  }

    projection_mass2_in_pT_bin[momentumBin-1]->Draw("same");

    /*subFunc[1]->SetLineColor(2);
    subFunc[1]->SetParameter(0,f1->GetParameter(0));
    subFunc[1]->SetParameter(1,f1->GetParameter(1));
    subFunc[1]->SetParameter(2,f1->GetParameter(2));
    subFunc[1]->Draw("same");

    subFunc[2]->SetLineColor(8);
    subFunc[2]->SetParameter(0,f1->GetParameter(3));
    subFunc[2]->SetParameter(1,f1->GetParameter(4));
    subFunc[2]->SetParameter(2,f1->GetParameter(5));
    subFunc[2]->Draw("same");

    subFunc[3]->SetLineColor(4);
    subFunc[3]->SetParameter(0,f1->GetParameter(6));
    subFunc[3]->SetParameter(1,f1->GetParameter(7));
    subFunc[3]->SetParameter(2,f1->GetParameter(8));
    subFunc[3]->Draw("same");*/

    Int_t precision = 1e4;

    /*subFunc[4] = new TF1("subFitFunc4",GaussianFunction,-0.6,5.0,3);
    subFunc[4]->SetLineColor(8);
    subFunc[4]->SetParameter(0,subFunc[1]->GetParameter(0));
    subFunc[4]->SetParameter(1,subFunc[1]->GetParameter(1));
    subFunc[4]->SetParameter(2,subFunc[1]->GetParameter(2));
    subFunc[4]->SetParError(0,subFunc[1]->GetParError(0));
    subFunc[4]->SetParError(1,subFunc[1]->GetParError(1));
    subFunc[4]->SetParError(2,subFunc[1]->GetParError(2));
    subFunc[4]->SetNpx(precision);
    subFunc[4]->Draw("same");*/

    subFunc[5] = new TF1("subFitFunc5",GaussianFunction,-0.6,5.0,3);
    subFunc[5]->SetLineColor(2);
    subFunc[5]->SetParameter(0,subFunc[1]->GetParameter(0));
    subFunc[5]->SetParameter(1,subFunc[1]->GetParameter(1));
    subFunc[5]->SetParameter(2,subFunc[1]->GetParameter(2));
    subFunc[5]->SetParError(0,subFunc[1]->GetParError(0));
    subFunc[5]->SetParError(1,subFunc[1]->GetParError(1));
    subFunc[5]->SetParError(2,subFunc[1]->GetParError(2));
    subFunc[5]->SetNpx(precision);
    subFunc[5]->Draw("same");

    subFunc[6] = new TF1("subFitFunc6",GaussianFunction,-0.6,5.0,3);
    subFunc[6]->SetLineColor(4);
    subFunc[6]->SetParameter(0,subFunc[2]->GetParameter(0));
    subFunc[6]->SetParameter(1,subFunc[2]->GetParameter(1));
    subFunc[6]->SetParameter(2,subFunc[2]->GetParameter(2));
    subFunc[6]->SetParError(0,subFunc[2]->GetParError(0));
    subFunc[6]->SetParError(1,subFunc[2]->GetParError(1));
    subFunc[6]->SetParError(2,subFunc[2]->GetParError(2));
    subFunc[6]->SetNpx(precision);
    subFunc[6]->Draw("same");

    subFunc[7] = new TF1("subFitFunc7",GaussianFunction,-0.6,5.0,3);
    subFunc[7]->SetLineColor(6);
    subFunc[7]->SetParameter(0,subFunc[3]->GetParameter(0));
    subFunc[7]->SetParameter(1,subFunc[3]->GetParameter(1));
    subFunc[7]->SetParameter(2,subFunc[3]->GetParameter(2));
    subFunc[7]->SetParError(0,subFunc[3]->GetParError(0));
    subFunc[7]->SetParError(1,subFunc[3]->GetParError(1));
    subFunc[7]->SetParError(2,subFunc[3]->GetParError(2));
    subFunc[7]->SetNpx(precision);
    subFunc[7]->Draw("same");

    Double_t pionFit = 0.0, pionFitError = 0.0, pionFitX[500], pionFitXerror[500]; Int_t pionFitPointCounter = 0;
    Double_t kaonContamination = 0.0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(-0.1);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.12);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        pionFit += subFunc[5]->Eval(projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin));
        kaonContamination += subFunc[6]->Eval(projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin));
        pionFitX[pionFitPointCounter] = projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin); pionFitPointCounter++;
    }
    r[1]->GetConfidenceIntervals(pionFitPointCounter, 1, 1, pionFitX, pionFitXerror, 0.683, false);
    for(Int_t bin=0;bin<pionFitPointCounter;bin++) {
        pionFitError += pionFitXerror[bin];
    }
  }
    Double_t pionPurity = (pionFit - kaonContamination) / pionTotal * 100.0;
    Double_t pionPurityError = pionPurity - (pionFit - kaonContamination - pionFitError) / (pionTotal + pionTotalError) * 100.0;
    cout<<"Pion PID purity is "<<pionPurity<<" +- "<<pionPurityError<<endl;

    Double_t kaonFit = 0.0, kaonFitError = 0.0, kaonFitX[500], kaonFitXerror[500]; Int_t kaonFitPointCounter = 0;
    Double_t pionContamination = 0.0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.15);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.34);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        kaonFit += subFunc[6]->Eval(projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin));
        pionContamination += subFunc[5]->Eval(projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin));
        kaonFitX[kaonFitPointCounter] = projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin); kaonFitPointCounter++;
    }
    r[2]->GetConfidenceIntervals(kaonFitPointCounter, 1, 1, kaonFitX, kaonFitXerror, 0.683, false);
    for(Int_t bin=0;bin<kaonFitPointCounter;bin++) {
        kaonFitError += kaonFitXerror[bin];
    }
  }
    Double_t kaonPurity = (kaonFit - pionContamination) / kaonTotal * 100.0;
    Double_t kaonPurityError = kaonPurity - (kaonFit - pionContamination - kaonFitError) / (kaonTotal + kaonTotalError) * 100.0;
    cout<<"Kaon PID purity is "<<kaonPurity<<" +- "<<kaonPurityError<<endl;

    Double_t protonFit = 0.0, protonFitError = 0.0, protonFitX[500], protonFitXerror[500]; Int_t protonFitPointCounter = 0;
  {
    Int_t lowBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(0.6);
    Int_t highBin = projection_mass2_in_pT_bin[momentumBin-1]->FindBin(1.1);
    for(Int_t bin=lowBin;bin<=highBin;bin++) {
        protonFit += subFunc[7]->Eval(projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin));
        protonFitX[protonFitPointCounter] = projection_mass2_in_pT_bin[momentumBin-1]->GetBinCenter(bin); protonFitPointCounter++;
    }
    r[3]->GetConfidenceIntervals(protonFitPointCounter, 1, 1, protonFitX, protonFitXerror, 0.683, false);
    for(Int_t bin=0;bin<protonFitPointCounter;bin++) {
        protonFitError += protonFitXerror[bin];
    }
  }
    Double_t protonPurity = protonFit / protonTotal * 100.0;
    Double_t protonPurityError = protonPurity - (protonFit - protonFitError) / (protonTotal + protonTotalError) * 100.0;
    cout<<"Proton PID purity is "<<protonPurity<<" +- "<<protonPurityError<<endl;

    f1->SetLineColor(6);
    //f1->Draw("same");

    TLegend *leg = new TLegend(0.674,0.537,0.944,0.890);
    leg->SetFillStyle(0);
    leg->SetMargin(0.15);
    leg->SetTextSize(0.03);
    leg->AddEntry(projection_mass2_in_pT_bin[momentumBin-1],"m^{2} data","lpe");
    //leg->AddEntry(subFunc[4],"Gaussian fitting for electron","l");
    leg->AddEntry(subFunc[5],"Gaussian fitting for pion","l");
    leg->AddEntry(subFunc[6],"Gaussian fitting for kaon","l");
    leg->AddEntry(subFunc[7],"Gaussian fitting for proton","l");
    leg->AddEntry(pionBox,"pion m^{2} cut","f");
    leg->AddEntry(kaonBox,"kaon m^{2} cut","f");
    leg->AddEntry(protonBox,"proton m^{2} cut","f");
    leg->Draw();

    TPaveText *t=new TPaveText(0.146,0.783,0.465,0.871,"nbNDC");
    t->SetFillStyle(0);t->SetFillColor(0);t->SetLineWidth(0);
    t->SetTextFont(62);t->SetTextSize(0.03);
    t->SetTextAlign(13);
    t->AddText("Proton PID purity: ");
    sprintf(name,"%2.2f%% #pm %1.2f%%",protonPurity,protonPurityError);
    t->AddText(name);
    /*t->AddText("");
    t->AddText("Pion PID purity: ");
    sprintf(name,"%2.2f%% #pm %1.2f%%",pionPurity,pionPurityError);
    t->AddText(name);
    t->AddText("");
    t->AddText("Kaon PID purity: ");
    sprintf(name,"%2.2f%% #pm %1.2f%%",kaonPurity,kaonPurityError);
    t->AddText(name);
    t->AddText("");*/
    t->Draw();

    sprintf(name,"proton_pid_purity_check_pTbin%d.root",momentumBin);
    TFile *outputFile1 = new TFile(name,"recreate");
    TH1D *proton_purity = (TH1D*)outputFile1->Get("proton_purity");
    TH1D *pion_purity = (TH1D*)outputFile1->Get("pion_purity");
    TH1D *kaon_purity = (TH1D*)outputFile1->Get("kaon_purity");
    if(!proton_purity) {
        proton_purity = new TH1D("proton_purity","Proton PID purity",60,-3.0,3.0);
        proton_purity->GetXaxis()->SetTitle("Transverse momentum p_{T} [GeV/c]");
        proton_purity->GetYaxis()->SetTitle("PID purity [%]");
    }
    proton_purity->SetBinContent(momentumBin,protonPurity);
    proton_purity->SetBinError(momentumBin,protonPurityError);
    if(!pion_purity) {
        pion_purity = new TH1D("pion_purity","Pion PID purity",60,-3.0,3.0);
        pion_purity->GetXaxis()->SetTitle("Transverse momentum p_{T} [GeV/c]");
        pion_purity->GetYaxis()->SetTitle("PID purity [%]");
    }
    pion_purity->SetBinContent(momentumBin,pionPurity);
    pion_purity->SetBinError(momentumBin,pionPurityError);
    if(!kaon_purity) {
        kaon_purity = new TH1D("kaon_purity","Kaon PID purity",60,-3.0,3.0);
        kaon_purity->GetXaxis()->SetTitle("Transverse momentum p_{T} [GeV/c]");
        kaon_purity->GetYaxis()->SetTitle("PID purity [%]");
    }
    kaon_purity->SetBinContent(momentumBin,kaonPurity);
    kaon_purity->SetBinError(momentumBin,kaonPurityError);
    //outputFile1->Write();

    sprintf(name,"fxt_data_pid_mass2_pTbin%d.pdf",momentumBin);
    //c1->Print(name);
}
