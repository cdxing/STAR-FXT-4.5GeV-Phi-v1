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
double d_norm = 0.00286541;//sideTest//0.00296682;//primary//0.00299232;//invM_eff_corr//0.0159051;//rebin//0.0159;//5mx_test//0.57346;
double d_norm_y_bin2 = 0.00256758;//sideTest//0.00273469;//primary//0.00283043;//invM_eff_corr//0.0161716;//rebin//0.0161451;//5mx_test//0.524116;
double d_norm_y_bin3 = 0.00298416;//sideTest//0.00306014;//primary//0.00304982;//invM_eff_corr//0.0156764;//rebin//0.0156756;//5mx_test//0.589171;
double d_norm_y_bin4 = 0.00303191;//sideTest//0.00308287;//primary//0.00293368;//invM_eff_corr//0.0181312;//rebin//0.0181722;//5mx_test//0.680254;

double d_p0 = 117.73;//sideTest//121.597;//primary//353.15;//invM_eff_corr//125.267;//rebin//62.6622;//5mx_test//62.4311;//secondfit//63.8139;//prefit//
double d_mean  = 1.0208;//sideTest//1.02087;//primary//1.02104;//invM_eff_corr//1.02094;//rebin//1.02094;//5mx_test//1.02017;//secondfit//1.0202;//prefit//
double d_sigma = 0.00435207;//sideTest//0.00444335;//primary//0.00492315;//invM_eff_corr//0.00452375;//rebin//0.00452732;//5mx_test//0.00463761;//secondfit//0.00483768;//prefit//
double d_p3 = 0;//2.03306//sideTest////2.22375//primary////7.57194;//invM_eff_corr////3.88171//rebin////1.96872//5mx_test////1.90684;//secondfit//

double d_p0_y_bin2 = 33.1623;//sideTest//36.8197;//primary//83.9407;//invM_eff_corr//37.6058;//rebin//18.8268;//5mx_test//20.918;//secondfit//20.8365;//prefit//
double d_mean_y_bin2 = 1.02249;//sideTest//1.02248;//primary//1.02267;//invM_eff_corr//1.02211;//rebin//1.02211;//5mx_test//1.0227;//secondfit//1.02278;//prefit//
double d_sigma_y_bin2 = 0.00420342;//sideTest//0.0045509;//primary//0.00471844;//invM_eff_corr//0.0047948;//rebin//0.004826;//5mx_test//0.0035232;//secondfit//0.00339162;//prefit//
double d_p3_y_bin2 = 0;//-1.09335//sideTest////-1.20777//primary////-2.46834;//invM_eff_corr////-1.55451//rebin////-0.75103//5mx_test////-0.332631;//secondfit//

double d_p0_y_bin3 = 81.2566;//sideTest//81.8149;//primary//252.679;//invM_eff_corr//82.6144;//rebin//41.3101;//5mx_test//45.3233;//secondfit//47.0577;//prefit//
double d_mean_y_bin3 = 1.02004;//sideTest//1.02005;//primary//1.0204;//invM_eff_corr//1.02019;//rebin//1.02019;//5mx_test//1.01942;//secondfit//1.01951;//prefit//
double d_sigma_y_bin3 = 0.0039893;//sideTest//0.00397084;//primary//0.00435745;//invM_eff_corr//0.00400454;//rebin//0.00400505;//5mx_test//0.00403223;//secondfit//0.00435069;//prefit//
double d_p3_y_bin3 = 0;//2.34521//sideTest////2.7601//primary////10.36;//invM_eff_corr////5.68407//rebin////2.84558//5mx_test////2.4756;//secondfit//

double d_p0_y_bin4 = 9.77643;//sideTest//9.64978;//primary//34.7556;//invM_eff_corr//10.2513;//rebin//5.10356;//5mx_test//4.89736;//secondfit//4.08391;//prefit//
double d_mean_y_bin4 = 1.02195;//sideTest//1.02201;//primary//1.02046;//invM_eff_corr//1.02269;//rebin//1.0227;//5mx_test//1.01759;//secondfit//1.01721;//prefit//
double d_sigma_y_bin4 = 0.00547513;//sideTest//0.00551411;//primary//0.00651522;//invM_eff_corr//0.00524815;//rebin//0.00522871;//5mx_test//0.0031585;//secondfit//0.00240508;//prefit//
double d_p3_y_bin4 = 0;//-0.46391//sideTest////-0.471533//primary////-1.80907;//invM_eff_corr////-0.623697//rebin////-0.312837//5mx_test////-0.983175;//secondfit//



double d_pol_scale_p0 = -60933.1-(2.22375);//primary//-142085-(7.57194);//invM_eff_corr//-323.52-(3.88171);//rebin//-161.708-(1.96872);//5mx_test//-381.415-(1.90684);//-54570.5;//
double d_pol_scale_p1 = 118766;//primary//276859;//invM_eff_corr//132.954;//rebin//66.4555;//5mx_test//66.4461;//104032;//
double d_pol_scale_p2 = -57476.4;//primary//-133800;//invM_eff_corr//568.54;//rebin//284.179;//5mx_test//495.591;//-49362.2;//

double d_pol_scale_p0_y_bin2 = -10.6493-(-1.20777);//primary//-34063.3-(-2.46834);//invM_eff_corr//-69.986-(-1.55451);//rebin//-34.9357-(-0.75103);//5mx_test//-166.539-(-0.332631);//-60.0655;//
double d_pol_scale_p1_y_bin2 = 35.0802;//primary//66632.8;//invM_eff_corr//34.9195;//rebin//17.4311;//5mx_test//17.6437;//17.4833;//
double d_pol_scale_p2_y_bin2 = 78.4821;//primary//-32355.3;//invM_eff_corr//134.999;//rebin//67.3891;//5mx_test//194.316;//91.3819;//

double d_pol_scale_p0_y_bin3 = -349.919-(2.7601);//primary//-663.671-(10.36);//invM_eff_corr//-297.011-(5.68407);//rebin//-148.498-(2.84558);//5mx_test//-93.4751-(2.4756);//-240.585;//
double d_pol_scale_p1_y_bin3 = 91.6685;//primary//258.427;//invM_eff_corr//91.192;//rebin//45.5937;//5mx_test//44.8648;//44.8336;//
double d_pol_scale_p2_y_bin3 = 1143.29;//invM_eff_corr//515.638;//primary//463.317;//rebin//231.647;//5mx_test//178.059;//319.881;//

double d_pol_scale_p0_y_bin4 = -4829.64-(-0.471533);//primary//-17279.9-(-1.80907);//invM_eff_corr//-32.4335-(-0.623697);//rebin//-16.2534-(-0.312837);//5mx_test//-118.027-(-0.983175);//-73.7018;//
double d_pol_scale_p1_y_bin4 = 9367.01;//primary//33270.5;//invM_eff_corr//6.66661;//rebin//3.34084;//5mx_test//3.49542;//3.37993;//
double d_pol_scale_p2_y_bin4 = -4521.97;//primary//-15939.1;//invM_eff_corr//43.9982;//rebin//22.0489;//5mx_test//121.323;//78.3075;//

// double d_bg_v_p0 = 1.61018;
// double d_bg_v_p1 = -3.05911;
// double d_bg_v_p2 = 1.45088;



double bg_func(double *x, double *p)
{
  if( x[0] >= (d_mean - (3*d_sigma)) && x[0] <= (d_mean + (3*d_sigma)))
   {
     return (((d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))
     /((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)+(d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double bg_func_y_bin2(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin2 - (3*d_sigma_y_bin2)) && x[0] <= (d_mean_y_bin2 + (3*d_sigma_y_bin2)))
   {
     return (((d_pol_scale_p0_y_bin2 + d_pol_scale_p1_y_bin2 * x[0] + d_pol_scale_p2_y_bin2 * (x[0]**2))
     /((d_p0_y_bin2*exp(-0.5*((x[0]-d_mean_y_bin2)/d_sigma_y_bin2)**2)+d_p3_y_bin2)+(d_pol_scale_p0_y_bin2 + d_pol_scale_p1_y_bin2 * x[0] + d_pol_scale_p2_y_bin2 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double bg_func_y_bin3(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin3 - (3*d_sigma_y_bin3)) && x[0] <= (d_mean_y_bin3 + (3*d_sigma_y_bin3)))
   {
     return (((d_pol_scale_p0_y_bin3 + d_pol_scale_p1_y_bin3 * x[0] + d_pol_scale_p2_y_bin3 * (x[0]**2))
     /((d_p0_y_bin3*exp(-0.5*((x[0]-d_mean_y_bin3)/d_sigma_y_bin3)**2)+d_p3_y_bin3)+(d_pol_scale_p0_y_bin3 + d_pol_scale_p1_y_bin3 * x[0] + d_pol_scale_p2_y_bin3 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double bg_func_y_bin4(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin4 - (3*d_sigma_y_bin4)) && x[0] <= (d_mean_y_bin4 + (3*d_sigma_y_bin4)))
   {
     return (((d_pol_scale_p0_y_bin4 + d_pol_scale_p1_y_bin4 * x[0] + d_pol_scale_p2_y_bin4 * (x[0]**2))
     /((d_p0_y_bin4*exp(-0.5*((x[0]-d_mean_y_bin4)/d_sigma_y_bin4)**2)+d_p3_y_bin4)+(d_pol_scale_p0_y_bin4 + d_pol_scale_p1_y_bin4 * x[0] + d_pol_scale_p2_y_bin4 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/);
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double sig_bg_func(double *x, double *p)
{
  if( x[0] >= (d_mean - (3*d_sigma)) && x[0] <= (d_mean + (3*d_sigma)))
   {
     return (
       (p[/*1*//*2*/3/*4*/]*((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)
     /((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)+(d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2)))))
     +
       (((d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))
     /((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)+(d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double sig_bg_func_y_bin2(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin2 - (3*d_sigma_y_bin2)) && x[0] <= (d_mean_y_bin2 + (3*d_sigma_y_bin2)))
   {
     return (
       (p[/*1*//*2*/3/*4*/]*((d_p0_y_bin2*exp(-0.5*((x[0]-d_mean_y_bin2)/d_sigma_y_bin2)**2)+d_p3_y_bin2)
     /((d_p0_y_bin2*exp(-0.5*((x[0]-d_mean_y_bin2)/d_sigma_y_bin2)**2)+d_p3_y_bin2)+(d_pol_scale_p0_y_bin2 + d_pol_scale_p1_y_bin2 * x[0] + d_pol_scale_p2_y_bin2 * (x[0]**2)))))
     +
       (((d_pol_scale_p0_y_bin2 + d_pol_scale_p1_y_bin2 * x[0] + d_pol_scale_p2_y_bin2 * (x[0]**2))
     /((d_p0_y_bin2*exp(-0.5*((x[0]-d_mean_y_bin2)/d_sigma_y_bin2)**2)+d_p3_y_bin2)+(d_pol_scale_p0_y_bin2 + d_pol_scale_p1_y_bin2 * x[0] + d_pol_scale_p2_y_bin2 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double sig_bg_func_y_bin3(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin3 - (3*d_sigma_y_bin3)) && x[0] <= (d_mean_y_bin3 + (3*d_sigma_y_bin3)))
   {
     return (
      ( p[/*1*//*2*/3/*4*/]*((d_p0_y_bin3*exp(-0.5*((x[0]-d_mean_y_bin3)/d_sigma_y_bin3)**2)+d_p3_y_bin3)
     /((d_p0_y_bin3*exp(-0.5*((x[0]-d_mean_y_bin3)/d_sigma_y_bin3)**2)+d_p3_y_bin3)+(d_pol_scale_p0_y_bin3 + d_pol_scale_p1_y_bin3 * x[0] + d_pol_scale_p2_y_bin3 * (x[0]**2)))))
     +
       (((d_pol_scale_p0_y_bin3 + d_pol_scale_p1_y_bin3 * x[0] + d_pol_scale_p2_y_bin3 * (x[0]**2))
     /((d_p0_y_bin3*exp(-0.5*((x[0]-d_mean_y_bin3)/d_sigma_y_bin3)**2)+d_p3_y_bin3)+(d_pol_scale_p0_y_bin3 + d_pol_scale_p1_y_bin3 * x[0] + d_pol_scale_p2_y_bin3 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

double sig_bg_func_y_bin4(double *x, double *p)
{
  if( x[0] >= (d_mean_y_bin4 - (3*d_sigma_y_bin4)) && x[0] <= (d_mean_y_bin4 + (3*d_sigma_y_bin4)))
   {
     return (
       (p[/*1*//*2*/3/*4*/]*((d_p0_y_bin4*exp(-0.5*((x[0]-d_mean_y_bin4)/d_sigma_y_bin4)**2)+d_p3_y_bin4)
     /((d_p0_y_bin4*exp(-0.5*((x[0]-d_mean_y_bin4)/d_sigma_y_bin4)**2)+d_p3_y_bin4)+(d_pol_scale_p0_y_bin4 + d_pol_scale_p1_y_bin4 * x[0] + d_pol_scale_p2_y_bin4 * (x[0]**2)))))
     +
       (((d_pol_scale_p0_y_bin4 + d_pol_scale_p1_y_bin4 * x[0] + d_pol_scale_p2_y_bin4 * (x[0]**2))
     /((d_p0_y_bin4*exp(-0.5*((x[0]-d_mean_y_bin4)/d_sigma_y_bin4)**2)+d_p3_y_bin4)+(d_pol_scale_p0_y_bin4 + d_pol_scale_p1_y_bin4 * x[0] + d_pol_scale_p2_y_bin4 * (x[0]**2))))
     */*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/));
   }
   else
   {
     return /*p[0]*//*(p[0] + p[1] * x[0])*/(p[0] + p[1] * x[0] + p[2] * (x[0]**2))/*(p[0] + p[1] * x[0] + p[2] * (x[0]**2) + p[3]*(x[0]**3))*/;
   }
}

// double sig_bg_func(double *x, double *p)
// {
//   if( x[0] >= (d_mean - (3*d_sigma)) && x[0] <= (d_mean + (3*d_sigma)))
//    {
//      return (
//        p[0]*((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)
//      /((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)+(d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))))
//        ((d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))
//      /((d_p0*exp(-0.5*((x[0]-d_mean)/d_sigma)**2)+d_p3)+(d_pol_scale_p0 + d_pol_scale_p1 * x[0] + d_pol_scale_p2 * (x[0]**2))))
//      *(d_bg_v_p0 + d_bg_v_p1 * x[0] + d_bg_v_p2 * (x[0]**2)));
//    }
//    else
//    {
//      return (d_bg_v_p0 + d_bg_v_p1 * x[0] + d_bg_v_p2 * (x[0]**2));
//    }
// }

void phiFlowFitting(
  string FileName,
  TString outFile = "test"
)
{
  // int beginIDx = FileName.find_last_of("\\/");
  // string outFileName = FileName.substr(beginIDx+1);
  // std::cout<< "outFileName: " <<outFileName << std::endl;

  outFile.Append(Form("_%s",FileName));
  outFile.Prepend("Phi_flow_fitting_out_");
  TFile * tf_out = new TFile(outFile,"RECREATE");
  // TFile * tf_mx = new TFile("h_mx_TPC2_TOF3_invM_out_tot.root","READ");
  // TFile * tf_mx  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_scheduler/merged_5B7367B5663C88E112AB01E56C95D38A.root","READ");
  TFile * tf_mx  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_systematic_error_analysis_invM/phi_meson_mixed_events_invariant_mass_primary.root","READ");
// TFile * tf_phi_out_single  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/Phi_out_single_tot.root","READ");//nm events
  TFile * tf_phi_out_single  = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/phiAna/phiInvMass/result_systematic_error_analysis_invM/phi_meson_normal_event_invariant_mass_primary.root","READ");//nm events

  TH1D  * h_mx_TPC2_TOF3_invM = (TH1D*) tf_mx -> Get("h_mx_TPC2_TOF3_invM");
  h_mx_TPC2_TOF3_invM->Rebin();
  // TFile * tf_mx_y_bin2 = new TFile("h_mx_TPC2_TOF3_invM_y_bin2_out_tot.root","READ");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin2 = (TH1D*) tf_mx -> Get("h_mx_TPC2_TOF3_invM_y_bin2");
  h_mx_TPC2_TOF3_invM_y_bin2->Rebin();

  // TFile * tf_mx_y_bin3 = new TFile("h_mx_TPC2_TOF3_invM_y_bin3_out_tot.root","READ");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin3 = (TH1D*) tf_mx -> Get("h_mx_TPC2_TOF3_invM_y_bin3");
  h_mx_TPC2_TOF3_invM_y_bin3->Rebin();

  // TFile * tf_mx_y_bin4 = new TFile("h_mx_TPC2_TOF3_invM_y_bin4_out_tot.root","READ");
  TH1D  * h_mx_TPC2_TOF3_invM_y_bin4 = (TH1D*) tf_mx -> Get("h_mx_TPC2_TOF3_invM_y_bin4");
  h_mx_TPC2_TOF3_invM_y_bin4->Rebin();


  // TFile * tf_nm = new TFile("h_TPC2_TOF3_invM_out_tot.root","READ");
  TH1D  * h_TPC2_TOF3_invM = (TH1D*) tf_phi_out_single -> Get("h_TPC2_TOF3_invM");
  // TFile * tf_nm_y_bin2 = new TFile("h_TPC2_TOF3_invM_y_bin2_out_tot.root","READ");
  TH1D  * h_TPC2_TOF3_invM_y_bin2 = (TH1D*) tf_phi_out_single -> Get("h_TPC2_TOF3_invM_y_bin2");
  // TFile * tf_nm_y_bin3 = new TFile("h_TPC2_TOF3_invM_y_bin3_out_tot.root","READ");
  TH1D  * h_TPC2_TOF3_invM_y_bin3 = (TH1D*) tf_phi_out_single -> Get("h_TPC2_TOF3_invM_y_bin3");
  // TFile * tf_nm_y_bin4 = new TFile("h_TPC2_TOF3_invM_y_bin4_out_tot.root","READ");
  TH1D  * h_TPC2_TOF3_invM_y_bin4 = (TH1D*) tf_phi_out_single -> Get("h_TPC2_TOF3_invM_y_bin4");

  // TFile * tf_phi_v1_vs_invM = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/files/bbcEventPlane/phi_v1_vs_invM_first_run.root","READ");
  // TFile * tf_phi_v1_vs_invM = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_systematic_error_analysis_directed_flow/phi_meson_directed_flow_primary.root","READ");
  // if( FileName.find(".list") != string::npos ||
  //     FileName.find(".lis") != string::npos ) {
  //
  //   std::ifstream inputStream( FileName.c_str() );
  //
  //   if(!inputStream) {
  //     cout << "ERROR: Cannot open list file " << FileName << endl;
  //   }
  //
  //   Int_t nFile = 0;
  //   string file;
  //   while(getline(inputStream, file)) {
  //     if(file.find(".root") != string::npos) {
  //       TFile* ftmp = TFile::Open(file.c_str());
  //       if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) {
  //         cout << " Read in root file " << file << endl;
  //         // t_K->Add(file.c_str());
  //         ++nFile;
  //       } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())
  //
  //       if (ftmp) {
  //         ftmp->Close();
  //       } //if (ftmp)
  //     } //if(file.find(".picoDst.root") != string::npos)
  //   } //while (getline(inputStream, file))
  //
  //   cout << " Total " << nFile << " files have been read in. " << endl;
  // } //if(FileName.find(".list") != string::npos || FileName.find(".lis" != string::npos))
  // else if(FileName.find(".root") != string::npos) {
  //   // t_K->Add(FileName.c_str());
  // }
  // else {
  //   cout << " No good input file to read ... " << endl;
  // }


  // TFile * tf_phi_v1_vs_invM = new TFile(Form("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_version4/%s",FileName),"READ");
  // TFile * tf_phi_v1_vs_invM = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_systematic_error_analysis_directed_flow/phi_meson_directed_flow_primary.root","READ");
  TFile * tf_phi_v1_vs_invM = new TFile("/star/data01/pwg/dchen/Ana/fxtPicoAna/result_test/test_phi_v1_vs_invM.root","READ");

  if(tf_phi_v1_vs_invM && !tf_phi_v1_vs_invM->IsZombie()) {
    std::cout << " Read in root file " << tf_phi_v1_vs_invM << std::endl;
  }

  //  hadd_18_10_0_.root
  TProfile  * tp_phi_v1_vs_invM_pfx = (TProfile*) tf_phi_v1_vs_invM -> Get("h2_phi_v1_vs_invM_pfx");

  TProfile  * tp_phi_v1_vs_invM_pfx_y_bin2 = (TProfile*) tf_phi_v1_vs_invM -> Get("h2_phi_v1_vs_invM_bin2_pfx");
  TProfile  * tp_phi_v1_vs_invM_pfx_y_bin3 = (TProfile*) tf_phi_v1_vs_invM -> Get("h2_phi_v1_vs_invM_bin3_pfx");
  TProfile  * tp_phi_v1_vs_invM_pfx_y_bin4 = (TProfile*) tf_phi_v1_vs_invM -> Get("h2_phi_v1_vs_invM_bin4_pfx");

  TH1D  *h_phi_v1_vs_y = new TH1D("h_phi_v1_vs_y","v1 of #phi meson VS rapidity",6,-3.0,0);
  h_phi_v1_vs_y->GetXaxis()->SetTitle("y");
  h_phi_v1_vs_y->GetYaxis()->SetTitle("v1");

  // Fitting the backgraound
  TF1 * tf1_pol = new TF1("tf1_pol","[0] + [1]*x + [2]*x**2",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));
  TH1D  * h_mx_TPC2_TOF3_invM_same = h_mx_TPC2_TOF3_invM;

  TH1D  * h_mx_TPC2_TOF3_invM_rebin = h_mx_TPC2_TOF3_invM;//->Rebin();
  gStyle->SetOptDate(0);
  TCanvas * TC_inv = new TCanvas("TC_inv","TC_inv",1280,720);
  TC_inv->cd();
  h_mx_TPC2_TOF3_invM_rebin->Draw();
  h_mx_TPC2_TOF3_invM_rebin->Fit(tf1_pol,"0","R",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));

  TF1 * tf1_pol_y_bin2 = new TF1("tf1_pol_y_bin2","[0] + [1]*x + [2]*x**2",(d_mean_y_bin2 - (3*d_sigma_y_bin2)),(d_mean_y_bin2 + (3*d_sigma_y_bin2)));
  TH1D  * h_mx_TPC2_TOF3_invM_same_y_bin2 = h_mx_TPC2_TOF3_invM_y_bin2;

  TH1D  * h_mx_TPC2_TOF3_invM_rebin_y_bin2 = h_mx_TPC2_TOF3_invM_y_bin2;//->Rebin();
  h_mx_TPC2_TOF3_invM_y_bin2->SetFillColor(kYellow);
  TCanvas * TC_inv_y_bin2 = new TCanvas("TC_inv_y_bin2","TC_inv_y_bin2",1280,720);
  TC_inv_y_bin2->cd();
  h_mx_TPC2_TOF3_invM_rebin_y_bin2->Draw();
  h_mx_TPC2_TOF3_invM_rebin_y_bin2->Fit(tf1_pol_y_bin2,"E","R",(d_mean_y_bin2 - (3*d_sigma_y_bin2)),(d_mean_y_bin2 + (3*d_sigma_y_bin2)));
  // h_TPC2_TOF3_invM_y_bin2->Draw();
  // h_mx_TPC2_TOF3_invM_y_bin2 -> Scale(d_norm_y_bin2);
  //
  // h_mx_TPC2_TOF3_invM_y_bin2->Draw("same");

  TF1 * tf1_pol_y_bin3 = new TF1("tf1_pol_y_bin3","[0] + [1]*x + [2]*x**2",(d_mean_y_bin3 - (3*d_sigma_y_bin3)),(d_mean_y_bin3 + (3*d_sigma_y_bin3)));
  TH1D  * h_mx_TPC2_TOF3_invM_rebin_y_bin3 = h_mx_TPC2_TOF3_invM_y_bin3;//->Rebin();
  TCanvas * TC_inv_y_bin3 = new TCanvas("TC_inv_y_bin3","TC_inv_y_bin3",1280,720);
  TC_inv_y_bin3->cd();
  h_mx_TPC2_TOF3_invM_rebin_y_bin3->Draw();
  h_mx_TPC2_TOF3_invM_rebin_y_bin3->Fit(tf1_pol_y_bin3,"E","R",(d_mean_y_bin3 - (3*d_sigma_y_bin3)),(d_mean_y_bin3 + (3*d_sigma_y_bin3)));

  TF1 * tf1_pol_y_bin4 = new TF1("tf1_pol_y_bin4","[0] + [1]*x + [2]*x**2",(d_mean_y_bin4 - (3*d_sigma_y_bin4)),(d_mean_y_bin4 + (3*d_sigma_y_bin4)));
  TH1D  * h_mx_TPC2_TOF3_invM_rebin_y_bin4 = h_mx_TPC2_TOF3_invM_y_bin4;//->Rebin();
  TCanvas * TC_inv_y_bin4 = new TCanvas("TC_inv_y_bin4","TC_inv_y_bin4",1280,720);
  TC_inv_y_bin4->cd();
  h_mx_TPC2_TOF3_invM_rebin_y_bin4->Draw();
  h_mx_TPC2_TOF3_invM_rebin_y_bin4->Fit(tf1_pol_y_bin4,"E","R",(d_mean_y_bin4 - (3*d_sigma_y_bin4)),(d_mean_y_bin4 + (3*d_sigma_y_bin4)));
  // END Fitting the backgraound

  // Scale the fitting function
  TF1 * tf1_pol_scale = new TF1("tf1_pol_scale","[0] + [1]*x + [2]*x**2",(d_mean - (3*d_sigma)),(d_mean + (3*d_sigma)));

  double  d_pol_scale_p0 = (tf1_pol->GetParameter(0)) /** 0.5*/ * d_norm;
  double  d_pol_scale_p1 = (tf1_pol->GetParameter(1)) /** 0.5*/ * d_norm;
  double  d_pol_scale_p2 = (tf1_pol->GetParameter(2)) /** 0.5*/ * d_norm;

  cout << "d_pol_scale_p0 = "<< d_pol_scale_p0 << endl;
  cout << "d_pol_scale_p1 = "<< d_pol_scale_p1 << endl;
  cout << "d_pol_scale_p2 = "<< d_pol_scale_p2 << endl;

  tf1_pol_scale -> SetParameter(0,d_pol_scale_p0);
  tf1_pol_scale -> SetParameter(1,d_pol_scale_p1);
  tf1_pol_scale -> SetParameter(2,d_pol_scale_p2);
  TCanvas * TC_inv_bg = new TCanvas("TC_inv_bg","TC_inv_bg",1280,720);
  TC_inv_bg->cd();
  // TC_inv->cd();
  h_TPC2_TOF3_invM->Rebin();
  h_TPC2_TOF3_invM->Draw();
  h_mx_TPC2_TOF3_invM_same->Scale(/*0.5**/d_norm);
  h_mx_TPC2_TOF3_invM_same->SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_same-> SetFillStyle(3001);
  h_mx_TPC2_TOF3_invM_same->Draw("same");
  tf1_pol_scale -> SetLineColor(kBlue);
  tf1_pol_scale->Draw("same");

  // y_bin2
  TF1 * tf1_pol_scale_y_bin2 = new TF1("tf1_pol_scale_y_bin2","[0] + [1]*x + [2]*x**2",(d_mean_y_bin2 - (3*d_sigma_y_bin2)),(d_mean_y_bin2 + (3*d_sigma_y_bin2)));

  double  d_pol_scale_p0_y_bin2 = (tf1_pol_y_bin2->GetParameter(0)) /** 0.5*/ * d_norm_y_bin2;
  double  d_pol_scale_p1_y_bin2 = (tf1_pol_y_bin2->GetParameter(1)) /** 0.5*/ * d_norm_y_bin2;
  double  d_pol_scale_p2_y_bin2 = (tf1_pol_y_bin2->GetParameter(2)) /** 0.5*/ * d_norm_y_bin2;

  cout << "d_pol_scale_p0_y_bin2 = "<< d_pol_scale_p0_y_bin2 << endl;
  cout << "d_pol_scale_p1_y_bin2 = "<< d_pol_scale_p1_y_bin2 << endl;
  cout << "d_pol_scale_p2_y_bin2 = "<< d_pol_scale_p2_y_bin2 << endl;

  tf1_pol_scale_y_bin2 -> SetParameter(0,d_pol_scale_p0_y_bin2);
  tf1_pol_scale_y_bin2 -> SetParameter(1,d_pol_scale_p1_y_bin2);
  tf1_pol_scale_y_bin2 -> SetParameter(2,d_pol_scale_p2_y_bin2);
  TCanvas * TC_inv_bg_y_bin2 = new TCanvas("TC_inv_bg_y_bin2","TC_inv_bg_y_bin2",1280,720);
  TC_inv_bg_y_bin2->cd();
  h_TPC2_TOF3_invM_y_bin2->Rebin();

  h_TPC2_TOF3_invM_y_bin2->Draw();
  h_mx_TPC2_TOF3_invM_y_bin2->Scale(/*0.5**/d_norm_y_bin2);
  h_mx_TPC2_TOF3_invM_y_bin2->SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin2-> SetFillStyle(3001);

  h_mx_TPC2_TOF3_invM_y_bin2->Draw("same");
  // TC_inv_y_bin2->cd();
  tf1_pol_scale_y_bin2 -> SetLineColor(kBlue);
  tf1_pol_scale_y_bin2->Draw("same");

  // y_bin3
  TF1 * tf1_pol_scale_y_bin3 = new TF1("tf1_pol_scale_y_bin3","[0] + [1]*x + [2]*x**2",(d_mean_y_bin3 - (3*d_sigma_y_bin3)),(d_mean_y_bin3 + (3*d_sigma_y_bin3)));

  double  d_pol_scale_p0_y_bin3 = (tf1_pol_y_bin3->GetParameter(0)) /** 0.5*/ * d_norm_y_bin3;
  double  d_pol_scale_p1_y_bin3 = (tf1_pol_y_bin3->GetParameter(1)) /** 0.5*/ * d_norm_y_bin3;
  double  d_pol_scale_p2_y_bin3 = (tf1_pol_y_bin3->GetParameter(2)) /** 0.5*/ * d_norm_y_bin3;

  cout << "d_pol_scale_p0_y_bin3 = "<< d_pol_scale_p0_y_bin3 << endl;
  cout << "d_pol_scale_p1_y_bin3 = "<< d_pol_scale_p1_y_bin3 << endl;
  cout << "d_pol_scale_p2_y_bin3 = "<< d_pol_scale_p2_y_bin3 << endl;

  tf1_pol_scale_y_bin3 -> SetParameter(0,d_pol_scale_p0_y_bin3);
  tf1_pol_scale_y_bin3 -> SetParameter(1,d_pol_scale_p1_y_bin3);
  tf1_pol_scale_y_bin3 -> SetParameter(2,d_pol_scale_p2_y_bin3);
  TCanvas * TC_inv_bg_y_bin3 = new TCanvas("TC_inv_bg_y_bin3","TC_inv_bg_y_bin3",1280,720);
  TC_inv_bg_y_bin3->cd();
  h_TPC2_TOF3_invM_y_bin3->Rebin();

  h_TPC2_TOF3_invM_y_bin3->Draw();
  h_mx_TPC2_TOF3_invM_y_bin3->Scale(/*0.5**/d_norm_y_bin3);
  h_mx_TPC2_TOF3_invM_y_bin3->SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin3-> SetFillStyle(3001);

  h_mx_TPC2_TOF3_invM_y_bin3->Draw("same");
  // TC_inv_y_bin3->cd();
  tf1_pol_scale_y_bin3 -> SetLineColor(kBlue);
  tf1_pol_scale_y_bin3->Draw("same");

  // y_bin4
  TF1 * tf1_pol_scale_y_bin4 = new TF1("tf1_pol_scale_y_bin4","[0] + [1]*x + [2]*x**2",(d_mean_y_bin4 - (3*d_sigma_y_bin4)),(d_mean_y_bin4 + (3*d_sigma_y_bin4)));

  double  d_pol_scale_p0_y_bin4 = (tf1_pol_y_bin4->GetParameter(0)) /** 0.5*/ * d_norm_y_bin4;
  double  d_pol_scale_p1_y_bin4 = (tf1_pol_y_bin4->GetParameter(1)) /** 0.5*/ * d_norm_y_bin4;
  double  d_pol_scale_p2_y_bin4 = (tf1_pol_y_bin4->GetParameter(2)) /** 0.5*/ * d_norm_y_bin4;

  cout << "d_pol_scale_p0_y_bin4 = "<< d_pol_scale_p0_y_bin4 << endl;
  cout << "d_pol_scale_p1_y_bin4 = "<< d_pol_scale_p1_y_bin4 << endl;
  cout << "d_pol_scale_p2_y_bin4 = "<< d_pol_scale_p2_y_bin4 << endl;

  tf1_pol_scale_y_bin4 -> SetParameter(0,d_pol_scale_p0_y_bin4);
  tf1_pol_scale_y_bin4 -> SetParameter(1,d_pol_scale_p1_y_bin4);
  tf1_pol_scale_y_bin4 -> SetParameter(2,d_pol_scale_p2_y_bin4);
  TCanvas * TC_inv_bg_y_bin4 = new TCanvas("TC_inv_bg_y_bin4","TC_inv_bg_y_bin4",1280,720);
  TC_inv_bg_y_bin4->cd();
  h_TPC2_TOF3_invM_y_bin4->Rebin();

  h_TPC2_TOF3_invM_y_bin4->Draw();
  h_mx_TPC2_TOF3_invM_y_bin4->Scale(/*0.5**/d_norm_y_bin4);
  h_mx_TPC2_TOF3_invM_y_bin4->SetFillColor(kYellow);
  h_mx_TPC2_TOF3_invM_y_bin4-> SetFillStyle(3001);

  h_mx_TPC2_TOF3_invM_y_bin4->Draw("same");
  // TC_inv_y_bin4->cd();
  tf1_pol_scale_y_bin4 -> SetLineColor(kBlue);
  tf1_pol_scale_y_bin4->Draw("same");

  // END Scale the fitting function

  // Get (bg/(sig+bg))*Pol2 fitting function
  TF1 * tf1_bg_times_pol2 = new TF1("tf1_bg_times_pol2",bg_func,/*(d_mean - (5*d_sigma)),(d_mean + (5*d_sigma))*/0.99,1.1,/*1*//*2*/3/*4*/);
  TCanvas * TC_v1_invM = new TCanvas("TC_v1_invM","TC_v1_invM",1280,720);
  TC_v1_invM->cd();
  tp_phi_v1_vs_invM_pfx->Draw();
  tp_phi_v1_vs_invM_pfx->Fit(tf1_bg_times_pol2,"E","R",/*(d_mean - (5*d_sigma)),(d_mean + (5*d_sigma))*/0.99,1.1);

  TF1 * tf1_sig_bg_total = new TF1("tf1_sig_bg_total",sig_bg_func,/*(d_mean - (5*d_sigma)),(d_mean + (5*d_sigma))*/0.99,1.1,/*2*//*3*/4/*5*/);
  double d_bg_p0 = tf1_bg_times_pol2->GetParameter(0);
  double d_bg_p1 = tf1_bg_times_pol2->GetParameter(1);
  double d_bg_p2 = tf1_bg_times_pol2->GetParameter(2);
  // double d_bg_p3 = tf1_bg_times_pol2->GetParameter(3);

  cout << "d_bg_p0 = "<< d_bg_p0 <<endl;
  cout << "d_bg_p1 = "<< d_bg_p1 <<endl;
  cout << "d_bg_p2 = "<< d_bg_p2 <<endl;
  // cout << "d_bg_p3 = "<< d_bg_p3 <<endl;


  tf1_sig_bg_total->SetParameter(0,d_bg_p0);
  tf1_sig_bg_total->SetParameter(1,d_bg_p1);
  tf1_sig_bg_total->SetParameter(2,d_bg_p2);
  // tf1_sig_bg_total->SetParameter(3,d_bg_p3);

  // tf1_sig_bg_total -> SetParLimits(3,-1,-0.0000000001);

  tf1_sig_bg_total -> SetLineColor(kBlue);
  tp_phi_v1_vs_invM_pfx->Fit(tf1_sig_bg_total,"E+","R",/*(d_mean - (5*d_sigma)),(d_mean + (5*d_sigma))*/0.99,1.1);
  double d_v1_sig = tf1_sig_bg_total->GetParameter(/*1*//*2*/3/*4*/);
  double d_v1_sig_err = tf1_sig_bg_total->GetParError(/*1*//*2*/3/*4*/);

  // Extended fit
  // TF1 * tf1_sig_bg_total_extend = new TF1("tf1_sig_bg_total_extend",sig_bg_func,0.99,1.1,4);
  //
  // double d_p0_v1_sig = tf1_sig_bg_total->GetParameter(0);
  // double d_p1_v1_sig = tf1_sig_bg_total->GetParameter(1);
  // double d_p2_v1_sig = tf1_sig_bg_total->GetParameter(2);
  //
  // tf1_sig_bg_total_extend->SetParameter(0,d_p0_v1_sig);
  // tf1_sig_bg_total_extend->SetParameter(1,d_p1_v1_sig);
  // tf1_sig_bg_total_extend->SetParameter(2,d_p2_v1_sig);
  // tf1_sig_bg_total_extend->SetParameter(3,d_v1_sig);
  // tf1_sig_bg_total_extend -> SetLineColor(kRed);
  // tf1_sig_bg_total_extend->Draw("same");
  // End Extended fit

  TPaveText * ptext = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext -> AddText(Form("v1^{sig}: %.4f",d_v1_sig));
  ptext -> AddText(Form("v1^{sig} Error: %.4f",d_v1_sig_err));
  ptext->Draw("same");

  // bin2
  TF1 * tf1_bg_times_pol2_y_bin2 = new TF1("tf1_bg_times_pol2_y_bin2",bg_func_y_bin2,/*(d_mean_y_bin2 - (5*d_sigma_y_bin2)),(d_mean_y_bin2 + (5*d_sigma_y_bin2))*/0.99,1.1,/*1*//*2*/3/*4*/);
  TCanvas * TC_v1_invM_y_bin2 = new TCanvas("TC_v1_invM_y_bin2","TC_v1_invM_y_bin2",1280,720);
  TC_v1_invM_y_bin2->cd();
  tp_phi_v1_vs_invM_pfx_y_bin2->Draw();
  tp_phi_v1_vs_invM_pfx_y_bin2->Fit(tf1_bg_times_pol2_y_bin2,"E","R",/*(d_mean_y_bin2 - (5*d_sigma_y_bin2)),(d_mean_y_bin2 + (5*d_sigma_y_bin2))*/0.99,1.1);

  TF1 * tf1_sig_bg_total_y_bin2 = new TF1("tf1_sig_bg_total_y_bin2",sig_bg_func_y_bin2,/*(d_mean_y_bin2 - (5*d_sigma_y_bin2)),(d_mean_y_bin2 + (5*d_sigma_y_bin2))*/0.99,1.1,/*2*//*3*/4/*5*/);

  double d_bg_p0_y_bin2 = tf1_bg_times_pol2_y_bin2->GetParameter(0);
  double d_bg_p1_y_bin2 = tf1_bg_times_pol2_y_bin2->GetParameter(1);
  double d_bg_p2_y_bin2 = tf1_bg_times_pol2_y_bin2->GetParameter(2);
  // double d_bg_p3_y_bin2 = tf1_bg_times_pol2_y_bin2->GetParameter(3);

  cout << "d_bg_p0_y_bin2 = "<< d_bg_p0_y_bin2 <<endl;
  cout << "d_bg_p1_y_bin2 = "<< d_bg_p1_y_bin2 <<endl;
  cout << "d_bg_p2_y_bin2 = "<< d_bg_p2_y_bin2 <<endl;
  // cout << "d_bg_p3_y_bin2 = "<< d_bg_p3_y_bin2 <<endl;

  tf1_sig_bg_total_y_bin2->SetParameter(0,d_bg_p0_y_bin2);
  tf1_sig_bg_total_y_bin2->SetParameter(1,d_bg_p1_y_bin2);
  tf1_sig_bg_total_y_bin2->SetParameter(2,d_bg_p2_y_bin2);
  // tf1_sig_bg_total_y_bin2->SetParameter(3,d_bg_p3_y_bin2);

  // tf1_sig_bg_total -> SetParLimits(3,-1,-0.0000000001);

  tf1_sig_bg_total_y_bin2 -> SetLineColor(kBlue);
  tp_phi_v1_vs_invM_pfx_y_bin2->Fit(tf1_sig_bg_total_y_bin2,"E+","R",/*(d_mean_y_bin2 - (5*d_sigma_y_bin2)),(d_mean_y_bin2 + (5*d_sigma_y_bin2))*/0.99,1.1);

  double d_v1_sig_y_bin2    = tf1_sig_bg_total_y_bin2->GetParameter(/*1*//*2*/3/*4*/);
  double d_v1_sig_err_y_bin2 = tf1_sig_bg_total_y_bin2->GetParError(/*1*//*2*/3/*4*/);

  // Extended fit
  // TF1 * tf1_sig_bg_total_extend_y_bin2 = new TF1("tf1_sig_bg_total_extend_y_bin2",sig_bg_func_y_bin2,0.99,1.1,4);
  //
  // double d_p0_v1_sig_y_bin2 = tf1_sig_bg_total_y_bin2->GetParameter(0);
  // double d_p1_v1_sig_y_bin2 = tf1_sig_bg_total_y_bin2->GetParameter(1);
  // double d_p2_v1_sig_y_bin2 = tf1_sig_bg_total_y_bin2->GetParameter(2);
  //
  // tf1_sig_bg_total_extend_y_bin2->SetParameter(0,d_p0_v1_sig_y_bin2);
  // tf1_sig_bg_total_extend_y_bin2->SetParameter(1,d_p1_v1_sig_y_bin2);
  // tf1_sig_bg_total_extend_y_bin2->SetParameter(2,d_p2_v1_sig_y_bin2);
  // tf1_sig_bg_total_extend_y_bin2->SetParameter(3,d_v1_sig_y_bin2);
  // tf1_sig_bg_total_extend_y_bin2 -> SetLineColor(kRed);
  // tf1_sig_bg_total_extend_y_bin2->Draw("same");
  // End Extended fit

  TPaveText * ptext_y_bin2 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin2 -> AddText(Form("v1^{sig}: %.4f",d_v1_sig_y_bin2));
  ptext_y_bin2 -> AddText(Form("v1^{sig} Error: %.4f",d_v1_sig_err_y_bin2));
  ptext_y_bin2->Draw("same");

  // bin3
  TF1 * tf1_bg_times_pol2_y_bin3 = new TF1("tf1_bg_times_pol2_y_bin3",bg_func_y_bin3,/*(d_mean_y_bin3 - (5*d_sigma_y_bin3)),(d_mean_y_bin3 + (5*d_sigma_y_bin3))*/0.99,1.1,/*1*//*2*/3/*4*/);
  TCanvas * TC_v1_invM_y_bin3 = new TCanvas("TC_v1_invM_y_bin3","TC_v1_invM_y_bin3",1280,720);
  TC_v1_invM_y_bin3->cd();
  tp_phi_v1_vs_invM_pfx_y_bin3->Draw();
  tp_phi_v1_vs_invM_pfx_y_bin3->Fit(tf1_bg_times_pol2_y_bin3,"E","R",/*(d_mean_y_bin3 - (5*d_sigma_y_bin3)),(d_mean_y_bin3 + (5*d_sigma_y_bin3))*/0.99,1.1);

  TF1 * tf1_sig_bg_total_y_bin3 = new TF1("tf1_sig_bg_total_y_bin3",sig_bg_func_y_bin3,/*(d_mean_y_bin3 - (5*d_sigma_y_bin3)),(d_mean_y_bin3 + (5*d_sigma_y_bin3))*/0.99,1.1,/*2*//*3*/4/*5*/);
  double d_bg_p0_y_bin3 = tf1_bg_times_pol2_y_bin3->GetParameter(0);
  double d_bg_p1_y_bin3 = tf1_bg_times_pol2_y_bin3->GetParameter(1);
  double d_bg_p2_y_bin3 = tf1_bg_times_pol2_y_bin3->GetParameter(2);
  // double d_bg_p3_y_bin3 = tf1_bg_times_pol2_y_bin3->GetParameter(3);

  cout << "d_bg_p0_y_bin3 = "<< d_bg_p0_y_bin3 <<endl;
  cout << "d_bg_p1_y_bin3 = "<< d_bg_p1_y_bin3 <<endl;
  cout << "d_bg_p2_y_bin3 = "<< d_bg_p2_y_bin3 <<endl;
  // cout << "d_bg_p3_y_bin3 = "<< d_bg_p3_y_bin3 <<endl;

  tf1_sig_bg_total_y_bin3->SetParameter(0,d_bg_p0_y_bin3);
  tf1_sig_bg_total_y_bin3->SetParameter(1,d_bg_p1_y_bin3);
  tf1_sig_bg_total_y_bin3->SetParameter(2,d_bg_p2_y_bin3);
  // tf1_sig_bg_total_y_bin3->SetParameter(3,d_bg_p3_y_bin3);

  //
  // tf1_sig_bg_total -> SetParLimits(3,-1,-0.0000000001);

  tf1_sig_bg_total_y_bin3 -> SetLineColor(kBlue);
  tp_phi_v1_vs_invM_pfx_y_bin3->Fit(tf1_sig_bg_total_y_bin3,"E+","R",/*(d_mean_y_bin3 - (5*d_sigma_y_bin3)),(d_mean_y_bin3 + (5*d_sigma_y_bin3))*/0.99,1.1);
  double d_v1_sig_y_bin3 = tf1_sig_bg_total_y_bin3->GetParameter(/*1*//*2*/3/*4*/);
  double d_v1_sig_err_y_bin3 = tf1_sig_bg_total_y_bin3->GetParError(/*1*//*2*/3/*4*/);

  // Extended fit
  // TF1 * tf1_sig_bg_total_extend_y_bin3 = new TF1("tf1_sig_bg_total_extend_y_bin3",sig_bg_func_y_bin3,0.99,1.1,4);
  //
  // double d_p0_v1_sig_y_bin3 = tf1_sig_bg_total_y_bin3->GetParameter(0);
  // double d_p1_v1_sig_y_bin3 = tf1_sig_bg_total_y_bin3->GetParameter(1);
  // double d_p2_v1_sig_y_bin3 = tf1_sig_bg_total_y_bin3->GetParameter(2);
  //
  // tf1_sig_bg_total_extend_y_bin3->SetParameter(0,d_p0_v1_sig_y_bin3);
  // tf1_sig_bg_total_extend_y_bin3->SetParameter(1,d_p1_v1_sig_y_bin3);
  // tf1_sig_bg_total_extend_y_bin3->SetParameter(2,d_p2_v1_sig_y_bin3);
  // tf1_sig_bg_total_extend_y_bin3->SetParameter(3,d_v1_sig_y_bin3);
  // tf1_sig_bg_total_extend_y_bin3 -> SetLineColor(kRed);
  // tf1_sig_bg_total_extend_y_bin3->Draw("same");
  // End Extended fit

  TPaveText * ptext_y_bin3 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin3 -> AddText(Form("v1^{sig}: %.4f",d_v1_sig_y_bin3));
  ptext_y_bin3 -> AddText(Form("v1^{sig} Error: %.4f",d_v1_sig_err_y_bin3));
  ptext_y_bin3->Draw("same");

  // bin4
  TF1 * tf1_bg_times_pol2_y_bin4 = new TF1("tf1_bg_times_pol2_y_bin4",bg_func_y_bin4,/*(d_mean_y_bin4 - (5*d_sigma_y_bin4)),(d_mean_y_bin4 + (5*d_sigma_y_bin4))*/0.99,1.1,/*1*//*2*/3/*4*/);
  TCanvas * TC_v1_invM_y_bin4 = new TCanvas("TC_v1_invM_y_bin4","TC_v1_invM_y_bin4",1280,720);
  TC_v1_invM_y_bin4->cd();
  tp_phi_v1_vs_invM_pfx_y_bin4->Draw();
  tp_phi_v1_vs_invM_pfx_y_bin4->Fit(tf1_bg_times_pol2_y_bin4,"E","R",/*(d_mean_y_bin4 - (5*d_sigma_y_bin4)),(d_mean_y_bin4 + (5*d_sigma_y_bin4))*/0.99,1.1);

  TF1 * tf1_sig_bg_total_y_bin4 = new TF1("tf1_sig_bg_total_y_bin4",sig_bg_func_y_bin4,/*(d_mean_y_bin4 - (5*d_sigma_y_bin4)),(d_mean_y_bin4 + (5*d_sigma_y_bin4))*/0.99,1.1,/*2*//*3*/4/*5*/);
  double d_bg_p0_y_bin4 = tf1_bg_times_pol2_y_bin4->GetParameter(0);
  double d_bg_p1_y_bin4 = tf1_bg_times_pol2_y_bin4->GetParameter(1);
  double d_bg_p2_y_bin4 = tf1_bg_times_pol2_y_bin4->GetParameter(2);
  // double d_bg_p3_y_bin4 = tf1_bg_times_pol2_y_bin4->GetParameter(3);

  cout << "d_bg_p0_y_bin4 = "<< d_bg_p0_y_bin4 <<endl;
  cout << "d_bg_p1_y_bin4 = "<< d_bg_p1_y_bin4 <<endl;
  cout << "d_bg_p2_y_bin4 = "<< d_bg_p2_y_bin4 <<endl;
  // cout << "d_bg_p3_y_bin4 = "<< d_bg_p3_y_bin4 <<endl;

  tf1_sig_bg_total_y_bin4->SetParameter(0,d_bg_p0_y_bin4);
  tf1_sig_bg_total_y_bin4->SetParameter(1,d_bg_p1_y_bin4);
  tf1_sig_bg_total_y_bin4->SetParameter(2,d_bg_p2_y_bin4);
  // tf1_sig_bg_total_y_bin4->SetParameter(3,d_bg_p3_y_bin4);

  // tf1_sig_bg_total -> SetParLimits(3,-1,-0.0000000001);

  tf1_sig_bg_total_y_bin4 -> SetLineColor(kBlue);
  tp_phi_v1_vs_invM_pfx_y_bin4->Fit(tf1_sig_bg_total_y_bin4,"E+","R",/*(d_mean_y_bin4 - (5*d_sigma_y_bin4)),(d_mean_y_bin4 + (5*d_sigma_y_bin4))*/0.99,1.1);
  double d_v1_sig_y_bin4 = tf1_sig_bg_total_y_bin4->GetParameter(/*1*//*2*/3/*4*/);
  double d_v1_sig_err_y_bin4 = tf1_sig_bg_total_y_bin4->GetParError(/*1*//*2*/3/*4*/);

  // Extended fit
  // TF1 * tf1_sig_bg_total_extend_y_bin4 = new TF1("tf1_sig_bg_total_extend_y_bin4",sig_bg_func_y_bin4,0.99,1.1,4);
  //
  // double d_p0_v1_sig_y_bin4 = tf1_sig_bg_total_y_bin4->GetParameter(0);
  // double d_p1_v1_sig_y_bin4 = tf1_sig_bg_total_y_bin4->GetParameter(1);
  // double d_p2_v1_sig_y_bin4 = tf1_sig_bg_total_y_bin4->GetParameter(2);
  //
  // tf1_sig_bg_total_extend_y_bin4->SetParameter(0,d_p0_v1_sig_y_bin4);
  // tf1_sig_bg_total_extend_y_bin4->SetParameter(1,d_p1_v1_sig_y_bin4);
  // tf1_sig_bg_total_extend_y_bin4->SetParameter(2,d_p2_v1_sig_y_bin4);
  // tf1_sig_bg_total_extend_y_bin4->SetParameter(3,d_v1_sig_y_bin4);
  // tf1_sig_bg_total_extend_y_bin4 -> SetLineColor(kRed);
  // tf1_sig_bg_total_extend_y_bin4->Draw("same");
  // End Extended fit

  TPaveText * ptext_y_bin4 = new TPaveText(0.2,0.65,0.40,0.9,"NDCARC");
  ptext_y_bin4 -> AddText(Form("v1^{sig}: %.4f",d_v1_sig_y_bin4));
  ptext_y_bin4 -> AddText(Form("v1^{sig} Error: %.4f",d_v1_sig_err_y_bin4));
  ptext_y_bin4->Draw("same");
  // END Get (bg/(sig+bg))*Pol2 fitting function

  TLine * TL_3sigma_bot = new TLine((d_mean - (3*d_sigma)),-1.0,(d_mean - (3*d_sigma)),1.0);
  TLine * TL_3sigma_top = new TLine((d_mean + (3*d_sigma)),-1.0,(d_mean + (3*d_sigma)),1.0);

  TL_3sigma_bot -> SetLineColor(kRed);
  TL_3sigma_top -> SetLineColor(kRed);
  TC_v1_invM->cd();
  TL_3sigma_bot -> Draw();
  TL_3sigma_top -> Draw();

  // _y_bin2
  TLine * TL_3sigma_bot_y_bin2 = new TLine((d_mean_y_bin2 - (3*d_sigma_y_bin2)),-1.0,(d_mean_y_bin2 - (3*d_sigma_y_bin2)),1.0);
  TLine * TL_3sigma_top_y_bin2 = new TLine((d_mean_y_bin2 + (3*d_sigma_y_bin2)),-1.0,(d_mean_y_bin2 + (3*d_sigma_y_bin2)),1.0);

  TL_3sigma_bot_y_bin2 -> SetLineColor(kRed);
  TL_3sigma_top_y_bin2 -> SetLineColor(kRed);
  TC_v1_invM_y_bin2->cd();
  TL_3sigma_bot_y_bin2 -> Draw();
  TL_3sigma_top_y_bin2 -> Draw();

  // _y_bin3
  TLine * TL_3sigma_bot_y_bin3 = new TLine((d_mean_y_bin3 - (3*d_sigma_y_bin3)),-1.0,(d_mean_y_bin3 - (3*d_sigma_y_bin3)),1.0);
  TLine * TL_3sigma_top_y_bin3 = new TLine((d_mean_y_bin3 + (3*d_sigma_y_bin3)),-1.0,(d_mean_y_bin3 + (3*d_sigma_y_bin3)),1.0);

  TL_3sigma_bot_y_bin3 -> SetLineColor(kRed);
  TL_3sigma_top_y_bin3 -> SetLineColor(kRed);
  TC_v1_invM_y_bin3->cd();
  TL_3sigma_bot_y_bin3 -> Draw();
  TL_3sigma_top_y_bin3 -> Draw();

  // _y_bin4
  TLine * TL_3sigma_bot_y_bin4 = new TLine((d_mean_y_bin4 - (3*d_sigma_y_bin4)),-1.0,(d_mean_y_bin4 - (3*d_sigma_y_bin4)),1.0);
  TLine * TL_3sigma_top_y_bin4 = new TLine((d_mean_y_bin4 + (3*d_sigma_y_bin4)),-1.0,(d_mean_y_bin4 + (3*d_sigma_y_bin4)),1.0);

  TL_3sigma_bot_y_bin4 -> SetLineColor(kRed);
  TL_3sigma_top_y_bin4 -> SetLineColor(kRed);
  TC_v1_invM_y_bin4->cd();
  TL_3sigma_bot_y_bin4 -> Draw();
  TL_3sigma_top_y_bin4 -> Draw();

  h_phi_v1_vs_y->SetBinContent(4,d_v1_sig_y_bin2);
  h_phi_v1_vs_y->SetBinContent(5,d_v1_sig_y_bin3);
  h_phi_v1_vs_y->SetBinContent(6,d_v1_sig_y_bin4);

  h_phi_v1_vs_y->SetBinError(4,d_v1_sig_err_y_bin2);
  h_phi_v1_vs_y->SetBinError(5,d_v1_sig_err_y_bin3);
  h_phi_v1_vs_y->SetBinError(6,d_v1_sig_err_y_bin4);

  TH1D  *h_phi_v1_vs_y_flip = new TH1D("h_phi_v1_vs_y_flip","flipped v1 of #phi meson VS rapidity",305,-3.045+1.52,0.005+1.52);
  h_phi_v1_vs_y_flip->GetXaxis()->SetTitle("y-y_{cm}");
  h_phi_v1_vs_y_flip->GetYaxis()->SetTitle("v1");

  h_phi_v1_vs_y_flip->SetBinContent(180,d_v1_sig_y_bin2);
  h_phi_v1_vs_y_flip->SetBinContent(230,d_v1_sig_y_bin3);
  h_phi_v1_vs_y_flip->SetBinContent(280,d_v1_sig_y_bin4);

  h_phi_v1_vs_y_flip->SetBinContent(126,-d_v1_sig_y_bin2);
  h_phi_v1_vs_y_flip->SetBinContent(76,-d_v1_sig_y_bin3);
  h_phi_v1_vs_y_flip->SetBinContent(26,-d_v1_sig_y_bin4);

  h_phi_v1_vs_y_flip->SetBinError(180,d_v1_sig_err_y_bin2);
  h_phi_v1_vs_y_flip->SetBinError(230,d_v1_sig_err_y_bin3);
  h_phi_v1_vs_y_flip->SetBinError(280,d_v1_sig_err_y_bin4);

  h_phi_v1_vs_y_flip->SetBinError(126,d_v1_sig_err_y_bin2);
  h_phi_v1_vs_y_flip->SetBinError(76,d_v1_sig_err_y_bin3);
  h_phi_v1_vs_y_flip->SetBinError(26,d_v1_sig_err_y_bin4);

  TCanvas * TC_v1_y = new TCanvas("TC_v1_y","TC_v1_y",1280,720);
  gStyle->SetErrorX(0);
  TC_v1_y->cd();
  h_phi_v1_vs_y->SetMarkerColor(kRed);
  h_phi_v1_vs_y->SetMarkerStyle(29);
  h_phi_v1_vs_y->SetMarkerSize(2);
  h_phi_v1_vs_y->Draw();

  TCanvas * TC_v1_y_flip = new TCanvas("TC_v1_y_flip","TC_v1_y_flip",1280,720);
  TC_v1_y_flip->cd();
  gStyle->SetErrorX(0);
  h_phi_v1_vs_y_flip->SetMarkerColor(kRed);
  h_phi_v1_vs_y_flip->SetMarkerStyle(29);
  h_phi_v1_vs_y_flip->SetMarkerSize(2);
  TF1  *f0 = new TF1("f0","[0] + [1]*x",-1.53,1.53);
  h_phi_v1_vs_y_flip -> Fit(f0,"E","R",-1.53,1.53);

  h_phi_v1_vs_y_flip->Draw();

  TCanvas * TC_v1_y_overlap = new TCanvas("TC_v1_y_overlap","TC_v1_y_overlap",1280,720);
  TC_v1_y_overlap->cd();
  gStyle->SetErrorX(0);

  // TH1D  *h_phi_v1_vs_y_77 = new TH1D("h_phi_v1_vs_y_77","flipped v1 of #phi meson VS rapidity 7.7GeV",305,-1.525,1.525);
  // h_phi_v1_vs_y_77->GetXaxis()->SetTitle("y-y_{cm}");
  // h_phi_v1_vs_y_77->GetYaxis()->SetTitle("v1");
  //
  // h_phi_v1_vs_y_77->SetMarkerColor(kOrange);
  // h_phi_v1_vs_y_77->SetMarkerStyle(28);
  // h_phi_v1_vs_y_77->SetMarkerSize(2);
  //
  // h_phi_v1_vs_y_77->SetBinContent(83,-0.0212881);
  // h_phi_v1_vs_y_77->SetBinContent(103,-0.00989985);
  // h_phi_v1_vs_y_77->SetBinContent(123,-0.00434102);
  // h_phi_v1_vs_y_77->SetBinContent(143,-0.0012152);
  //
  // h_phi_v1_vs_y_77->SetBinContent(163,0.00122807);
  // h_phi_v1_vs_y_77->SetBinContent(183,0.00385395);
  // h_phi_v1_vs_y_77->SetBinContent(203,0.00983117);
  // h_phi_v1_vs_y_77->SetBinContent(223,0.0250622);
  //
  // h_phi_v1_vs_y_77->SetBinError(83,0.00153729);
  // h_phi_v1_vs_y_77->SetBinError(103,0.000939931);
  // h_phi_v1_vs_y_77->SetBinError(123,0.000790407);
  // h_phi_v1_vs_y_77->SetBinError(143,0.000809816);
  //
  // h_phi_v1_vs_y_77->SetBinError(163,0.000797943);
  // h_phi_v1_vs_y_77->SetBinError(183,0.000770219);
  // h_phi_v1_vs_y_77->SetBinError(203,0.000910965);
  // h_phi_v1_vs_y_77->SetBinError(223,0.0014815);
  //
  // h_phi_v1_vs_y_77->Draw();


  h_phi_v1_vs_y_flip->Draw("same");

  TCanvas * TC_resolution = new TCanvas("TC_resolution","TC_resolution",1280,720);
  TC_resolution->cd();
  // gStyle->SetErrorX(0);
  TH1D *h_resolution2 = new TH1D("h_resolution2","Event plane resolution",7,0.,35);
  h_resolution2->GetXaxis()->SetTitle("Centrality (%)");
  h_resolution2->GetYaxis()->SetTitle("Resolution ");
  h_resolution2->SetMarkerSize(2);
  h_resolution2->SetMarkerStyle(20);
  h_resolution2->SetBinContent(1,0.201086);
  h_resolution2->SetBinContent(2,0.312115);
  h_resolution2->SetBinContent(3,0.398419);
  h_resolution2->SetBinContent(4,0.454645);
  h_resolution2->SetBinContent(5,0.498949);
  h_resolution2->SetBinContent(6,0.596082);

  // h_resolution2->SetBinTitle(1,"5");
  // h_resolution2->SetBinTitle(2,"10");
  // h_resolution2->SetBinTitle(3,"15");
  // h_resolution2->SetBinTitle(4,"20");
  // h_resolution2->SetBinTitle(5,"25");
  // h_resolution2->SetBinTitle(6,"30");
  // h_resolution2->SetBinTitle(7,"> 30");

  h_resolution2->Draw();

  tf_out->cd();
  h_mx_TPC2_TOF3_invM->Write();
  h_mx_TPC2_TOF3_invM_y_bin2->Write();
  h_mx_TPC2_TOF3_invM_y_bin3->Write();
  h_mx_TPC2_TOF3_invM_y_bin4->Write();
  h_TPC2_TOF3_invM->Write();
  h_TPC2_TOF3_invM_y_bin2->Write();
  h_TPC2_TOF3_invM_y_bin3->Write();
  h_TPC2_TOF3_invM_y_bin4->Write();
  h_phi_v1_vs_y->Write();
  h_phi_v1_vs_y_flip->Write();

  h_mx_TPC2_TOF3_invM_same->Write();
  h_mx_TPC2_TOF3_invM_rebin->Write();
  h_mx_TPC2_TOF3_invM_same_y_bin2->Write();
  h_mx_TPC2_TOF3_invM_rebin_y_bin2->Write();
  h_mx_TPC2_TOF3_invM_rebin_y_bin3->Write();
  h_mx_TPC2_TOF3_invM_rebin_y_bin4->Write();

  tp_phi_v1_vs_invM_pfx->Write();
  tp_phi_v1_vs_invM_pfx_y_bin2->Write();
  tp_phi_v1_vs_invM_pfx_y_bin3->Write();
  tp_phi_v1_vs_invM_pfx_y_bin4->Write();
  TC_v1_invM->Write();
  TC_v1_invM_y_bin2->Write();
  TC_v1_invM_y_bin3->Write();
  TC_v1_invM_y_bin4->Write();
}
