#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>

#include "TMath.h"
#include "TGraph.h"
#include "TFrame.h"
#include "TMultiGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TRandom2.h"

using namespace std;



//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a Skewed Gaussian function for the fits
//



Double_t skewgaussian(Double_t *x, Double_t *par){
  double alpha = par[0], p = par[1], sigma = par[2], mu = par[3], I = par[4];
  
  double G = p/(sqrt(2*TMath::Pi())*sigma)*exp(-0.5*((x[0]-mu)/sigma)*((x[0]-mu)/sigma));
  double Corr = (1-p)*alpha/(2*sigma)*exp(-alpha*alpha/4)*exp(alpha*(x[0]-mu)/sigma) * (1 + TMath::Erf(-(x[0]-mu)/sigma));
  
  return I*(G+Corr);
}

TH1D * HistoCalibration(TH1D *h, double a, double b){
  int Nbin = h->GetNbinsX();
  cout << Nbin << endl;
  auto hcalib = new TH1D("","",Nbin*a/2,b,Nbin*a+b);
  
  for(int i=0;i<Nbin;i++){
    for(int k=0;k<h->GetBinContent(i);k++){
      hcalib->Fill(a*i+b);
    }
  }
  
  return(hcalib);
}
//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a global fit function for the skewed gaussian. This enable the fit of a whole
// spectrum with a series of skewed gaussian peaks and a quadratic background.
//
// The first parameter should be the number of peaks and should be fixed.
// The parameters 1, 2 and 3 are the parameters of the quadratic background.
// Then, for each peak there are 5 parameters :
//     _ a deformation parameter;
//     _ a parameter encoding the relative importance of the deformed gaussian and the regular one;
//     _ the mean of the gaussian;
//     _ the standard deviation of the gaussian;
//     _ the integral of the whole function.
//



Double_t fitfunctionskewgaussian(Double_t *x, Double_t *par){
  int npeaks = par[0];
  Double_t fit = 0;
  for(int peak = 0; peak<npeaks;peak++){
    Double_t peak_par[5] = {par[1],par[2],par[3],par[4+peak*2],par[5+peak*2]};
    fit+=skewgaussian(x,peak_par);
  }
  return par[4+npeaks*2]*fit;
}



//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Function fitting a spectrum in a certain range with a sum of skewed gaussians and a quadratic
// background (not always necessary. Fix parameters to 0 when not needed). The results are stored in
// a file and the histogram is plotted with the whole fit function and each individual contribution.
//
// The range of the fit, the fit parameters etc should be carefully selected before using the function
// Details about the different paramaters can be found as comments in the function.
//
// fname :    The name of the file containing the histogram
//
// outfname : The name of the file where the peaks information are to be written.
//

//
int FitSpectraSkewGaussian(TH1D *h1, TCanvas * cFit, vector<vector<double>> *Peaks, string peakfile, int type, double start = 3000, double end = 3600){
  double emin = start, emax = end; // min and max define the plot range.
  // min search define the starting range for the
  // peak search function  
  
  double binwidth = h1->GetBinWidth(1);
  
  ifstream fin(peakfile.c_str(), ios::in | ios::binary);
  
  string line;
  stringstream linestream;
  int npeaks = 0;
  bool UsePeak;
  double pos;
  
  vector<double> peakspos;
  vector<bool> peaksuse;
  bool first = true;
  while(fin.good()){
    getline(fin,line);
    linestream.clear();
    linestream.str(line);
    if(linestream.fail() or line.empty()) continue;
    if(line.find("#") != string::npos) continue;
    
    linestream >> pos >> UsePeak;
    if(pos<start or pos>end) continue;
    cout << start << " " << end << endl;
    cout << pos << endl;
    peakspos.push_back(pos);
    peaksuse.push_back(UsePeak);
    npeaks++;
  }
  
  fin.close();
  
  
  h1->GetXaxis()->SetRangeUser(emin, emax);
  
  
  // definition of the fit function and of an individual Skew Gaussian for each peak. The individual
  // Skew Gaussian is only used for plotting
  vector <TF1*> SkewGaussian;
  TF1 *FitFunction = new TF1("FitFunction", fitfunctionskewgaussian, emin, emax, 5+npeaks*2);
  
  // peakpos and peakmax will store the peak positions given by the search function.
  // sigmamin and sigmamax are the boundaries of the allowed range for the sigma of the gaussian
  // sigma is the initial sigma value
  // counts is to be computed using the peakmax and the initial sigma value. It will be used for
  // initialisation but also for integral boundaries definition
  // width is half of the range allowed to move the peak position
  double peakpos, peakmax, sigma = 5, sigmamin = 0, sigmamax = 40, counts, width = sigma;
  
  FitFunction->FixParameter(0, npeaks); // parameter 0 is the number of peaks. It is given by the
  // search function and must be fixed
  
  FitFunction->SetNpx((emax-emin)*10);
  
  
  FitFunction->SetParameter(1, 0.3); // alpha parameter
  FitFunction->SetParLimits(1, 0.01, 3); // The current boundaries should work
  //FitFunction->FixParameter(1,1);
  
  FitFunction->SetParameter(2, 1); // p parameter
  FitFunction->SetParLimits(2, 0.4, 1); // by definition p should be between 0 and 1
  
  //FitFunction->FixParameter(2,1);
  
  FitFunction->SetParameter(3, sigma); // sigma parameter
  //FitFunction->SetParLimits(3, sigmamin, sigmamax); //
  
  double prop[3];
  
  if(type==0){
    prop[0] = 0.1194;
    prop[1] = 0.1711;
    prop[2] = 0.7077;
  }
  else if(type==1){
    prop[0] = 0.0166;
    prop[1] = 0.1310;
    prop[2] = 0.8480;
  }
  else if(type==2){
    prop[0] = 0.2310;
    prop[1] = 0.7690;
    prop[2] = 0.0000;
  }
  
  stringstream ss;
  // creating all the Skew Gaussians
  for(int i=0; i<npeaks; i++){
    ss.str("");
    ss << "SkewGaussian" << i;
    SkewGaussian.push_back(new TF1(ss.str().c_str(),skewgaussian, emin, emax, 5));
    SkewGaussian[i]->SetParNames("alpha", "p", "Std_Dev", "Mean", "Counts");
    SkewGaussian[i]->SetNpx(10000);
  }
  
  
  TFitResultPtr fitRes;
  
  cout << "Npeaks : " << npeaks << endl;
  counts = 0;
  
  for(int peak=0; peak<npeaks; peak++){
    
    peakpos = peakspos[peak];
    cout << peakpos << endl;
    peakmax = h1->GetBinContent(h1->FindBin(peakpos));
    counts += peakmax * sigma * TMath::Sqrt(2*TMath::Pi());
    
    FitFunction->SetParameter(4 + peak * 2, peakpos); // position of the peak
    FitFunction->SetParLimits(4 + peak * 2, peakpos-10, peakpos+10);
    
    FitFunction->FixParameter(5 + peak * 2, prop[peak]); // Integral
    //FitFunction->SetParLimits(5 + peak * 2, 0.01*counts, 100*counts);
    
  }
  
  
  FitFunction->SetParameter(4 + 2 * npeaks, counts);
  
  // fitting the histogram
  fitRes = h1->Fit("FitFunction", "SMQE", "SAME", emin, emax);
  
  
  cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
  
  
  
  
  //cout << endl << "Bin width : " << binwidth << " keV" << endl;
  
  //cout << endl << "Fit function integral : " << FitFunction->Integral(emin,emax)/binwidth << endl;
  //cout << "Spectrum integral : " << h1->Integral(emin/binwidth,emax/binwidth) << endl;
  //cout << "Ratio : " << FitFunction->Integral(emin,emax)/(binwidth*h1->Integral(emin/binwidth,emax/binwidth)) << endl;
  
  double alpha = fitRes->Parameter(1), p = fitRes->Parameter(2), STD = fitRes->Parameter(3);
  double alphaerr = fitRes->ParError(1), perr = fitRes->ParError(2), STDerr = fitRes->ParError(3);
  /*
  cout << "alpha : " << alpha << " +/- " << alphaerr << endl;
  cout << "p : " << p << " +/- " << perr << endl;
  cout << "Sigma : " << STD << " +/- " << STDerr << endl;
  */
  
  // plotting each individual skew gaussian and printing the peaks informations
  for(int peak = 0; peak<npeaks; peak++){
    double Mean = fitRes->Parameter(4 + peak * 2), Counts = prop[peak]*fitRes->Parameter(4 + npeaks * 2)/binwidth;
    double Meanerr = fitRes->ParError(4 + peak * 2), Countserr = prop[peak]*fitRes->ParError(4 + npeaks * 2)/binwidth;
    
    SkewGaussian[peak]->SetParameter(0, alpha);
    SkewGaussian[peak]->SetParameter(1, p);
    SkewGaussian[peak]->SetParameter(2, STD);
    SkewGaussian[peak]->SetParameter(3, Mean);
    SkewGaussian[peak]->SetParameter(4, Counts*binwidth);
    SkewGaussian[peak]->SetLineColor(3);
    SkewGaussian[peak]->Draw("SAME");
    
    /*
    cout << "////////////////////////////////////////" << endl;
    cout << "Peak number : " << peak << endl;
    cout << "Peak Position : " << Mean << " +/- " << Meanerr << endl << endl;
    cout << "alpha : " << alpha << " +/- " << alphaerr << endl;
    cout << "p : " << p << " +/- " << perr << endl;
    cout << "Mean : " << Mean << " +/- " << Meanerr << endl;
    cout << "Sigma : " << STD << " +/- " << STDerr << endl;
    cout << "Counts : " << Counts << " +/- " << Countserr << endl;
    */
    
    vector<double> Peak;
    Peak.push_back(Mean);
    Peak.push_back(Meanerr);
    Peak.push_back(STD);
    Peak.push_back(STDerr);
    Peak.push_back(Counts);
    Peak.push_back(Countserr);
    Peak.push_back(peaksuse[peak]);
    Peaks->push_back(Peak);
    
  }
  
  // drawing background and fit function
  FitFunction->Draw("SAME");
  cFit->Update();
  //cFit->WaitPrimitive();
  
  //cFit->Close();
  //h1->Delete();
  
  return npeaks;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------
// This is a sort function for the sort method. It enables a sorting by decreasing energies
//
//
// peak1 :        peak 1
//
// peak2 :        peak 2
//
//



bool sorting_peaks_energy (vector<double> peak1, vector<double> peak2){
  if(peak1[0]<peak2[0]) return true;
  else return false;
}


int FitAllSkewGaussian(int Nfoil, int rebinwidth=1, bool V = false){
  TH1D * hfoil = new TH1D(), * hout = new TH1D();
  
  double min1 = 3000, max1 = 3600;
  double min2 = 3000, max2 = 3600;
  double min3 = 3000, max3 = 3600;
  
  double s = 5, ss = 0.005;
  string peakfile1, peakfile2;
  
  if(Nfoil==17){
    TFile *f = new TFile("Scan_FE_all_17.root");
    string hnames[13] = {"spec0", "spec1", "spec2","spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 17
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<10;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=10;i<13;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_17_0.txt";
    peakfile2 = "Peak_pos/C_17_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==15){
    TFile *f = new TFile("Scan_FE_all_15.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"}; // Foil 15
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<8;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=8;i<11;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_15_0.txt";
    peakfile2 = "Peak_pos/C_15_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==11){
    TFile *f = new TFile("Scan_FE_all_11.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec10", "spec11", "spec12"}; // Foil 11
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_11_0.txt";
    peakfile2 = "Peak_pos/C_11_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==10){
    TFile *f = new TFile("Scan_FE_all_10.root");
    string hnames[10] = {"spec0", "spec1", "spec2", "spec5", "spec6", "spec7", "spec8","spec9", "spec10", "spec11"}; // Foil 10
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<8;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=8;i<10;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_10_0.txt";
    peakfile2 = "Peak_pos/C_10_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==7){
    TFile *f = new TFile("Scan_FE_all_7.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 7
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_7_0.txt";
    peakfile2 = "Peak_pos/C_7_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==4){
    TFile *f = new TFile("Scan_FE_all_4.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 4
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<11;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    peakfile1 = "Peak_pos/C_4_0.txt";
    peakfile2 = "Peak_pos/C_4_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
  }
  else if(Nfoil==3){
    TFile *f = new TFile("Scan_FE_all_2_and_3.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10"}; // Foil 3
    
    hout->Add((TH1D*)f->Get(hnames[0].c_str()));
    hout->Add((TH1D*)f->Get(hnames[4].c_str()));
    hout->Add((TH1D*)f->Get(hnames[9].c_str()));
    hout->Add((TH1D*)f->Get(hnames[10].c_str()));
    
    hfoil->Add((TH1D*)f->Get(hnames[5].c_str()));
    hfoil->Add((TH1D*)f->Get(hnames[6].c_str()));
    hfoil->Add((TH1D*)f->Get(hnames[7].c_str()));
    hfoil->Add((TH1D*)f->Get(hnames[8].c_str()));
    
    peakfile1 = "Peak_pos/C_3_0.txt";
    peakfile2 = "Peak_pos/C_3_1.txt";
    
    min1 = 6240;
    max1 = 6400;
    
    min2 = 6600;
    max2 = 6800;
    
    min3 = 7040;
    max3 = 7200;
  }
  else if(Nfoil==2){
    TFile *f = new TFile("Scan_FE_all_2_and_3.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10"}; // Foil 2
    
    hout->Add((TH1D*)f->Get(hnames[0].c_str()));
    hout->Add((TH1D*)f->Get(hnames[4].c_str()));
    hout->Add((TH1D*)f->Get(hnames[9].c_str()));
    hout->Add((TH1D*)f->Get(hnames[10].c_str()));
    
    hfoil->Add((TH1D*)f->Get(hnames[1].c_str()));
    hfoil->Add((TH1D*)f->Get(hnames[2].c_str()));
    hfoil->Add((TH1D*)f->Get(hnames[3].c_str()));
    
    peakfile1 = "Peak_pos/C_2_0.txt";
    peakfile2 = "Peak_pos/C_2_1.txt";
    
    min1 = 6240;
    max1 = 6400;
    
    min2 = 6600;
    max2 = 6800;
    
    min3 = 7040;
    max3 = 7200;
  }
  else if(Nfoil==1){
    TFile *f = new TFile("Scan_FE_all_1.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 1
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    peakfile1 = "Peak_pos/C_1_0.txt";
    peakfile2 = "Peak_pos/C_1_1.txt";
    
    min1 = 3060;
    max1 = 3160;
    
    min2 = 3240;
    max2 = 3360;
    
    min3 = 3460;
    max3 = 3550;
    
  }
  else if(Nfoil==50){
    TFile *f = new TFile("Scan_FE_Si3N4_50nm.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec4", "spec5"}; // Foil Si3N4 50nm
    
    for(int i=0;i<1;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=1;i<5;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=5;i<6;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    peakfile1 = "Peak_pos/50nm_0.txt";
    peakfile2 = "Peak_pos/50nm_1.txt";
    
    min1 = 6270;
    max1 = 6410;
    
    min2 = 6620;
    max2 = 6810;
    
    min3 = 7070;
    max3 = 7200;
  }
  else if(Nfoil==31){
    TFile *f = new TFile("Scan_FE_Si3N4_30nm_1.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec4", "spec5"}; // Foil Si3N4 30nm
    
    for(int i=0;i<1;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=1;i<4;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=4;i<6;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    peakfile1 = "Peak_pos/30nm_1_0.txt";
    peakfile2 = "Peak_pos/30nm_1_1.txt";
    
    min1 = 6280;
    max1 = 6420;
    
    min2 = 6640;
    max2 = 6820;
    
    min3 = 7080;
    max3 = 7200;
  }
  else if(Nfoil==32){
    TFile *f = new TFile("Scan_FE_Si3N4_30nm_2.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec4", "spec5"}; // Foil Si3N4 30nm
    
    for(int i=0;i<1;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=1;i<4;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=4;i<6;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    peakfile1 = "Peak_pos/30nm_2_0.txt";
    peakfile2 = "Peak_pos/30nm_2_1.txt";
    
    min1 = 6280;
    max1 = 6420;
    
    min2 = 6640;
    max2 = 6820;
    
    min3 = 7080;
    max3 = 7200;
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  hout->Rebin(rebinwidth);
  hfoil->Rebin(rebinwidth);
  hout->SetTitle("All Data without the foil;Channel (ADC);Counts/keV");
  hfoil->SetTitle("All Data with the foil;Channel (ADC);Counts/keV");
  
  
  vector<vector<double>> Peaksout, Peaksfoil;
  if(V) cout << "Results out of the foil : " << endl;
  TCanvas *c1 =  new TCanvas("c1", "c1", 800, 600);
  FitSpectraSkewGaussian(hout,c1,&Peaksout,peakfile1,0,min1,max1);
  FitSpectraSkewGaussian(hout,c1,&Peaksout,peakfile1,1,min2,max2);
  FitSpectraSkewGaussian(hout,c1,&Peaksout,peakfile1,2,min3,max3);
  if(V) cout << "Results in the foil : " << endl;
  TCanvas *c2 =  new TCanvas("c2", "c2", 800, 600);
  FitSpectraSkewGaussian(hfoil,c2,&Peaksfoil,peakfile2,0,min1,max1);
  FitSpectraSkewGaussian(hfoil,c2,&Peaksfoil,peakfile2,1,min2,max2);
  FitSpectraSkewGaussian(hfoil,c2,&Peaksfoil,peakfile2,2,min3,max3);
  
  hout->GetXaxis()->SetRangeUser(min1,max3);
  hfoil->GetXaxis()->SetRangeUser(min1,max3);
  
  
  double Litt[8] = {5105.5,5144.3,5156.59,5388,5442.80,5485.56, 5762.64, 5804.77};
  double Litterr[8] = {0.8,0.8,0.14,1,0.13,0.12, 0.03, 0.05};
  
  double PeaksoutPos[8];
  double PeaksoutPoserr[8];
  double PeaksoutCounts[8];
  double PeaksoutCountserr[8];
  double FWHMout1 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][2], sFWHMout1 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][3], FWHMfoil1 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][2], sFWHMfoil1 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][3];
  double FWHMout2 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[3][2], sFWHMout2 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[3][3], FWHMfoil2 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[3][2], sFWHMfoil2 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[3][3];
  double FWHMout3 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[6][2], sFWHMout3 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[6][3], FWHMfoil3 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[6][2], sFWHMfoil3 = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[6][3];
  double PeaksfoilPos[8];
  double PeaksfoilPoserr[8];
  double PeaksfoilCounts[8];
  double PeaksfoilCountserr[8];
  double LittUse[8];
  double LitterrUse[8];
  bool UsePeak[8];
  int npeaksUse = 0;
  double SE[8] = {0.1798,0.1790,0.1785,0.1731,0.1718,0.1708,0.1650,0.1641}; // Stopping power carbon
  if(Nfoil==50 or Nfoil==31 or Nfoil == 32){ // Stopping power Si3N4
    SE[0] = 0.2300;
    SE[1] = 0.2292;
    SE[2] = 0.2286;
    SE[3] = 0.2222;
    SE[4] = 0.2207;
    SE[5] = 0.2195;
    SE[6] = 0.2126;
    SE[7] = 0.2116;
  }
  double SEUse[8];
  
  for(int i=0;i<8;i++){
    if(Peaksout[i][6]==false) continue;
    PeaksoutPos[npeaksUse] = Peaksout[i][0];
    PeaksoutPoserr[npeaksUse] = Peaksout[i][1];
    PeaksoutCounts[npeaksUse] = Peaksout[i][4];
    PeaksoutCountserr[npeaksUse] = Peaksout[i][5];
    PeaksfoilPos[npeaksUse] = Peaksfoil[i][0];
    PeaksfoilPoserr[npeaksUse] = Peaksfoil[i][1];
    PeaksfoilCounts[npeaksUse] = Peaksfoil[i][4];
    PeaksfoilCountserr[npeaksUse] = Peaksfoil[i][5];
    UsePeak[npeaksUse] = Peaksout[i][6];
    LittUse[npeaksUse] = Litt[i];
    LitterrUse[npeaksUse] = Litterr[i];
    SEUse[npeaksUse] = SE[i];
    npeaksUse++;
  }
  
  double residue[8];
  double residue_error[8];
  
  
  auto c3 = new TCanvas("c3", "c3", 200, 10, 700, 500);
  auto gr = new TGraphErrors(npeaksUse, PeaksoutPos, LittUse, PeaksoutPoserr, LitterrUse);
  
  c3->Divide(1,2);
  c3->cd(1);
  gr->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  TFitResultPtr FitRes = gr->Fit(line, "QES", "SAME E");
  TMatrixDSym cov = FitRes->GetCovarianceMatrix();
  double a = FitRes->Parameter(1), sa = FitRes->ParError(1), b = FitRes->Parameter(0), sb = FitRes->ParError(0), sab = cov[0][1];
  
  
  auto gr3 = new TGraphErrors(npeaksUse, PeaksfoilPos, LittUse, PeaksfoilPoserr, LitterrUse);
  TF1 *line3 = new TF1("Linear_regression3", "pol1", 0, 10000);
  TFitResultPtr FitRes3 = gr3->Fit(line3, "QES", "");

  
  
  auto c4 = new TCanvas("Calibrated histogram");
  TH1D *hcalibout = HistoCalibration(hout,a,b);
  TH1D *hcalibfoil = HistoCalibration(hfoil,FitRes3->Parameter(1),FitRes3->Parameter(0));
  
  hcalibout->SetLineColor(1);
  hcalibfoil->SetLineColor(2);
  hcalibout->Scale(hcalibfoil->GetEntries()/hcalibout->GetEntries());
  hcalibfoil->Draw("hist");
  hcalibout->Draw("hist same");
  
  
  
  for(int i = 0; i<npeaksUse; i++){
    double C = PeaksoutPos[i], sC = PeaksoutPoserr[i];
    residue[i] = a * C + b - LittUse[i];
    if(V) cout << "Residue n°" << i+1 << " : " << residue[i] << endl;
    residue_error[i] = TMath::Sqrt(a*a*sC*sC+sa*sa*C*C+sb*sb+2*C*sab);
  }
  
  
  TF1 *line2 = new TF1("Line", "pol0", 0, 10000);
  line2->SetParameter(0,0);
  line2->SetLineColor(2);
  auto gr2 = new TGraphErrors(npeaksUse, PeaksoutPos, residue, PeaksoutPoserr, residue_error);
  c3->cd(2);
  gr2->Draw("ap");
  line2->Draw("same");
  
  if(V){
    cout << " Chi2/NdF : " << FitRes->Chi2() << "/" << FitRes->Ndf() << " = " << FitRes->Chi2()/FitRes->Ndf() << " => Prob : " << FitRes->Prob() << endl << endl;
    
    cout << "Slope : " << a << " Slope error : " << sa << endl;
    cout << "Const : " << b << " Const error : " << sb << endl;
  }
  
  sFWHMout1 = TMath::Sqrt(a*a*sFWHMout1*sFWHMout1+sa*sa*FWHMout1*FWHMout1);
  FWHMout1 *= a;
  sFWHMfoil1 = TMath::Sqrt(a*a*sFWHMfoil1*sFWHMfoil1+sa*sa*FWHMfoil1*FWHMfoil1);
  FWHMfoil1 *= a;
  
  sFWHMout2 = TMath::Sqrt(a*a*sFWHMout2*sFWHMout2+sa*sa*FWHMout2*FWHMout2);
  FWHMout2 *= a;
  sFWHMfoil2 = TMath::Sqrt(a*a*sFWHMfoil2*sFWHMfoil2+sa*sa*FWHMfoil2*FWHMfoil2);
  FWHMfoil2 *= a;
  
  sFWHMout3 = TMath::Sqrt(a*a*sFWHMout3*sFWHMout3+sa*sa*FWHMout3*FWHMout3);
  FWHMout3 *= a;
  sFWHMfoil3 = TMath::Sqrt(a*a*sFWHMfoil3*sFWHMfoil3+sa*sa*FWHMfoil3*FWHMfoil3);
  FWHMfoil3 *= a;
  
  double dE[8];
  double sdE[8];
  double dx[8];
  double sdx[8];
  
  
  double sSE = 0.05;
  double w[8];
  double wE[8];
  double Sw = 0, mdx = 0, smdx, SwE=0, mdE=0, smdE, K1, K2;
  double meandx = 0, dxstd = 0;
  
  double stragg1 = TMath::Sqrt(FWHMfoil1*FWHMfoil1-FWHMout1*FWHMout1);
  double sstragg1 = TMath::Sqrt(TMath::Power(sFWHMfoil1*FWHMfoil1/stragg1,2)+TMath::Power(sFWHMout1*FWHMout1/stragg1,2));
  
  double stragg2 = TMath::Sqrt(FWHMfoil2*FWHMfoil2-FWHMout2*FWHMout2);
  double sstragg2 = TMath::Sqrt(TMath::Power(sFWHMfoil2*FWHMfoil2/stragg2,2)+TMath::Power(sFWHMout2*FWHMout2/stragg2,2));
  
  double stragg3 = TMath::Sqrt(FWHMfoil3*FWHMfoil3-FWHMout3*FWHMout3);
  double sstragg3 = TMath::Sqrt(TMath::Power(sFWHMfoil3*FWHMfoil3/stragg3,2)+TMath::Power(sFWHMout3*FWHMout3/stragg3,2));
  
  double stragg = (stragg1/(sstragg1*sstragg1)+stragg2/(sstragg2*sstragg2)+stragg3/(sstragg3*sstragg3))/(1/(sstragg1*sstragg1)+1/(sstragg2*sstragg2)+1/(sstragg3*sstragg3));
  double sstragg = 1/TMath::Sqrt(1/(sstragg1*sstragg1)+1/(sstragg2*sstragg2)+1/(sstragg3*sstragg3));
  
  double mstragg = (stragg1+stragg2+stragg3)/3;
  double stdstragg = TMath::Sqrt(((stragg1-mstragg)*(stragg1-mstragg)+(stragg2-mstragg)*(stragg2-mstragg)+(stragg3-mstragg)*(stragg3-mstragg))/3);
  
  K1 = TMath::Sqrt(13.52*13.52-10.9*10.9);
  K2 = TMath::Sqrt(TMath::Power(13.52*0.05/K1,2)+TMath::Power(10.9*0.08/K1,2));
  
  for(int i=0;i<npeaksUse;i++){
    double Cout = PeaksoutPos[i], sCout = PeaksoutPoserr[i], Cfoil = PeaksfoilPos[i], sCfoil = PeaksfoilPoserr[i];
    PeaksoutPos[i] = a*Cout+b;
    PeaksoutPoserr[i] = TMath::Sqrt(a*a*sCout*sCout+sa*sa*Cout*Cout+sb*sb+2*Cout*sab);
    PeaksfoilPos[i] = a*Cfoil+b;
    PeaksfoilPoserr[i] = TMath::Sqrt(a*a*sCfoil*sCfoil+sa*sa*Cfoil*Cfoil+sb*sb+2*Cfoil*sab);
    dE[i] = (Cout-Cfoil)*a;
    sdE[i] = TMath::Sqrt(a*a*(sCout*sCout+sCfoil*sCfoil)+(Cout-Cfoil)*(Cout-Cfoil)*sa*sa);
    dx[i] = dE[i]/SEUse[i];
    sdx[i] = dx[i]*TMath::Sqrt(sdE[i]*sdE[i]/(dE[i]*dE[i])+sSE*sSE);
    w[i] = 1/(sdx[i]*sdx[i]);
    wE[i] = 1/(sdE[i]*sdE[i]);
    Sw += w[i];
    SwE += wE[i];
    mdx += dx[i]*w[i];
    meandx += dx[i]/npeaksUse;
    mdE += dE[i]*wE[i];
    
    if(V){
      cout << "//////////////////////////////////////////////" << endl;
      cout << "Peak position in litterature : " << LittUse[i] << " keV" << endl;
      cout << "Measured peak position without the foil : " << PeaksoutPos[i] << " +/- " << PeaksoutPoserr[i] << " keV" << endl;
      cout << "Measured peak position with the foil : " << PeaksfoilPos[i] << " +/- " << PeaksfoilPoserr[i] << " keV" << endl;
      cout << "Measured peak integral without the foil : " << PeaksoutCounts[i] << " +/- " << PeaksoutCountserr[i] << endl;
      cout << "Measured peak integral with the foil : " << PeaksfoilCounts[i] << " +/- " << PeaksfoilCountserr[i] << endl;
      cout << "Delta E : " << dE[i] << " +/- " << sdE[i] << " keV" << endl;
      cout << "Delta x : " << dx[i] << " +/- " << sdx[i] << " nm" << endl << endl;
    }
  }
  
  for(int i=0;i<npeaksUse;i++) dxstd += (dx[i]-meandx)*(dx[i]-meandx)/npeaksUse;
  
  dxstd = TMath::Sqrt(dxstd);
  
  mdx =  mdx/Sw;
  smdx = 1/TMath::Sqrt(Sw);
  mdE = mdE/SwE;
  smdE = 1/TMath::Sqrt(SwE);
  
  
  if(V){
    cout << "//////////////////////////////////////////////" << endl;
    cout << "FWHM1 without the foil : " << FWHMout1 << " +/- " << sFWHMout1 << " keV" << endl;
    cout << "FWHM1 with the foil : " << FWHMfoil1 << " +/- " << sFWHMfoil1 << " keV" << endl;
    cout << "FWHM2 without the foil : " << FWHMout2 << " +/- " << sFWHMout2 << " keV" << endl;
    cout << "FWHM2 with the foil : " << FWHMfoil2 << " +/- " << sFWHMfoil2 << " keV" << endl;
    cout << "FWHM3 without the foil : " << FWHMout3 << " +/- " << sFWHMout3 << " keV" << endl;
    cout << "FWHM3 with the foil : " << FWHMfoil3 << " +/- " << sFWHMfoil3 << " keV" << endl << endl;
    
    cout << "Stragg1 : " << stragg1 << " +/- " << sstragg1 << endl;
    cout << "Stragg2 : " << stragg2 << " +/- " << sstragg2 << endl;
    cout << "Stragg3 : " << stragg3 << " +/- " << sstragg3 << endl;
  }
  
  
  cout << "//////////////////////////////////////////////" << endl;
  cout << "Weighted average delta E : " << mdE << " +/- " << smdE << " keV" << endl;
  cout << "Weighted average delta x : " << mdx << " +/- " << smdx << " nm" << endl;
  cout << "mean delta x : " << meandx << " +/- " << dxstd << " nm" << endl;
  cout << "Weighted average Straggling : " << stragg << " +/- " << sstragg << " keV" << endl;
  cout << "mean Straggling : " << mstragg << " +/- " << stdstragg << " keV" << endl;
  cout << "Simulated straggling : " << K1 << " +/- " << K2 << " keV" << endl;
  
  cout << endl << "Binwidth out : " << hout->GetBinWidth(0) << endl;
  cout << "Binwidth foil : " << hfoil->GetBinWidth(0) << endl;
  
  
  cout << stragg1 << " " << sstragg1 << " " << stragg2 << " " << sstragg2 << " " << stragg3 << " " << sstragg3;
  cout << " " << stragg << " " << sstragg << " " << mstragg << " " << stdstragg;
  cout << " " << mdx << " " << smdx << " " << meandx << " " << dxstd << endl;
  return 0;
}


int CheckAlignementSE(int Nfoil){
  
  TH1D * hstart = new TH1D();
  TH1D * hfoil = new TH1D();
  TH1D * hend = new TH1D();
  
  if(Nfoil==17){
    TFile *f = new TFile("Scan_FE_all_17.root");
    string hnames[13] = {"spec0", "spec1", "spec2","spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11","spec12"}; // Foil 17
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<10;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=10;i<13;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==15){
    TFile *f = new TFile("Scan_FE_all_15.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"}; // Foil 15
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<8;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=8;i<11;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==11){
    TFile *f = new TFile("Scan_FE_all_11.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec10", "spec11", "spec12"}; // Foil 11
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==10){
    TFile *f = new TFile("Scan_FE_all_10.root");
    string hnames[10] = {"spec0", "spec1", "spec2", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"}; // Foil 10
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<8;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=8;i<10;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==7){
    TFile *f = new TFile("Scan_FE_all_7.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 7
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==4){
    TFile *f = new TFile("Scan_FE_all_4.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 4
    /*
    hstart->Add((TH1D*)f->Get("spec0"));
    hfoil->Add((TH1D*)f->Get("spec4"));
    hend->Add((TH1D*)f->Get("spec4"));
    */
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<11;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==1){
    TFile *f = new TFile("Scan_FE_all_1.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 4
    
    
    for(int i=0;i<3;i++){
      hstart->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<11;i++){
      hend->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  hstart->GetXaxis()->SetRangeUser(3000, 3600);
  hfoil->GetXaxis()->SetRangeUser(3000, 3600);
  hend->GetXaxis()->SetRangeUser(3000, 3600);
  
  auto c1 = new TCanvas("c1","c1",700,500);
  hfoil->Draw("hist");
  hstart->SetLineColor(2);
  hend->SetLineColor(3);
  hstart->Draw("hist same");
  hend->Draw("hist same");
  c1->Update();
  return 0;
}

int CheckAlignementF(int Nfoil){
  
  vector<TH1D*> h;
  
  int colors[] = {1,2,3,4,6,7,8,9,10,11,12};
  if(Nfoil==17){
    TFile *f = new TFile("Scan_FE_all_17.root");
    string hnames[7] = {"spec3","spec4", "spec5", "spec6", "spec7", "spec8","spec9"}; // Foil 17
    
    for(int i=0;i<7;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==15){
    TFile *f = new TFile("Scan_FE_all_15.root");
    string hnames[5] = {"spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 15
    
    for(int i=0;i<5;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==11){
    TFile *f = new TFile("Scan_FE_all_11.root");
    string hnames[6] = {"spec3", "spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 11
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==10){
    TFile *f = new TFile("Scan_FE_all_10.root");
    string hnames[5] = {"spec5", "spec6", "spec7", "spec8","spec9"}; // Foil 10
    
    for(int i=0;i<5;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==7){
    TFile *f = new TFile("Scan_FE_all_7.root");
    string hnames[6] = {"spec4", "spec5", "spec6", "spec7", "spec8", "spec9"}; // Foil 7
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==4){
    TFile *f = new TFile("Scan_FE_all_4.root");
    string hnames[6] = {"spec5", "spec6", "spec7", "spec8", "spec9", "spec10"}; // Foil 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==1){
    TFile *f = new TFile("Scan_FE_all_1.root");
    string hnames[6] = {"spec4", "spec5", "spec6", "spec7", "spec8", "spec9"}; // Foil 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  double max = 0;
  for(int i=0;i<h.size();i++){
    if(h[i]->GetMaximum()>max) max = h[i]->GetMaximum();
  }
  
  for(int i=0;i<h.size();i++){
    h[i]->GetXaxis()->SetRangeUser(3000, 3600);
    h[i]->GetYaxis()->SetRangeUser(0, max*1.1);
  }
  
  auto c1 = new TCanvas("c1","c1",700,500);
  h[0]->SetLineColor(colors[0]);
  h[0]->Draw("hist");
  for(int i=1;i<h.size();i++){
    h[i]->SetLineColor(colors[i]);
    h[i]->Draw("hist same");
  }
  c1->Update();
  return 0;
}


int CheckAlignementAll(int Nfoil){
  
  vector<TH1D*> h;
  
  int colors[] = {1,2,3,4,6,7,8,9,10,11,12,13};
  
  if(Nfoil==17){
    TFile *f = new TFile("Scan_FE_all_17.root");
    string hnames[7] = {"spec3","spec4", "spec5", "spec6", "spec7", "spec8","spec9"}; // Foil 17
    
    for(int i=0;i<7;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==15){
    TFile *f = new TFile("Scan_FE_all_15.root");
    string hnames[5] = {"spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 15
    
    for(int i=0;i<5;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==11){
    TFile *f = new TFile("Scan_FE_all_11.root");
    string hnames[6] = {"spec3", "spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 11
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==10){
    TFile *f = new TFile("Scan_FE_all_10.root");
    string hnames[5] = {"spec5", "spec6", "spec7", "spec8","spec9"}; // Foil 10
    
    for(int i=0;i<5;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==7){
    TFile *f = new TFile("Scan_FE_all_7.root");
    string hnames[6] = {"spec4", "spec5", "spec6", "spec7", "spec8", "spec9"}; // Foil 7
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==4){
    TFile *f = new TFile("Scan_FE_all_4.root");
    string hnames[6] = {"spec5", "spec6", "spec7", "spec8", "spec9", "spec10"}; // Foil 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==1){
    TFile *f = new TFile("Scan_FE_all_1.root");
    string hnames[6] = {"spec4", "spec5", "spec6", "spec7", "spec8", "spec9"}; // Foil 1
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==2){
    TFile *f = new TFile("Scan_FE_all_2_and_3.root");
    string hnames[11] = {"spec0", "spec1", "spec2", "spec3","spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10"}; // Foil 2 and 3
    
    for(int i=0;i<11;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==50){
    TFile *f = new TFile("Scan_FE_Si3N4_50nm.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec3","spec4", "spec5"}; // Foil 2 and 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==31){
    TFile *f = new TFile("Scan_FE_Si3N4_30nm_1.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec3","spec4", "spec5"}; // Foil 2 and 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==32){
    TFile *f = new TFile("Scan_FE_Si3N4_30nm_2.root");
    string hnames[6] = {"spec0", "spec1", "spec2", "spec3","spec4", "spec5"}; // Foil 2 and 3
    
    for(int i=0;i<6;i++){
      h.push_back((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  double tmax = 0;
  for(int i=0;i<h.size();i++){
    if(h[i]->GetMaximum()>tmax) tmax = h[i]->GetMaximum();
  }
  
  for(int i=0;i<h.size();i++){
    //h[i]->GetXaxis()->SetRangeUser(3000, 3600);
    h[i]->GetYaxis()->SetRangeUser(0, tmax*1.1);
  }
  
  auto c1 = new TCanvas("c1","c1",700,500);
  h[0]->SetLineColor(colors[0]);
  h[0]->Draw("hist");
  for(int i=1;i<h.size();i++){
    h[i]->SetLineColor(colors[i]);
    h[i]->Draw("hist same");
  }
  c1->BuildLegend();
  c1->Update();
  return 0;
}




void plot(){
  double straggC[] = {12.2,8.1,8.7,11.5,11.5,10.8,9.4,10.6,10.9};
  double straggCerr[] = {0.8,0.7,0.5,0.9,1.1,0.9,0.9,0.9,1.0};
  double dxC[] = {151.5,102.4,99.4,132.3,148.8,122.0,88.3,90.2,112.3};
  double dxCerr[] = {3.2,2.2,2.1,3.1,3.4,3.0,2.3,2.5,2.7};
  
  double straggS[] = {10.5,7.8,9.2};
  double straggSerr[] = {0.5,0.7,0.7};
  double dxS[] = {28.2,12.3,18.7};
  double dxSerr[] = {1.2,1.1,1.1};
  
  auto c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);
  auto gr1 = new TGraphErrors(9, dxC, straggC, dxCerr, straggCerr);
  auto gr2 = new TGraphErrors(3, dxS, straggS, dxSerr, straggSerr);
  
  gr1->SetTitle("Carbon Foils;Thickness (nm);Straggling (keV)");
  
  gr1->GetXaxis()->SetRangeUser(0,160);
  
  c1->Divide(1,2);
  c1->cd(1);
  gr1->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  
  
  
  
  TFitResultPtr fitRes = gr1->Fit(line, "QES", "SAME E");
  TMatrixDSym cov = fitRes->GetCovarianceMatrix();
  double a = fitRes->Parameter(1), sa = fitRes->ParError(1), b = fitRes->Parameter(0), sb = fitRes->ParError(0), sab = cov[0][1];
  
  
  cout << " Chi2/NdF C : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl << endl;
  
  
  
  
  gr2->SetTitle("Si3N4 Foils;Thickness (nm);Straggling (keV)");
  
  c1->cd(2);
  gr2->Draw("ap");
  
  TF1 *line2 = new TF1("Linear_regression2", "pol1", 0, 10000);
  
  TFitResultPtr fitRes2 = gr2->Fit(line2, "QES", "SAME E");
  TMatrixDSym cov2 = fitRes2->GetCovarianceMatrix();
  double a2 = fitRes2->Parameter(1), sa2 = fitRes2->ParError(1), b2 = fitRes2->Parameter(0), sb2 = fitRes2->ParError(0), sab2 = cov2[0][1];
  
  
  cout << " Chi2/NdF S : " << fitRes2->Chi2() << "/" << fitRes2->Ndf() << " = " << fitRes2->Chi2()/fitRes2->Ndf() << " => Prob : " << fitRes2->Prob() << endl;
  
  return 0;
}
