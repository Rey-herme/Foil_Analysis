#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>

#include "Event.cpp"

#include "TGraph.h"
#include "TFrame.h"
#include "TMultiGraph.h"
#include "TGraph2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TRandom2.h"

using namespace std;


TTree * GetTree(string fname){
  TFile * f = new TFile(fname.c_str(), "read");
  TTree * tree = (TTree*)f->Get("TreeMaster");

  return tree;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a Skewed Gaussian function for the fits
//



Double_t skewgaussian2(Double_t *x, Double_t *par){
  //par[0] =  alpha // par[1] = p // par[2] = mu // par[3] = sigma // par[4] = C

  double tau = par[0], s = par[1], mu = par[2], A = par[3];
    
  double G = A/(2*tau)*TMath::Exp((x[0]-mu)/tau+s*s/(2*tau*tau))*TMath::Erfc(((x[0]-mu)/s+s/tau)/TMath::Sqrt(2));
 
  return G;
}


//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a Skewed Gaussian function for the fits
//



Double_t skewgaussian(Double_t *x, Double_t *par){
  //par[0] =  alpha // par[1] = p // par[2] = mu // par[3] = sigma // par[4] = C
    
  double G = par[1]/(sqrt(2*TMath::Pi())*par[3])*exp(-0.5*((x[0]-par[2])/par[3])*((x[0]-par[2])/par[3]));
  double Corr = (1-par[1])*par[0]/(2*par[3])*exp(-par[0]*par[0]/4)*exp(par[0]*(x[0]-par[2])/par[3]) * (1 + TMath::Erf(-(x[0]-par[2])/par[3]));
 
  return par[4]*(G+Corr);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a Skewed Gaussian function for the fits
//



Double_t gaussian(Double_t *x, Double_t *par){
  double s = par[0], c = par[1], A = par[2];
    
  double G = A*exp(-0.5*((x[0]-c)/s)*((x[0]-c)/s));
 
  return G;
}





//-------------------------------------------------------------------------------------------------------------------------------------------------
//
// Definition of a quadratic background for the fits
//


Double_t background(Double_t *x, Double_t *par){
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
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



Double_t fitfunctiongaussian(Double_t *x, Double_t *par){
  int npeaks = par[0];
  Double_t fit = 0;
  for(int peak = 0; peak<npeaks;peak++){
    Double_t peak_par[3] = {par[1],par[2+peak*2],par[3+peak*2]};
    fit+=gaussian(x,peak_par);
  }
  return fit;
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
  Double_t fit = background(x,&par[1]);
  for(int peak = 0; peak<npeaks;peak++){
    Double_t peak_par[5] = {par[4],par[5],par[7+peak*2],par[6],par[8+peak*2]};
    //Double_t peak_par[5] = {par[4],par[5],par[7+peak*2],par[6],par[8+peak*2]};
    fit+=skewgaussian(x,peak_par);
  }
  return fit;
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



Double_t fitfunctionskewgaussian2(Double_t *x, Double_t *par){
  int npeaks = par[0];
  Double_t fit = 0;
  for(int peak = 0; peak<npeaks;peak++){
    Double_t peak_par[4] = {par[1],par[2],par[3+peak*2],par[4+peak*2]};
    fit+=skewgaussian2(x,peak_par);
  }
  return fit;
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


int FitSpectraGaussian(TH1D *h1, TCanvas * cFit, vector<vector<double>> *Peaks, int search = 5, double s = 0.005, double start = 3000, double end = 3600){

  //cout << "////////////////////////////////////////" << endl;
  double min = start, max = end, minsearch = start; // min and max define the plot range.
                                                   // min search define the starting range for the
                                                   // peak search function

  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(minsearch), h1->GetXaxis()->FindFixBin(max));
  double binwidth = h1->GetBinWidth(1);


  // Peak search function
  TSpectrum *spec = new TSpectrum();
  spec->Search(h1, search, "", s); // the parameters of the search should be carefully examined.
  int npeaks=spec->GetNPeaks();
  //cout << "npeaks =  " << npeaks << endl;
  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(min), h1->GetXaxis()->FindFixBin(max));
  cFit->Update();
  double emin = start, emax = end;// Range of channels to be fitted

  // definition of the fit function and of an individual Skew Gaussian for each peak. The individual
  // Skew Gaussian is only used for plotting
  vector <TF1*> Gaussian;
  
  // peakpos and peakmax will store the peak positions given by the search function.
  // sigmamin and sigmamax are the boundaries of the allowed range for the sigma of the gaussian
  // sigma is the initial sigma value
  // counts is to be computed using the peakmax and the initial sigma value. It will be used for
  // initialisation but also for integral boundaries definition
  // width is half of the range allowed to move the peak position
  double peakpos, peakmax, sigma = 5, sigmamin = 3, sigmamax = 12, width = sigma;
  
  // creating a background function to plot it on the histogram
  stringstream ss;

  // creating all the Gaussians
  /*
  for(int i=0; i<npeaks; i++){
    
    ss.str("");
    ss << "Gaussian" << i;
    Gaussian.push_back(new TF1(ss.str().c_str(),gaussian, emin, emax, 3));
    Gaussian[i]->SetParNames("Std_Dev", "Mean", "Amp");
    Gaussian[i]->SetNpx(10000);
    }*/


  TFitResultPtr fitRes;
  vector<double> vSTD, vSTDerr, vMean, vMeanerr, vAmp, vAmperr;

  //cout << "Npeaks : " << npeaks << endl;

  for(int peak=0; peak<npeaks; peak++){

    
    peakpos = spec->GetPositionX()[peak];
    peakmax = spec->GetPositionY()[peak];
    ss.str("");
    ss << "Gaussian" << peak;
    Gaussian.push_back(new TF1(ss.str().c_str(),gaussian, peakpos-sigma, peakpos+2*sigma, 3));
    Gaussian[peak]->SetParNames("Std_Dev", "Mean", "Amp");
    Gaussian[peak]->SetNpx(10000);

    Gaussian[peak]->SetParameter(0, sigma); // position of the peak
    Gaussian[peak]->SetParLimits(0, sigmamin, sigmamax);
    
    Gaussian[peak]->SetParameter(1, peakpos); // position of the peak
    Gaussian[peak]->SetParLimits(1, peakpos-width, peakpos+width);
    
    Gaussian[peak]->SetParameter(2, peakmax); // Amplitude
    Gaussian[peak]->SetParLimits(2, 0.8*peakmax, 1.2*peakmax);
    

    fitRes = h1->Fit(ss.str().c_str(), "LSMQE", "SAME", peakpos-sigma, peakpos+2*sigma);
    /*
      cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
    */

    vSTD.push_back(fitRes->Parameter(0));
    vMean.push_back(fitRes->Parameter(1));
    vAmp.push_back(fitRes->Parameter(2));
    vSTDerr.push_back(fitRes->ParError(0));
    vMeanerr.push_back(fitRes->ParError(1));
    vAmperr.push_back(fitRes->ParError(2));

    
  }
  
  
  //cout << endl << "Bin width : " << binwidth << " keV" << endl;

  //cout << endl << "Fit function integral : " << FitFunction->Integral(emin,emax)/binwidth << endl;
  //cout << "Spectrum integral : " << h1->Integral(emin/binwidth,emax/binwidth) << endl;
  //cout << "Ratio : " << FitFunction->Integral(emin,emax)/(binwidth*h1->Integral(emin/binwidth,emax/binwidth)) << endl;


  // plotting each individual skew gaussian and printing the peaks informations
  for(int peak = 0; peak<npeaks; peak++){
    double STD = vSTD[peak], Mean = vMean[peak], Amp = vAmp[peak]/binwidth;
    double STDerr = vSTDerr[peak], Meanerr = vMeanerr[peak], Amperr = vAmperr[peak]/binwidth;

    Gaussian[peak]->SetParameter(0, STD);
    Gaussian[peak]->SetParameter(1, Mean);
    Gaussian[peak]->SetParameter(2, Amp*binwidth);
    Gaussian[peak]->SetLineColor(3);
    Gaussian[peak]->Draw("SAME");

    /*
    cout << "////////////////////////////////////////" << endl;
    cout << "Peak number : " << peak << endl;
    cout << "Sigma : " << STD << " +/- " << STDerr << endl << endl;
    cout << "Peak Position : " << Mean << " +/- " << Meanerr << endl;
    cout << "Counts : " << Amp*TMath::Sqrt(2*TMath::Pi())*STD << " +/- " << Amp*TMath::Sqrt(Amperr*Amperr/(Amp*Amp)+STDerr*STDerr/(STD*STD)) << endl << endl;
    */

    vector<double> Peak;
    Peak.push_back(Mean);
    Peak.push_back(Meanerr);
    Peak.push_back(STD);
    Peak.push_back(STDerr);
    Peaks->push_back(Peak);
  }
  
  cFit->Update();
  //cFit->WaitPrimitive();

  //cFit->Close();
  //h1->Delete();
  
  return npeaks;
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


int FitSpectraSkewGaussian(TH1D *h1, TCanvas * cFit, vector<vector<double>> *Peaks, int search = 5, double s = 0.005, double start = 3000, double end = 3600){

  //cout << "////////////////////////////////////////" << endl;
  double min = start, max = end, minsearch = start; // min and max define the plot range.
                                                   // min search define the starting range for the
                                                   // peak search function

  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(minsearch), h1->GetXaxis()->FindFixBin(max));
  double binwidth = h1->GetBinWidth(1);


  // Peak search function
  TSpectrum *spec = new TSpectrum();
  spec->Search(h1, search, "", s); // the parameters of the search should be carefully examined.
  int npeaks=spec->GetNPeaks();
  //cout << "npeaks =  " << npeaks << endl;
  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(min), h1->GetXaxis()->FindFixBin(max));
  cFit->Update();
  double emin = start, emax = end;// Range of channels to be fitted

  // definition of the fit function and of an individual Skew Gaussian for each peak. The individual
  // Skew Gaussian is only used for plotting
  vector <TF1*> SkewGaussian;
  TF1 *FitFunction = new TF1("FitFunction", fitfunctionskewgaussian, emin, emax, 7+npeaks*2);
  
  FitFunction->FixParameter(0, npeaks); // parameter 0 is the number of peaks. It is given by the
                                        // search function and must be fixed

  // Background initialisation. Use set parameter to enable background. Use fix parameter to have no
  // background
  FitFunction->FixParameter(1, 0);
  FitFunction->FixParameter(2, 0);
  FitFunction->FixParameter(3, 0);
  
  
  FitFunction->SetNpx((emax-emin)*10);

  FitFunction->SetParameter(4, 0.3); // alpha parameter
  FitFunction->SetParLimits(4, 0, 3); // The current boundaries should work
  
  FitFunction->SetParameter(5, 0.3); // p parameter
  FitFunction->SetParLimits(5, 0, 1); // by definition p should be between 0 and 1

  FitFunction->SetParameter(6, 5); // sigma parameter
  //FitFunction->SetParLimits(6, 5, 40); // 
  
  // creating a background function to plot it on the histogram
  TF1 *Background = new TF1("Background", background, emin, emax, 3);
  stringstream ss;

  // creating all the Skew Gaussians
  for(int i=0; i<npeaks; i++){
    
    ss.str("");
    ss << "SkewGaussian" << i;
    SkewGaussian.push_back(new TF1(ss.str().c_str(),skewgaussian, emin, emax, 5));
    SkewGaussian[i]->SetParNames("alpha", "p", "Mean", "Std_Dev", "Counts");
    SkewGaussian[i]->SetNpx(10000);
    
   }


  // peakpos and peakmax will store the peak positions given by the search function.
  // sigmamin and sigmamax are the boundaries of the allowed range for the sigma of the gaussian
  // sigma is the initial sigma value
  // counts is to be computed using the peakmax and the initial sigma value. It will be used for
  // initialisation but also for integral boundaries definition
  // width is half of the range allowed to move the peak position
  double peakpos, peakmax, sigma = 11, counts, width = sigma;
  TFitResultPtr fitRes;

  //cout << "Npeaks : " << npeaks << endl;

  for(int peak=0; peak<npeaks; peak++){
    
    peakpos = spec->GetPositionX()[peak];
    peakmax = spec->GetPositionY()[peak];
    counts = peakmax * sigma * TMath::Sqrt(2*TMath::Pi());
    
    FitFunction->SetParameter(7 + peak * 2, peakpos); // position of the peak
    FitFunction->SetParLimits(7 + peak * 2, peakpos-width, peakpos+width);

    FitFunction->SetParameter(8 + peak * 2, counts); // Integral
    FitFunction->SetParLimits(8 + peak * 2, 0.01*counts, 100*counts);
    
  }

  // fitting the histogram
  fitRes = h1->Fit("FitFunction", "LSMQE", "SAME", emin, emax);


  //cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
  
  // Setting the background function for the plot
  Background->SetParameter(0, fitRes->Parameter(1));
  Background->SetParameter(1, fitRes->Parameter(2));
  Background->SetParameter(2, fitRes->Parameter(3));
  Background->SetLineColor(4);
  

  // drawing background and fit function
  Background->Draw("SAME");
  FitFunction->Draw("SAME");

  //cout << endl << "Bin width : " << binwidth << " keV" << endl;

  //cout << endl << "Fit function integral : " << FitFunction->Integral(emin,emax)/binwidth << endl;
  //cout << "Spectrum integral : " << h1->Integral(emin/binwidth,emax/binwidth) << endl;
  //cout << "Ratio : " << FitFunction->Integral(emin,emax)/(binwidth*h1->Integral(emin/binwidth,emax/binwidth)) << endl;

  double alpha = fitRes->Parameter(4), p = fitRes->Parameter(5), STD = fitRes->Parameter(6);
  double alphaerr = fitRes->ParError(4), perr = fitRes->ParError(5), STDerr = fitRes->ParError(6);

  // plotting each individual skew gaussian and printing the peaks informations
  for(int peak = 0; peak<npeaks; peak++){
    double Mean = fitRes->Parameter(7 + peak * 2), Counts = fitRes->Parameter(8 + peak * 2)/binwidth;
    double Meanerr = fitRes->ParError(7 + peak * 2), Countserr = fitRes->ParError(8 + peak * 2)/binwidth;

    SkewGaussian[peak]->SetParameter(0, alpha);
    SkewGaussian[peak]->SetParameter(1, p);
    SkewGaussian[peak]->SetParameter(2, Mean);
    SkewGaussian[peak]->SetParameter(3, STD);
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
    Peaks->push_back(Peak);
    
  }


  cFit->Update();
  //cFit->WaitPrimitive();

  //cFit->Close();
  //h1->Delete();
  
  return npeaks;
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


int FitSpectraSkewGaussian2(TH1D *h1, TCanvas * cFit, vector<vector<double>> *Peaks, int search = 5, double s = 0.005, double start = 3000, double end = 3600){

  //cout << "////////////////////////////////////////" << endl;
  double min = start, max = end, minsearch = start; // min and max define the plot range.
                                                   // min search define the starting range for the
                                                   // peak search function

  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(minsearch), h1->GetXaxis()->FindFixBin(max));
  double binwidth = h1->GetBinWidth(1);


  // Peak search function
  TSpectrum *spec = new TSpectrum();
  spec->Search(h1, search, "", s); // the parameters of the search should be carefully examined.
  int npeaks=spec->GetNPeaks();
  //cout << "npeaks =  " << npeaks << endl;
  h1->GetXaxis()->SetRange(h1->GetXaxis()->FindFixBin(min), h1->GetXaxis()->FindFixBin(max));
  cFit->Update();
  double emin = start, emax = end;// Range of channels to be fitted

  // definition of the fit function and of an individual Skew Gaussian for each peak. The individual
  // Skew Gaussian is only used for plotting
  vector <TF1*> SkewGaussian;
  TF1 *FitFunction = new TF1("FitFunction", fitfunctionskewgaussian2, emin, emax, 3+npeaks*2);
  
  FitFunction->FixParameter(0, npeaks); // parameter 0 is the number of peaks. It is given by the
                                        // search function and must be fixed
  
  FitFunction->SetNpx((emax-emin)*10);

  FitFunction->SetParameter(1, 5); // tau parameter
  //FitFunction->SetParLimits(1, 2, 8); // The current boundaries should work
  
  FitFunction->SetParameter(2, 5); // sigma parameter
  //FitFunction->SetParLimits(2, 4.5, 8.5); // 
  
  stringstream ss;

  // creating all the Skew Gaussians
  for(int i=0; i<npeaks; i++){
    
    ss.str("");
    ss << "SkewGaussian" << i;
    SkewGaussian.push_back(new TF1(ss.str().c_str(),skewgaussian2, emin, emax, 4));
    SkewGaussian[i]->SetParNames("tau", "Std_Dev", "Mean", "Counts");
    SkewGaussian[i]->SetNpx(10000);
    
   }


  // peakpos and peakmax will store the peak positions given by the search function.
  // sigmamin and sigmamax are the boundaries of the allowed range for the sigma of the gaussian
  // sigma is the initial sigma value
  // counts is to be computed using the peakmax and the initial sigma value. It will be used for
  // initialisation but also for integral boundaries definition
  // width is half of the range allowed to move the peak position
  double peakpos, peakmax, sigma = 5, counts, width = 3*sigma;
  TFitResultPtr fitRes;

  //cout << "Npeaks : " << npeaks << endl;

  for(int peak=0; peak<npeaks; peak++){
    
    peakpos = spec->GetPositionX()[peak];
    peakmax = spec->GetPositionY()[peak];
    counts = peakmax * sigma * TMath::Sqrt(2*TMath::Pi());
    
    FitFunction->SetParameter(3 + peak * 2, peakpos); // position of the peak
    FitFunction->SetParLimits(3 + peak * 2, peakpos-width, peakpos+width);

    FitFunction->SetParameter(4 + peak * 2, counts); // Integral
    //FitFunction->SetParLimits(4 + peak * 2, 0.01*counts, 100*counts);
    
  }

  // fitting the histogram
  fitRes = h1->Fit("FitFunction", "LSQ", "SAME", emin, emax);


  //cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;
  
  
  // drawing fit function
  FitFunction->Draw("SAME");

  /*
  cout << endl << "Bin width : " << binwidth << " keV" << endl;

  cout << endl << "Fit function integral : " << FitFunction->Integral(emin,emax)/binwidth << endl;
  cout << "Spectrum integral : " << h1->Integral(emin/binwidth,emax/binwidth) << endl;
  cout << "Ratio : " << FitFunction->Integral(emin,emax)/(binwidth*h1->Integral(emin/binwidth,emax/binwidth)) << endl;
  */
  
  double tau = fitRes->Parameter(1), STD = fitRes->Parameter(2);
  double tauerr = fitRes->ParError(1), STDerr = fitRes->ParError(2);
  
  // plotting each individual skew gaussian and printing the peaks informations
  for(int peak = 0; peak<npeaks; peak++){
    double Mean = fitRes->Parameter(3 + peak * 2), Counts = fitRes->Parameter(4 + peak * 2)/binwidth;
    double Meanerr = fitRes->ParError(3 + peak * 2), Countserr = fitRes->ParError(4 + peak * 2)/binwidth;

    SkewGaussian[peak]->SetParameter(0, tau);
    SkewGaussian[peak]->SetParameter(1, STD);
    SkewGaussian[peak]->SetParameter(2, Mean);
    SkewGaussian[peak]->SetParameter(3, Counts*binwidth);
    SkewGaussian[peak]->SetLineColor(3);
    SkewGaussian[peak]->Draw("SAME");



    /*
    cout << "////////////////////////////////////////" << endl;
    cout << "Peak number : " << peak << endl;
    cout << "Peak Position : " << Mean << " +/- " << Meanerr << endl << endl;
    cout << "tau : " << tau << " +/- " << tauerr << endl;
    cout << "Mean : " << Mean << " +/- " << Meanerr << endl;
    cout << "Sigma : " << STD << " +/- " << STDerr << endl;
    cout << "Counts : " << Counts << " +/- " << Countserr << endl;
    */
    vector<double> Peak;
    Peak.push_back(Mean);
    Peak.push_back(Meanerr);
    Peak.push_back(STD);
    Peak.push_back(STDerr);
    Peaks->push_back(Peak);
    
  }
  
  cFit->Update();

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



int FitAllGaussian(double min = 3000, double max = 3600){
  string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"};

  TH1D * hout = new TH1D();
  hout->Add((TH1D*)gDirectory->Get("spec0"));
  TH1D * hfoil = new TH1D();
  hfoil->Add((TH1D*)gDirectory->Get("spec4"));
  
  for(int i=1;i<3;i++){
    hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=3;i<8;i++){
    hfoil->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=8;i<11;i++){
    hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  vector<vector<double>> Peaksout, Peaksfoil;
  //cout << "Results out of the foil : " << endl;
  TCanvas *c1 =  new TCanvas("c1", "c1", 800, 600);
  FitSpectraGaussian(hout,c1,&Peaksout,5,0.005,min,max);
  //cout << "Results in the foil : " << endl;
  TCanvas *c2 =  new TCanvas("c2", "c2", 800, 600);
  FitSpectraGaussian(hfoil,c2,&Peaksfoil,5,0.005,min,max);

  sort(Peaksout.begin(),Peaksout.end(),sorting_peaks_energy);
  sort(Peaksfoil.begin(),Peaksfoil.end(),sorting_peaks_energy);

  double Litt[5] = {5388,5442.80,5485.56, 5762.64, 5804.77};
  double Litterr[5] = {1,0.13,0.12, 0.03, 0.05};
  double PeaksoutPos[5];
  double PeaksoutPoserr[5];
  double PeaksfoilPos[5];
  double PeaksfoilPoserr[5];
  for(int i=0;i<5;i++){
    PeaksoutPos[i] = Peaksout[i][0];
    PeaksoutPoserr[i] = Peaksout[i][1];
    PeaksfoilPos[i] = Peaksfoil[i][0];
    PeaksfoilPoserr[i] = Peaksfoil[i][1];
  }

  double residue[5];
  double residue_error[5];

  auto c3 = new TCanvas("c3", "c3", 200, 10, 700, 500);
  auto gr = new TGraphErrors(5, PeaksoutPos, Litt, PeaksoutPoserr, Litterr);

  c3->Divide(1,2);
  c3->cd(1);
  gr->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  TFitResultPtr FitRes = gr->Fit(line, "QES", "SAME E");

  double a = FitRes->Parameter(1), sa = FitRes->ParError(1), b = FitRes->Parameter(0), sb = FitRes->ParError(0);

  for(int i = 0; i<5; i++){
    double C = PeaksoutPos[i], sC = PeaksoutPoserr[i];
    residue[i] = a * C + b - Litt[i];
    cout << "Residue n°" << i+1 << " : " << residue[i] << endl;
    residue_error[i] = TMath::Sqrt(a*a*sC*sC+sa*sa*C*C+sb*sb);
  }

  TF1 *line2 = new TF1("Line", "pol0", 0, 10000);
  line2->SetParameter(0,0);
  line2->SetLineColor(2);
  auto gr2 = new TGraphErrors(5, PeaksoutPos, residue, PeaksoutPoserr, residue_error);
  c3->cd(2);
  gr2->Draw("ap");
  line2->Draw("same");

  cout << " Chi2/NdF : " << FitRes->Chi2() << "/" << FitRes->Ndf() << " = " << FitRes->Chi2()/FitRes->Ndf() << " => Prob : " << FitRes->Prob() << endl << endl;

  cout << "Slope : " << a << " Slope error : " << sa << endl;
  cout << "Const : " << b << " Const error : " << sb << endl;

  double dE[5];
  double sdE[5];
  double dx[5];
  double sdx[5];
  double SE[5] = {0.1730,0.1717,0.1706,0.1649,0.1641};
  double w[5];
  double Sw = 0, mdx = 0, smdx;

  for(int i=0;i<5;i++){
    double Cout = PeaksoutPos[i], sCout = PeaksoutPoserr[i], Cfoil = PeaksfoilPos[i], sCfoil = PeaksfoilPoserr[i];
    PeaksoutPos[i] = a*Cout+b;
    PeaksoutPoserr[i] = TMath::Sqrt(a*a*sCout*sCout+sa*sa*Cout*Cout+sb*sb);
    PeaksfoilPos[i] = a*Cfoil+b;
    PeaksfoilPoserr[i] = TMath::Sqrt(a*a*sCfoil*sCfoil+sa*sa*Cfoil*Cfoil+sb*sb);
    dE[i] = (Cout-Cfoil)*a;
    sdE[i] = TMath::Sqrt(a*a*(sCout*sCout+sCfoil*sCfoil)+(Cout-Cfoil)*(Cout-Cfoil)*sa*sa);
    dx[i] = dE[i]/SE[i];
    sdx[i] = dx[i]*TMath::Sqrt(sdE[i]*sdE[i]/(dE[i]*dE[i])+0.05*0.05);
    w[i] = 1/(sdx[i]*sdx[i]);
    Sw += w[i];
    mdx += dx[i]*w[i];

    cout << "//////////////////////////////////////////////" << endl;
    cout << "Peak position in litterature : " << Litt[i] << endl;
    cout << "Measured peak position without the foil : " << PeaksoutPos[i] << " +/- " << PeaksoutPoserr[i] << endl;
    cout << "Measured peak position with the foil : " << PeaksfoilPos[i] << " +/- " << PeaksfoilPoserr[i] << endl;
    cout << "Delta E : " << dE[i] << " +/- " << sdE[i] << endl;
    cout << "Delta x : " << dx[i] << " +/- " << sdx[i] << endl << endl;
  }

  mdx =  mdx/Sw;
  smdx = 1/TMath::Sqrt(Sw);

  cout << "//////////////////////////////////////////////" << endl;
  cout << "Mean delta x : " << mdx << " +/- " << smdx << endl;


  return 0;
}


int FitAllSkewGaussian(int Nfoil, double min = 3000, double max = 3600){

  TH1D * hout = new TH1D();
  TH1D * hfoil = new TH1D();
  
  if(Nfoil==15){
    string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"}; // Foil 15
    hout->Add((TH1D*)gDirectory->Get("spec0"));
    hfoil->Add((TH1D*)gDirectory->Get("spec4"));
    
    for(int i=1;i<3;i++){
      hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
    
    for(int i=4;i<8;i++){
      hfoil->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
    
    for(int i=8;i<11;i++){
      hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
  }
  else if(Nfoil==11){
    string hnames[12] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec10", "spec11", "spec12"}; // Foil 11
    hout->Add((TH1D*)gDirectory->Get("spec0"));
    hfoil->Add((TH1D*)gDirectory->Get("spec3"));
    
    for(int i=1;i<3;i++){
      hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
    
    for(int i=4;i<9;i++){
      hfoil->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
    }
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  hout->SetTitle("All Data without the foil;Energy (keV);Counts/keV");
  hfoil->SetTitle("All Data with the foil;Energy (keV);Counts/keV");
  
  vector<vector<double>> Peaksout, Peaksfoil;
  cout << "Results out of the foil : " << endl;
  TCanvas *c1 =  new TCanvas("c1", "c1", 800, 600);
  FitSpectraSkewGaussian(hout,c1,&Peaksout,5,0.005,min,max);
  cout << "Results in the foil : " << endl;
  TCanvas *c2 =  new TCanvas("c2", "c2", 800, 600);
  FitSpectraSkewGaussian(hfoil,c2,&Peaksfoil,5,0.005,min,max);

  sort(Peaksout.begin(),Peaksout.end(),sorting_peaks_energy);
  sort(Peaksfoil.begin(),Peaksfoil.end(),sorting_peaks_energy);

  double Litt[5] = {5388,5442.80,5485.56, 5762.64, 5804.77};
  double Litterr[5] = {1,0.13,0.12, 0.03, 0.05};
  double PeaksoutPos[5];
  double PeaksoutPoserr[5];
  double FWHMout = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][2], sFWHMout = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][3], FWHMfoil = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][2], sFWHMfoil = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][3];
  double PeaksfoilPos[5];
  double PeaksfoilPoserr[5];
  
  for(int i=0;i<5;i++){
    PeaksoutPos[i] = Peaksout[i][0];
    PeaksoutPoserr[i] = Peaksout[i][1];
    PeaksfoilPos[i] = Peaksfoil[i][0];
    PeaksfoilPoserr[i] = Peaksfoil[i][1];
  }

  double residue[5];
  double residue_error[5];

  auto c3 = new TCanvas("c3", "c3", 200, 10, 700, 500);
  auto gr = new TGraphErrors(5, PeaksoutPos, Litt, PeaksoutPoserr, Litterr);

  c3->Divide(1,2);
  c3->cd(1);
  gr->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  TFitResultPtr FitRes = gr->Fit(line, "QES", "SAME E");
  double a = FitRes->Parameter(1), sa = FitRes->ParError(1), b = FitRes->Parameter(0), sb = FitRes->ParError(0);

  
  for(int i = 0; i<5; i++){
    double C = PeaksoutPos[i], sC = PeaksoutPoserr[i];
    residue[i] = a * C + b - Litt[i];
    cout << "Residue n°" << i+1 << " : " << residue[i] << endl;
    residue_error[i] = TMath::Sqrt(a*a*sC*sC+sa*sa*C*C+sb*sb);
  }


  TF1 *line2 = new TF1("Line", "pol0", 0, 10000);
  line2->SetParameter(0,0);
  line2->SetLineColor(2);
  auto gr2 = new TGraphErrors(5, PeaksoutPos, residue, PeaksoutPoserr, residue_error);
  c3->cd(2);
  gr2->Draw("ap");
  line2->Draw("same");

  cout << " Chi2/NdF : " << FitRes->Chi2() << "/" << FitRes->Ndf() << " = " << FitRes->Chi2()/FitRes->Ndf() << " => Prob : " << FitRes->Prob() << endl << endl;

  cout << "Slope : " << a << " Slope error : " << sa << endl;
  cout << "Const : " << b << " Const error : " << sb << endl;

  sFWHMout = TMath::Sqrt(a*a*sFWHMout*sFWHMout+sa*sa*FWHMout*FWHMout);
  FWHMout *= a;
  sFWHMfoil = TMath::Sqrt(a*a*sFWHMfoil*sFWHMfoil+sa*sa*FWHMfoil*FWHMfoil);
  FWHMfoil *= a;

  double dE[5];
  double sdE[5];
  double dx[5];
  double sdx[5];
  //double SE[5] = {0.1730,0.1717,0.1706,0.1649,0.1641};
  double SE[5] = {0.1662,0.1651,0.1641,0.1586,0.1578};
  double sSE = 0.05;
  double w[5];
  double wE[5];
  double Sw = 0, mdx = 0, smdx, SwE=0, mdE=0, smdE, stragg, sstragg, K1, K2;
  double meandx = 0, dxstd = 0;

  stragg = TMath::Sqrt(FWHMfoil*FWHMfoil-FWHMout*FWHMout);
  sstragg = TMath::Sqrt(TMath::Power(sFWHMfoil*FWHMfoil/stragg,2)+TMath::Power(sFWHMout*FWHMout/stragg,2));

  K1 = TMath::Sqrt(13.52*13.52-10.9*10.9);
  K2 = TMath::Sqrt(TMath::Power(13.52*0.05/K1,2)+TMath::Power(10.9*0.08/K1,2));

  for(int i=0;i<5;i++){
    double Cout = PeaksoutPos[i], sCout = PeaksoutPoserr[i], Cfoil = PeaksfoilPos[i], sCfoil = PeaksfoilPoserr[i];
    PeaksoutPos[i] = a*Cout+b;
    PeaksoutPoserr[i] = TMath::Sqrt(a*a*sCout*sCout+sa*sa*Cout*Cout+sb*sb);
    PeaksfoilPos[i] = a*Cfoil+b;
    PeaksfoilPoserr[i] = TMath::Sqrt(a*a*sCfoil*sCfoil+sa*sa*Cfoil*Cfoil+sb*sb);
    dE[i] = (Cout-Cfoil)*a;
    sdE[i] = TMath::Sqrt(a*a*(sCout*sCout+sCfoil*sCfoil)+(Cout-Cfoil)*(Cout-Cfoil)*sa*sa);
    dx[i] = dE[i]/SE[i];
    sdx[i] = dx[i]*TMath::Sqrt(sdE[i]*sdE[i]/(dE[i]*dE[i])+sSE*sSE);
    w[i] = 1/(sdx[i]*sdx[i]);
    wE[i] = 1/(sdE[i]*sdE[i]);
    Sw += w[i];
    SwE += wE[i];
    mdx += dx[i]*w[i];
    meandx += dx[i]/5;
    mdE += dE[i]*wE[i];

    cout << "//////////////////////////////////////////////" << endl;
    cout << "Peak position in litterature : " << Litt[i] << " keV" << endl;
    cout << "Measured peak position without the foil : " << PeaksoutPos[i] << " +/- " << PeaksoutPoserr[i] << " keV" << endl;
    cout << "Measured peak position with the foil : " << PeaksfoilPos[i] << " +/- " << PeaksfoilPoserr[i] << " keV" << endl;
    cout << "Delta E : " << dE[i] << " +/- " << sdE[i] << " keV" << endl;
    cout << "Delta x : " << dx[i] << " +/- " << sdx[i] << " nm" << endl << endl;
  }

  for(int i=0;i<5;i++) dxstd += (dx[i]-meandx)*(dx[i]-meandx)/5;

  dxstd = TMath::Sqrt(dxstd);

  mdx =  mdx/Sw;
  smdx = 1/TMath::Sqrt(Sw);
  mdE = mdE/SwE;
  smdE = 1/TMath::Sqrt(SwE);

  cout << "//////////////////////////////////////////////" << endl;
  cout << "FWHM without the foil : " << FWHMout << " +/- " << sFWHMout << " keV" << endl;
  cout << "FWHM with the foil : " << FWHMfoil << " +/- " << sFWHMfoil << " keV" << endl << endl;;

  cout << "//////////////////////////////////////////////" << endl;
  cout << "Weighted average mean delta E : " << mdE << " +/- " << smdE << " keV" << endl;
  cout << "Weighted average mean delta x : " << mdx << " +/- " << smdx << " nm" << endl;
  cout << "Straggling : " << stragg << " +/- " << sstragg << " keV" << endl;
  cout << "Simulated straggling : " << K1 << " +/- " << K2 << " keV" << endl;


  /*
  for(int i=0;i<5;i++){
    cout << sdE[i]*sdE[i]/(dE[i]*dE[i])/(sSE*sSE) << " " << sdE[i]*sdE[i]/(dE[i]*dE[i])/(0.01*0.01) << endl;
  }
  */

  
  return 0;
}




int FitAllSkewGaussian2(double min = 3000, double max = 3600){
  string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"};
  //TCanvas *cFit =  new TCanvas("cFit", "cFit", 800, 600);

  //TH1D * h = new TH1D("h","h",3000,1000,4000);
  TH1D * hout = new TH1D();
  hout->Add((TH1D*)gDirectory->Get("spec0"));
  TH1D * hfoil = new TH1D();
  hfoil->Add((TH1D*)gDirectory->Get("spec4"));

  //hprev->SetTitle("");
  
  for(int i=1;i<3;i++){
    hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=3;i<8;i++){
    hfoil->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=8;i<11;i++){
    hout->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  vector<vector<double>> Peaksout, Peaksfoil;
  cout << "Results out of the foil : " << endl;
  TCanvas *c1 =  new TCanvas("c1", "c1", 800, 600);
  FitSpectraSkewGaussian2(hout,c1,&Peaksout,5,0.005,min,max);
  cout << "Results in the foil : " << endl;
  TCanvas *c2 =  new TCanvas("c2", "c2", 800, 600);
  FitSpectraSkewGaussian2(hfoil,c2,&Peaksfoil,5,0.005,min,max);

  sort(Peaksout.begin(),Peaksout.end(),sorting_peaks_energy);
  sort(Peaksfoil.begin(),Peaksfoil.end(),sorting_peaks_energy);

  double Litt[5] = {5388,5442.80,5485.56, 5762.64, 5804.77};
  double Litterr[5] = {1,0.13,0.12, 0.03, 0.05};
  double PeaksoutPos[5];
  double PeaksoutPoserr[5];
  double PeaksfoilPos[5];
  double PeaksfoilPoserr[5];
  for(int i=0;i<5;i++){
    PeaksoutPos[i] = Peaksout[i][0];
    PeaksoutPoserr[i] = Peaksout[i][1];
    PeaksfoilPos[i] = Peaksfoil[i][0];
    PeaksfoilPoserr[i] = Peaksfoil[i][1];
  }

  double residue[5];
  double residue_error[5];

  auto c3 = new TCanvas("c3", "c3", 200, 10, 700, 500);
  auto gr = new TGraphErrors(5, PeaksoutPos, Litt, PeaksoutPoserr, Litterr);

  c3->Divide(1,2);
  c3->cd(1);
  gr->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  TFitResultPtr FitRes = gr->Fit(line, "QES", "SAME E");
  double a = FitRes->Parameter(1), sa = FitRes->ParError(1), b = FitRes->Parameter(0), sb = FitRes->ParError(0);
  
  for(int i = 0; i<5; i++){
    double C = PeaksoutPos[i], sC = PeaksoutPoserr[i];
    residue[i] = a * C + b - Litt[i];
    cout << "Residue n°" << i+1 << " : " << residue[i] << endl;
    residue_error[i] = TMath::Sqrt(a*a*sC*sC+sa*sa*C*C+sb*sb);
  }

  TF1 *line2 = new TF1("Line", "pol0", 0, 10000);
  line2->SetParameter(0,0);
  line2->SetLineColor(2);
  auto gr2 = new TGraphErrors(5, PeaksoutPos, residue, PeaksoutPoserr, residue_error);
  c3->cd(2);
  gr2->Draw("ap");
  line2->Draw("same");

  cout << " Chi2/NdF : " << FitRes->Chi2() << "/" << FitRes->Ndf() << " = " << FitRes->Chi2()/FitRes->Ndf() << " => Prob : " << FitRes->Prob() << endl << endl;

  cout << "Slope : " << a << " Slope error : " << sa << endl;
  cout << "Const : " << b << " Const error : " << sb << endl;


  double dE[5];
  double sdE[5];
  double dx[5];
  double sdx[5];
  double SE[5] = {0.1730,0.1717,0.1706,0.1649,0.1641};
  double w[5];
  double Sw = 0, mdx = 0, smdx;

  for(int i=0;i<5;i++){
    double Cout = PeaksoutPos[i], sCout = PeaksoutPoserr[i], Cfoil = PeaksfoilPos[i], sCfoil = PeaksfoilPoserr[i];
    PeaksoutPos[i] = a*Cout+b;
    PeaksoutPoserr[i] = TMath::Sqrt(a*a*sCout*sCout+sa*sa*Cout*Cout+sb*sb);
    PeaksfoilPos[i] = a*Cfoil+b;
    PeaksfoilPoserr[i] = TMath::Sqrt(a*a*sCfoil*sCfoil+sa*sa*Cfoil*Cfoil+sb*sb);
    dE[i] = (Cout-Cfoil)*a;
    sdE[i] = TMath::Sqrt(a*a*(sCout*sCout+sCfoil*sCfoil)+(Cout-Cfoil)*(Cout-Cfoil)*sa*sa);
    dx[i] = dE[i]/SE[i];
    sdx[i] = dx[i]*TMath::Sqrt(sdE[i]*sdE[i]/(dE[i]*dE[i])+0.05*0.05);
    w[i] = 1/(sdx[i]*sdx[i]);
    Sw += w[i];
    mdx += dx[i]*w[i];

    cout << "//////////////////////////////////////////////" << endl;
    cout << "Peak position in litterature : " << Litt[i] << endl;
    cout << "Measured peak position without the foil : " << PeaksoutPos[i] << " +/- " << PeaksoutPoserr[i] << endl;
    cout << "Measured peak position with the foil : " << PeaksfoilPos[i] << " +/- " << PeaksfoilPoserr[i] << endl;
    cout << "Delta E : " << dE[i] << " +/- " << sdE[i] << endl;
    cout << "Delta x : " << dx[i] << " +/- " << sdx[i] << endl << endl;
  }

  mdx =  mdx/Sw;
  smdx = 1/TMath::Sqrt(Sw);

  cout << "//////////////////////////////////////////////" << endl;
  cout << "Mean delta x : " << mdx << " +/- " << smdx << endl;
  
  return 0;
}



int CheckAlignementSE(){
  //string hnames[11] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11"}; // Foil 15
  //string hnames[12] = {"spec0", "spec1", "spec2", "spec3", "spec4", "spec5", "spec6", "spec7", "spec8", "spec10", "spec11", "spec12"}; // Foil 11
  string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 7

  TH1D * hstart = new TH1D();
  hstart->Add((TH1D*)gDirectory->Get("spec0"));
  TH1D * hfoil = new TH1D();
  hfoil->Add((TH1D*)gDirectory->Get("spec4"));
  TH1D * hend = new TH1D();
  hend->Add((TH1D*)gDirectory->Get("spec10"));

  //hprev->SetTitle("");
  
  for(int i=1;i<3;i++){
    hstart->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=3;i<9;i++){
    hfoil->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  for(int i=10;i<12;i++){
    hend->Add((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  auto c1 = new TCanvas("c1","c1",700,500);
  hfoil->Draw("hist");
  hstart->SetLineColor(2);
  hend->SetLineColor(3);
  hstart->Draw("hist same");
  hend->Draw("hist same");
  c1->Update();
  return 0;
}

int CheckAlignementF(){
  //string hnames[5] = {"spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 15
  //string hnames[6] = {"spec3", "spec4", "spec5", "spec6", "spec7", "spec8"}; // Foil 11
  string hnames[5] = {"spec5", "spec6", "spec7", "spec8", "spec9"}; // Foil 7

  vector<TH1D*> h;

  //hprev->SetTitle("");
  
  for(int i=0;i<5;i++){
    h.push_back((TH1D*)gDirectory->Get(hnames[i].c_str()));
  }

  auto c1 = new TCanvas("c1","c1",700,500);
  h[0]->Draw("hist");
  for(int i=1;i<5;i++) h[i]->Draw("hist same");
  c1->Update();
  return 0;
}

