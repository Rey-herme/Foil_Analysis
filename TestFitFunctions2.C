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


int FitSpectraSkewGaussian(TH1D *h1, TCanvas * cFit, vector<vector<double>> *Peaks, int search = 5, double s = 0.005, double start = 3000, double end = 3600){

  //cout << "////////////////////////////////////////" << endl;
  double emin = start, emax = end; // min and max define the plot range.
                                                   // min search define the starting range for the
                                                   // peak search function

  h1->GetXaxis()->SetRangeUser(emin, emax);
  double binwidth = h1->GetBinWidth(1);


  // Peak search function
  TSpectrum *spec = new TSpectrum();
  spec->Search(h1, search, "", s); // the parameters of the search should be carefully examined.
  int npeaks=spec->GetNPeaks();
  cFit->Update();

  // definition of the fit function and of an individual Skew Gaussian for each peak. The individual
  // Skew Gaussian is only used for plotting
  vector <TF1*> SkewGaussian;
  TF1 *FitFunction = new TF1("FitFunction", fitfunctionskewgaussian, emin, emax, 4+(npeaks+1)*2);

  // peakpos and peakmax will store the peak positions given by the search function.
  // sigmamin and sigmamax are the boundaries of the allowed range for the sigma of the gaussian
  // sigma is the initial sigma value
  // counts is to be computed using the peakmax and the initial sigma value. It will be used for
  // initialisation but also for integral boundaries definition
  // width is half of the range allowed to move the peak position
  double peakpos, peakmax, sigma = 5, sigmamin = 3, sigmamax = 10, counts, width = sigma;
  
  FitFunction->FixParameter(0, npeaks+1); // parameter 0 is the number of peaks. It is given by the
                                        // search function and must be fixed
  
  FitFunction->SetNpx((emax-emin)*10);
  
  
  FitFunction->SetParameter(1, 0.3); // alpha parameter
  FitFunction->SetParLimits(1, 0, 3); // The current boundaries should work
    
  FitFunction->SetParameter(2, 1); // p parameter
  FitFunction->SetParLimits(2, 0, 1); // by definition p should be between 0 and 1
    
  FitFunction->SetParameter(3, sigma); // sigma parameter
  //FitFunction->SetParLimits(3, sigmamin, sigmamax); // 
  
  stringstream ss;
  // creating all the Skew Gaussians
  for(int i=0; i<npeaks+1; i++){
    ss.str("");
    ss << "SkewGaussian" << i;
    SkewGaussian.push_back(new TF1(ss.str().c_str(),skewgaussian, emin, emax, 5));
    SkewGaussian[i]->SetParNames("alpha", "p", "Std_Dev", "Mean", "Counts");
    SkewGaussian[i]->SetNpx(10000);
   }


  TFitResultPtr fitRes;

  //cout << "Npeaks : " << npeaks << endl;
  double lastcounts= 0;
  double peakposmin = 4000;
  for(int peak=0; peak<npeaks; peak++){
    
    peakpos = spec->GetPositionX()[peak];
    peakmax = spec->GetPositionY()[peak];
    counts = peakmax * sigma * TMath::Sqrt(2*TMath::Pi());

    if(peakpos<peakposmin) {
      peakposmin=peakpos;
      lastcounts = counts;
    }
    
    FitFunction->SetParameter(4 + peak * 2, peakpos); // position of the peak
    FitFunction->SetParLimits(4 + peak * 2, peakpos-width, peakpos+width);

    FitFunction->SetParameter(5 + peak * 2, counts); // Integral
    FitFunction->SetParLimits(5 + peak * 2, 0.01*counts, 100*counts);
    
  }

  double lastpeakpos=peakposmin+25;
  
  FitFunction->SetParameter(4 + npeaks * 2, lastpeakpos); // position of the peak
  FitFunction->SetParLimits(4 + npeaks * 2, lastpeakpos-2*width, lastpeakpos+2*width);
  
  FitFunction->SetParameter(5 + npeaks * 2, lastcounts); // Integral
  FitFunction->SetParLimits(5 + npeaks * 2, 0.01*lastcounts, 100*lastcounts);
  // fitting the histogram
  fitRes = h1->Fit("FitFunction", "LSMQE", "SAME", emin, emax);


  cout << " Chi2/NdF : " << fitRes->Chi2() << "/" << fitRes->Ndf() << " = " << fitRes->Chi2()/fitRes->Ndf() << " => Prob : " << fitRes->Prob() << endl;


  // drawing background and fit function
  FitFunction->Draw("SAME");

  //cout << endl << "Bin width : " << binwidth << " keV" << endl;

  //cout << endl << "Fit function integral : " << FitFunction->Integral(emin,emax)/binwidth << endl;
  //cout << "Spectrum integral : " << h1->Integral(emin/binwidth,emax/binwidth) << endl;
  //cout << "Ratio : " << FitFunction->Integral(emin,emax)/(binwidth*h1->Integral(emin/binwidth,emax/binwidth)) << endl;

  double alpha = fitRes->Parameter(1), p = fitRes->Parameter(2), STD = fitRes->Parameter(3);
  double alphaerr = fitRes->ParError(1), perr = fitRes->ParError(2), STDerr = fitRes->ParError(3);


   cout << "////////////////////////////////////////" << endl;
   cout << "alpha : " << alpha << " +/- " << alphaerr << endl;
   cout << "p : " << p << " +/- " << perr << endl;
   cout << "Sigma : " << STD << " +/- " << STDerr << endl;

  // plotting each individual skew gaussian and printing the peaks informations
  for(int peak = 0; peak<npeaks+1; peak++){
    double Mean = fitRes->Parameter(4 + peak * 2), Counts = fitRes->Parameter(5 + peak * 2)/binwidth;
    double Meanerr = fitRes->ParError(4 + peak * 2), Countserr = fitRes->ParError(5 + peak * 2)/binwidth;

    SkewGaussian[peak]->SetParameter(0, alpha);
    SkewGaussian[peak]->SetParameter(1, p);
    SkewGaussian[peak]->SetParameter(2, STD);
    SkewGaussian[peak]->SetParameter(3, Mean);
    SkewGaussian[peak]->SetParameter(4, Counts*binwidth);
    SkewGaussian[peak]->SetLineColor(3);
    SkewGaussian[peak]->Draw("SAME");

    
    cout << "////////////////////////////////////////" << endl;
    cout << "Peak number : " << peak << endl;
    cout << "Peak Position : " << Mean << " +/- " << Meanerr << endl;
    cout << "Counts : " << Counts << " +/- " << Countserr << endl;
    
    
    vector<double> Peak;
    Peak.push_back(Mean);
    Peak.push_back(Meanerr);
    Peak.push_back(STD);
    Peak.push_back(STDerr);
    Peaks->push_back(Peak);
    
  }

  cout << endl;

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



int FitAllSkewGaussian(int Nfoil, double min = 3000, double max = 3600, bool V = false){

  TH1D * hout = new TH1D();
  TH1D * hfoil = new TH1D();
  
  if(Nfoil==15){
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
  }
  else if(Nfoil==1){
    TFile *f = new TFile("Scan_FE_all_1.root");
    string hnames[12] = {"spec0", "spec1", "spec2", "spec4", "spec5", "spec6", "spec7", "spec8", "spec9", "spec10", "spec11", "spec12"}; // Foil 4
    
    for(int i=0;i<3;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=3;i<9;i++){
      hfoil->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
    
    for(int i=9;i<12;i++){
      hout->Add((TH1D*)f->Get(hnames[i].c_str()));
    }
  }
  else{
    cout << "Foil number not valid !" << endl;
    return 0;
  }
  
  hout->SetTitle("All Data without the foil;Channel (ADC);Counts/keV");
  hfoil->SetTitle("All Data with the foil;Channel (ADC);Counts/keV");
  
  vector<vector<double>> Peaksout, Peaksfoil;
  if(V) cout << "Results out of the foil : " << endl;
  TCanvas *c1 =  new TCanvas("c1", "c1", 800, 600);
  FitSpectraSkewGaussian(hout,c1,&Peaksout,5,0.003,min,max);
  if(V) cout << "Results in the foil : " << endl;
  TCanvas *c2 =  new TCanvas("c2", "c2", 800, 600);
  FitSpectraSkewGaussian(hfoil,c2,&Peaksfoil,5,0.003,min,max);

  sort(Peaksout.begin(),Peaksout.end(),sorting_peaks_energy);
  sort(Peaksfoil.begin(),Peaksfoil.end(),sorting_peaks_energy);
  
  double Litt[] = {5105.81,5143.82,5156.59,5388.25,5442.86,5485.56,5762.65,5804.77};
  double Litterr[] = {0.21,0.21,0.14,0.13,0.12,0.12,0.03,0.05};
  double PeaksoutPos[8];
  double PeaksoutPoserr[8];
  double FWHMout = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][2], sFWHMout = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksout[0][3], FWHMfoil = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][2], sFWHMfoil = 2*TMath::Sqrt(2*TMath::Log(2))*Peaksfoil[0][3];
  double PeaksfoilPos[8];
  double PeaksfoilPoserr[8];

  for(int i=0;i<8;i++){
    PeaksoutPos[i] = Peaksout[i][0];
    PeaksoutPoserr[i] = Peaksout[i][1];
    PeaksfoilPos[i] = Peaksfoil[i][0];
    PeaksfoilPoserr[i] = Peaksfoil[i][1];
  }


  double residue[8];
  double residue_error[8];

  auto c3 = new TCanvas("c3", "c3", 200, 10, 700, 500);
  auto gr = new TGraphErrors(8, PeaksoutPos, Litt, PeaksoutPoserr, Litterr);

  c3->Divide(1,2);
  c3->cd(1);
  gr->Draw("ap");
  TF1 *line = new TF1("Linear_regression", "pol1", 0, 10000);
  TFitResultPtr FitRes = gr->Fit(line, "QES", "SAME E");
  TMatrixDSym cov = FitRes->GetCovarianceMatrix();
  double a = FitRes->Parameter(1), sa = FitRes->ParError(1), b = FitRes->Parameter(0), sb = FitRes->ParError(0), sab = cov[0][1];

  cout << endl << sab << endl;
  for(int i = 0; i<8; i++){
    double C = PeaksoutPos[i], sC = PeaksoutPoserr[i];
    residue[i] = a * C + b - Litt[i];
    if(V) cout << "Residue n°" << i+1 << " : " << residue[i] << endl;
    residue_error[i] = TMath::Sqrt(a*a*sC*sC+sa*sa*C*C+sb*sb+2*C*sab);
  }


  TF1 *line2 = new TF1("Line", "pol0", 0, 10000);
  line2->SetParameter(0,0);
  line2->SetLineColor(2);
  auto gr2 = new TGraphErrors(8, PeaksoutPos, residue, PeaksoutPoserr, residue_error);
  c3->cd(2);
  gr2->Draw("ap");
  line2->Draw("same");

  if(V){
    cout << " Chi2/NdF : " << FitRes->Chi2() << "/" << FitRes->Ndf() << " = " << FitRes->Chi2()/FitRes->Ndf() << " => Prob : " << FitRes->Prob() << endl << endl;

    cout << "Slope : " << a << " Slope error : " << sa << endl;
    cout << "Const : " << b << " Const error : " << sb << endl;
  }
  
  sFWHMout = TMath::Sqrt(a*a*sFWHMout*sFWHMout+sa*sa*FWHMout*FWHMout);
  FWHMout *= a;
  sFWHMfoil = TMath::Sqrt(a*a*sFWHMfoil*sFWHMfoil+sa*sa*FWHMfoil*FWHMfoil);
  FWHMfoil *= a;

  double dE[8];
  double sdE[8];
  double dx[8];
  double sdx[8];
  double SE[8] = {0.180,0.179,0.178,0.173,0.172,0.171,0.1650,0.1641};
  double sSE = 0.05;
  double w[8];
  double wE[8];
  double Sw = 0, mdx = 0, smdx, SwE=0, mdE=0, smdE, stragg, sstragg, K1, K2;
  double meandx = 0, dxstd = 0;

  stragg = TMath::Sqrt(FWHMfoil*FWHMfoil-FWHMout*FWHMout);
  sstragg = TMath::Sqrt(TMath::Power(sFWHMfoil*FWHMfoil/stragg,2)+TMath::Power(sFWHMout*FWHMout/stragg,2));

  K1 = TMath::Sqrt(13.52*13.52-10.9*10.9);
  K2 = TMath::Sqrt(TMath::Power(13.52*0.05/K1,2)+TMath::Power(10.9*0.08/K1,2));

  for(int i=0;i<8;i++){
    double Cout = PeaksoutPos[i], sCout = PeaksoutPoserr[i], Cfoil = PeaksfoilPos[i], sCfoil = PeaksfoilPoserr[i];
    PeaksoutPos[i] = a*Cout+b;
    PeaksoutPoserr[i] = TMath::Sqrt(a*a*sCout*sCout+sa*sa*Cout*Cout+sb*sb+2*Cout*sab);
    PeaksfoilPos[i] = a*Cfoil+b;
    PeaksfoilPoserr[i] = TMath::Sqrt(a*a*sCfoil*sCfoil+sa*sa*Cfoil*Cfoil+sb*sb+2*Cfoil*sab);
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

    if(V){
      cout << "//////////////////////////////////////////////" << endl;
      cout << "Peak position in litterature : " << Litt[i] << " keV" << endl;
      cout << "Measured peak position without the foil : " << PeaksoutPos[i] << " +/- " << PeaksoutPoserr[i] << " keV" << endl;
      cout << "Measured peak position with the foil : " << PeaksfoilPos[i] << " +/- " << PeaksfoilPoserr[i] << " keV" << endl;
      cout << "Delta E : " << dE[i] << " +/- " << sdE[i] << " keV" << endl;
      cout << "Delta x : " << dx[i] << " +/- " << sdx[i] << " nm" << endl << endl;
    }
  }

  for(int i=0;i<8;i++) dxstd += (dx[i]-meandx)*(dx[i]-meandx)/5;

  dxstd = TMath::Sqrt(dxstd);

  mdx =  mdx/Sw;
  smdx = 1/TMath::Sqrt(Sw);
  mdE = mdE/SwE;
  smdE = 1/TMath::Sqrt(SwE);

  if(V){
    cout << "//////////////////////////////////////////////" << endl;
    cout << "FWHM without the foil : " << FWHMout << " +/- " << sFWHMout << " keV" << endl;
    cout << "FWHM with the foil : " << FWHMfoil << " +/- " << sFWHMfoil << " keV" << endl << endl;;
  }
    
  cout << "//////////////////////////////////////////////" << endl;
  cout << "Weighted average mean delta E : " << mdE << " +/- " << smdE << " keV" << endl;
  cout << "Weighted average mean delta x : " << mdx << " +/- " << smdx << " nm" << endl;
  cout << "Straggling : " << stragg << " +/- " << sstragg << " keV" << endl;
  cout << "Simulated straggling : " << K1 << " +/- " << K2 << " keV" << endl;

  
  
  return 0;
}


int CheckAlignementSE(int Nfoil){

  TH1D * hstart = new TH1D();
  TH1D * hfoil = new TH1D();
  TH1D * hend = new TH1D();

  if(Nfoil==15){
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
  if(Nfoil==15){
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

