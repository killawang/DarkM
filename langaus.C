
//-----------------------------------------------------------------------
//
// Convoluted Landau and Gaussian Fitting Function for a certain interval
// for selected data
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//   Markus Friedl (Markus.Friedl@cern.ch)
//  Modified by Yuhan Wang to be used on the data analysis for UM DM group
//  to execute this code, do:
//  root > .x langaus.C
// or
//  root > .x langaus.C++
//
//-----------------------------------------------------------------------

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
using namespace std;



Double_t length; // the length of the fitting interval
Double_t fitlowbd; // select the lower bound for the fit
Double_t fituppbd; // select the upper bound for the fit
Double_t peaklocation;//read the fit parameter

// enter the length of the fitting interval!!
// Enter the interval of your fit
// doing this by choosing the lowerboound and upperbound by entering the fitlowbd and fituppbd






Double_t langaufun(Double_t *x, Double_t *par) {
    
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.
    
    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;
    
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0];
    
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    
    
    ////////not sure about this
    xlow = fitlowbd;
    xupp=fituppbd;
    
    step = (xupp-xlow) / np;
    
    peaklocation = xx;
    
    
    
    
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau((length
                             -xx),mpc,par[0]) / par[0]; // negative landau function
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau((length-xx),mpc,par[0]) / par[0]; // negative landau function
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    
    return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
    // Once again, here are the Landau * Gaussian parameters:
    //   par[0]=Width (scale) parameter of Landau density
    //   par[1]=Most Probable (MP, location) parameter of Landau density
    //   par[2]=Total area (integral -inf to inf, normalization constant)
    //   par[3]=Width (sigma) of convoluted Gaussian function
    //
    // Variables for langaufit call:
    //   his             histogram to fit
    //   fitrange[2]     lo and hi boundaries of fit range
    //   startvalues[4]  reasonable start values for the fit
    //   parlimitslo[4]  lower parameter limits
    //   parlimitshi[4]  upper parameter limits
    //   fitparams[4]    returns the final fit parameters
    //   fiterrors[4]    returns the final fit errors
    //   ChiSqr          returns the chi square
    //   NDF             returns ndf
    
    Int_t i;
    Char_t FunName[100];
    
    sprintf(FunName,"Fitfcn_%s",his->GetName());
    
    TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold) delete ffitold;
    
    TF1 *ffit = new TF1(FunName,langaufun,fitlowbd ,fituppbd,4);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Landau width","MP","Area","Width langgau fit","peak");
    
    for (i=0; i<4; i++) {
        ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
    }
    
    his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
    
    ffit->GetParameters(fitparams);    // obtain fit parameters
    for (i=0; i<4; i++) {
        fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
    }
    ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
    NDF[0] = ffit->GetNDF();           // obtain ndf
    
    return (ffit);              // return fit function
    
}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
    
    // Seaches for the location (x value) at the maximum of the
    // Landau-Gaussian convolute and its full width at half-maximum.
    //
    // The search is probably not very efficient, but it's a first try.
    
    Double_t p,x,fy,fxr,fxl;
    Double_t step;
    Double_t l,lold;
    Int_t i = 0;
    Int_t MAXCALLS = 10000;
    
    
    // Search for maximum
    
    p = params[1] - 0.1 * params[0];
    step = 0.05 * params[0];
    lold = -2.0;
    l    = -1.0;
    
    
    while ( (l != lold) && (i < MAXCALLS) ) {
        i++;
        
        lold = l;
        x = p + step;
        l = langaufun(&x,params);
        
        if (l < lold)
            step = -step/10;
        
        p += step;
    }
    
    if (i == MAXCALLS)
        return (-1);
    
    maxx = x;
    
    fy = l/2;
    
    
    // Search for right x location of fy
    
    p = maxx + params[0];
    step = params[0];
    lold = -2.0;
    l    = -1e300;
    i    = 0;
    
    
    while ( (l != lold) && (i < MAXCALLS) ) {
        i++;
        
        lold = l;
        x = p + step;
        l = TMath::Abs(langaufun(&x,params) - fy);
        
        if (l > lold)
            step = -step/10;
        
        p += step;
    }
    
    if (i == MAXCALLS)
        return (-2);
    
    fxr = x;
    
    
    // Search for left x location of fy
    
    p = maxx - 0.5 * params[0];
    step = -params[0];
    lold = -2.0;
    l    = -1e300;
    i    = 0;
    
    while ( (l != lold) && (i < MAXCALLS) ) {
        i++;
        
        lold = l;
        x = p + step;
        l = TMath::Abs(langaufun(&x,params) - fy);
        
        if (l > lold)
            step = -step/10;
        
        p += step;
    }
    
    if (i == MAXCALLS)
        return (-3);
    
    
    fxl = x;
    
    FWHM = fxr - fxl;
    return (0);
    
}



void langaus() {
    // Fill Histogram
    //
    TCanvas *se = new TCanvas;
    ///////setting the file name
    string name;
    string filetitle;
    string filetale;
    string filenum;
    int fnum;
    
    ofstream out_data("peak.txt");
    for(fnum=0;fnum<23;fnum++){
    filetitle="/Users/Yuhan/Desktop/Reflectivity/RM03/RM3_run";
    
    // convert int into string
    ostringstream convert;
    convert << fnum;
    filenum = convert.str();
    filetale=".mca";
    // here is the file name
    name=filetitle + filenum + filetale;
    //
    cout << "fiting " << endl;
    cout << name<<endl;
    //read the file
    int data[1024];
    int i=0;
    string line;
    ifstream myfile; //("/Users/Yuhan/Desktop/Reflectivity/RM02/RM2_run04.mca");
    myfile.open(name.c_str());
    bool bStartData = false;
    int temp=0;
    int max=0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            if (bStartData && line.find("<END>") == string::npos){
                
                data[i] = atoi(line.c_str());
                if (data[i]>max){
                    max=data[i];
                    temp = i;
                }
                //cout << data[i] << endl;
                i=i+1;
            }
            if (line.find("<DATA>") != string::npos){
                bStartData = true;
            }
        }
        myfile.close();
    }

    
    // setting the interval
    fitlowbd = temp-3;
    fituppbd = temp+6;
    length = fituppbd;
    //begin the fitting procese
    
    TH1F *hSNR = new TH1F("reflectivity","reflectivity",600,0,600);
    
    for (Int_t i=0; i<1024; i++) hSNR->Fill(i,data[i]);
    
    // Fitting SNR histo
    printf("Fitting...\n");
    
    // Setting fit range and start values
    Double_t fr[2];
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    //xlow = fitlowbd;
    //xupp=fituppbd;
    
    pllo[0]=0.5; pllo[1]=5.0; pllo[2]=1.0; pllo[3]=0.4;
    plhi[0]=5.0; plhi[1]=50.0; plhi[2]=1000000.0; plhi[3]=5.0;
    sv[0]=1.8; sv[1]=20.0; sv[2]=50000.0; sv[3]=3.0;
    
    Double_t chisqr;
    Int_t    ndf;
    TF1 *fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
    
    Double_t SNRPeak, SNRFWHM;
    langaupro(fp,SNRPeak,SNRFWHM);
    
    printf("Fitting done\nPlotting results...\n");
    
    // Global style settings
    gStyle->SetOptStat();
    gStyle->SetOptFit(1);
    
   // gStyle->SetLabelSize(0.03,"x");
    //gStyle->SetLabelSize(0.03,"y");
    
    hSNR->GetXaxis()->SetRange(0,1100);
    hSNR->Draw();
    fitsnr->Draw("lsame");
    
    
    /////add a new linr to box
    se->Update();
    // Retrieve the stat box
    TPaveStats *ps = (TPaveStats*)se->GetPrimitive("stats");
    ps->SetName("mystats");
    TList *list = ps->GetListOfLines();
    // Remove the RMS line
    TText *tconst = ps->GetLineWith("RMS");
    list->Remove(tconst);
    // Add a new line in the stat box.
    // Note that "=" is a control character
    
    TLatex *myt = new TLatex(0, 0, TString::Format("Peak = %g", SNRPeak));
    myt=myt;
    myt ->SetTextFont(42);
    myt ->SetTextSize(0.04);
    myt ->SetTextColor(kRed);
    list->Add(myt);
    // the following line is needed to avoid that the automatic redrawing of stats
    hSNR->SetStats(0);
    se->Modified();
    
    
    
    
    
    
    ////
    
    
    //out put the result to a txt
    
    out_data<< "run "<< fnum <<"Peak location "<< SNRPeak<<endl;
        
    
    }
    
    
    
    
    
    
    cout << endl<< endl;

    cout <<"SNRpeak: "<<SNRPeak << endl;
    
    cout <<"SNRFWHM: "<< SNRFWHM << endl;
    
    cout <<"chisqr: "<< chisqr << endl;
    
    cout << "peak from the Gaussian fit for reference " << peaklocation<<endl;
    
    cout << "biggest num" << max<<"("<<temp<<endl;
    
    
    
}
