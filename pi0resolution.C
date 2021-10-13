#include "../common/binningheader.h"
#include "../common/plottingheader.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>

#define NELEMS(arr) (sizeof(arr)/sizeof(arr[0])) 


double mysine(double* x, double* par) { double 

Amplitude = par[0];
double wavelength = par[1];
double phase = par[2];
return Amplitude*sin(2*TMath::Pi()/wavelength*x[0]+phase); 

}


double GaussExpLinear(double* x, double* par) {

    double p0 = par[0];
    double p1 = par[1];
    double p2 = par[2];
    double p3 = par[3];
    double p4 = par[4];
    double p5 = par[5];
    double p6 = par[6];

    double func;
    if(x[0]<p1){
      func = (p0*(TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))+TMath::Exp((x[0]-p1)/p3)*(1.-TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))))+p4+p5*x[0]);

    }else{
      func = (p0*(TMath::Exp(-0.5*pow(((x[0]-p1)/p2),2))+TMath::Exp(-(x[0]-p1)/p6)*(1.-TMath::Exp(-0.5*pow(((x[0]-p1)/p1),2))))+p4+p5*x[0]);
    }
    

    return func;


}

void pi0resolution(
    TString inputFileName     = "central_allculsters.root",
    TString caloName          = "BECAL",
    //TString inputFileName     = "fwd_allculsters.root",
    //TString caloName          = "FEMC",
    //TString inputFileName     = "bck_allculsters.root",
    //TString caloName          = "EEMC",
    TString clusterizerName   =  "MA",
    TString suffix            = "pdf",
    Bool_t doFitting          = kTRUE
){  

  
  Double_t mass_pi0PDG = 0.1349770;

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  StyleSettingsThesis();
  TString collisionSystem           = "e-p, #sqrt{#it{s}} = 100 GeV";
  TString perfLabel                 = "ECCE G4 simulation";
  
  Double_t textSizeSinglePad        = 0.05;
  Double_t textSizeLabelsPixel      = 35;
  Double_t textSizeLabelsRel        = 58./1300;

  Double_t minEtaMax  = -4;
  Double_t maxEtaMax  = 4;
  Double_t minPhiMax  = -4;
  Double_t maxPhiMax  = 4;
  Double_t maxReso    = 1.25;
  TString detLabel    = "";
  Bool_t isEMCal      = 0;
  Int_t rebinRes      = 1;
  Float_t maxResoPer  = 31;
  Float_t maxFitRes1oE  = -1;
  Float_t minFitRes1oE  = -1;
  Float_t maxFitResE    = -1;
  Float_t minFitResE    = -1;
  Float_t maxFitMeanE    = -1;
  Float_t minFitMeanE    = 10;
  TString caloNameRead = caloName;
  Bool_t isComb         = 0;
  
  Int_t exEbin          = 15;
  Int_t exEtabin        = 2;

  
  // fwd Hcal numbers
  Float_t sqrtETerm[2]  = {35, 50};
  Float_t linTerm[2]    = {7, 10};

  int etabine = 0;
  Double_t etabinfwd[7] = { 1.2,1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
  Double_t etabincentral[5] = { -1.5, -1.2, -0.4, 0.4, 1.2};
  Double_t etabinbck[6] = {-4.0, -3.5, -3.0, -2.5, -2.0, -1.5};

  
  if (caloName.CompareTo("EEMC") == 0){
    detLabel  = "EEMC (PbWO_{4} crystal)";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    isEMCal   = 1;
    maxResoPer  = 12.5;
    maxFitRes1oE  = 1.45;
    minFitRes1oE  = 1/TMath::Sqrt(20.);
    maxFitResE    = 20;
    minFitResE    = 0.2;
    exEtabin      = 2;
    sqrtETerm[0]  = 2;
    sqrtETerm[1]  = 3;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
    etabine       = 5;
  } else if (caloName.CompareTo("FEMC") == 0){
    detLabel  = "FEMC (PHENIX re-use)";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    isEMCal   = 1;
    maxResoPer  = 20.5;
    exEtabin      = 11;
    sqrtETerm[0]  = 7;
    sqrtETerm[1]  = 10;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
    maxFitMeanE   = 20;
    etabine =6;


  } else if (caloName.CompareTo("BECAL") == 0){
    detLabel  = "BECAL (Sci-glass)";
    minEtaMax = -1.8;
    maxEtaMax = 1.3;
    isEMCal   = 1;
    maxResoPer  = 17.5;

    exEtabin      = 7;
    sqrtETerm[0]  = 7;
    sqrtETerm[1]  = 10;
    linTerm[0]    = 1;
    linTerm[1]    = 3;
    etabine =5;
  } else if (caloName.CompareTo("LFHCAL") == 0){
    detLabel  = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    exEtabin      = 11;
  } else if (caloName.CompareTo("LFHCAL-wMat") == 0){
    detLabel  = "LFHCAL, FEMC infront";
    caloNameRead = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    exEtabin      = 11;
  } else if (caloName.CompareTo("LFHCAL-FEMC") == 0){
    detLabel  = "LFHCAL+FEMC";
    caloNameRead = "LFHCAL";
    minEtaMax = 1.1;
    maxEtaMax = 4;
    maxReso   = 1.65;
    rebinRes    = 2;
    maxResoPer  = 72;
    isComb      = 1;
    exEtabin      = 11;
  } else if (caloName.CompareTo("EHCAL") == 0){
    detLabel  = "EHCAL (STAR re-use)";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;

  } else if (caloName.CompareTo("EHCAL-wMat") == 0){
    detLabel  = "EHCAL, EEMC infront";
    caloNameRead = "EHCAL";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("EHCAL-EEMC") == 0){
    detLabel  = "EHCAL+EEMC";
    caloNameRead = "EHCAL";
    minEtaMax = -4;
    maxEtaMax = -1.7;
    maxReso   = 1.25;
    maxResoPer  = 62;
    isComb      = 1;
    exEtabin      = 2;
    sqrtETerm[0]  = 45;
    sqrtETerm[1]  = 50;
    linTerm[0]    = 6;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("CHCAL") == 0){
    detLabel  = "oHCAL, iHCal infront";
    caloNameRead = "HCALOUT";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  } else if (caloName.CompareTo("CHCAL-comb") == 0){
    detLabel  = "oHCAL+iHCal";
    caloNameRead = "HCALOUT";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
    isComb      = 1;
  } else if (caloName.CompareTo("HCALIN") == 0){
    detLabel  = "iHCal";
    caloNameRead = "HCALIN";
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    maxReso   = 1.25;
    maxResoPer  = 82;
    isComb      = 0;
    exEtabin      = 2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  } else {
    minEtaMax = -1.8;
    maxEtaMax = 1.2;
    sqrtETerm[0]  = 85;
    sqrtETerm[1]  = 100;
    linTerm[0]    = 7;
    linTerm[1]    = 10;
  }
  // determine eta bins
  Int_t binEtaMin = 0;
  Int_t binEtaMax = 15;
  
  TString outputDir                 = Form("plotsResolution%s",caloName.Data());
  if(doFitting) outputDir+="Fitting";
  gSystem->Exec("mkdir -p "+outputDir);

  TFile* inputFile                      = new TFile(inputFileName.Data());


  //************************** Read data **************************************************
  const int nEne = 31;
  const int nEta = 40;
  const int nPhi = 25;
  
  const static Double_t partE[]   = {   0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5,
                                        10, 12.5, 15, 17.5, 20, 22.5, 25, 30, 35, 40, 50, 67.5
                                    };


  const static Double_t partPhi[]   = {   -3, 2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0,
                                          0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3
                                    };

  const static Double_t partEta[]   = {   -5, -4.75, -4.5, -4.25, -4, -3.75,  -3.5, -3.25, -3,
                                         -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0,
                                          0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 
                                          4.5, 4.75, 5
                                    };


  Int_t numlim = 11;
  const static Double_t limE[]   = {   0.5, 1.,  3., 5., 7.,10., 15.,20, 25, 40, 75 }; 

  Int_t numlimeta = 10;
  const static Double_t limEta[]   = {   -4, -3.5, -3, -2.5, -2, -1.5, -0.75, 0., 0.75, 1.5, 2, 2.5, 3, 3.5, 4  };                         

  double energymax = 40.;



  TGraphErrors* graphReq          = new TGraphErrors(nEne); 
  TGraphErrors* graphReq1oE       = new TGraphErrors(nEne); 


  
  
  Int_t useEta[7][nEta]             = {{ 0}};
  Int_t usePhi[7][nPhi]             = {{ 0}};
  Int_t usePart[7]                  = { 0};
  Int_t usePartPhi[7]               = { 0};
  Int_t usePartEta[7]               = { 0};

  Int_t useParEnetPhi[7]            = { 0};
  Int_t useEnePhi[7][nEne]         = {{ 0}};

  Int_t useParEnetEta[7]            = { 0};
  Int_t useEneEta[7][nEne]         = {{ 0}};
          


  TH2F* histEMinv;
  TH2F* histpTMinv;
  TH2F* histEreconMinv;

  TH1F* histEMinvbin[50]            = {NULL};
  TH2F* histM20p1;
  TH1F* histM20_photon1[50]                 = {NULL};
  TH2F* histM20p2;
  TH1F* histM20_photon2[50]                 = {NULL};
  TH2F* histM02p1;
  TH1F* histM02_photon1[50]                 = {NULL};
  TH2F* histM02p2;
  TH1F* histM02_photon2[50]                 = {NULL};

  TH2F* histM02SC;
  TH1F* histM02_SC[7]                 = {NULL};
  TH2F* histM20SC;
  TH1F* histM20_SC[7]                 = {NULL};

  TF1*  fithistEMinv[50]; 
  TF1*  fithistEMinv2[50];               
  TF1*  fithistEMinvFINAL[50];               

  TH1F* histResoEVsMinv = new TH1F("histResoEVsMinv","",100,0,41);
  TH1F* histMeanEVsMinv = new TH1F("histMeanEVsMinv","",100,0,41);
  TH1F* histSigmaEVsMinv = new TH1F("histSigmaEVsMinv","",100,0,41);
  
  TH1F* histEMinvbinSlice[50]            = {NULL};
  TF1*  fithistEMinvSlice[50];
  TF1*  fithistEMinvSliceFINAL[50];
  
  TH3F* histEtaEMinv;
  TH1F* histEEtaMinvbin[50][50]         = {{NULL}};
  TF1* fitEEtaMinvbin[50][50]           = {{NULL}};
  TF1* fitEEtaMinvbinFINAL[50][50]           = {{NULL}};
  TH1F* histResoEEtaVsMinv[50];
  TH1F* histMeanEEtaVsMinv[50];
  TH1F* histSigmaEEtaVsMinv[50];

  TH2F* histEClusters;
  TH1F* histEClusterbin[50];
  TH1F* histEProb[50];
  TH2F* histM02p1Pt; 
  TH2F* histM02p2Pt; 
  TH2F* histM02SCPt;
  TH2F* histM20p1Pt; 
  TH2F* histM20p2Pt;
  TH2F* histM20SCPt;


  Int_t nParticles                  = 7;
  TString readNames[7]              = {"gamma_", "electron_", "pion_", "proton_", "kaon_", "neutron_", ""};
  TString labelPart[7]              = {"#gamma", "e^{#pm}", "#pi^{#pm}", "p/#bar{p}", "K^{#pm}", "n", "all"};
  Color_t colorPart[7]              = {807, kBlue+1, kRed+2, kGreen+2, kViolet+2, kGray+2, kBlack};
  Style_t markerStylePart[7]        = {20, 46, 33, 30, 28, 42, 25};
  Style_t lineStylePart[7]          = {3, 4, 5, 7, 9, 3, 1};
  Size_t markerSizePart[7]          = {2.5*1.5, 2*1.5, 3*1.5, 2.7*1.5, 2*1.5, 3*1.5, 2*1.5};

  TH1F* hist1DResoEtabins3D[100][100] = {{NULL}};
  TH3F* hist_E_eta_Clusters;

  Int_t nbins = 100; 
  double ebin = 1.0;
  int padebin = 0; 
  double minmassfit = 0.05;
  double maxmassfit = 0.15;
  double mean, meanErr, sigma, sigmaErr;
  double sigmat = 2.3548;
  int clusmax = 10;

  double Ptbin = 1;
  double Ptmax = 30;
  double binpT = 0.5;
  double ptskip = 1.;
  TH1F* histPtMinvbin[50];
  TH3F* hist_Pt_eta_Clusters; 
  TH1F* histPtProb[50];
  TH1F* hist1DResoPttabins3D[50][50] = {{NULL}};

  
  int padptbin = 0;
  TF1* fithistPtMinv[50];
  TF1* fithistPtMinvFINAL[50];
  TH1F* histResoPtVsMinv = new TH1F("histResoEVsMinv","",100,0,31);
  TH1F* histMeanPtVsMinv = new TH1F("histResoEVsMinv","",100,0,31);
  TH1F* histSigmaPtVsMinv = new TH1F("histResoEVsMinv","",100,0,31);

  TString NameE_Minv = Form("%s/h_E_Minv_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histEMinv    = (TH2F*)inputFile->Get(NameE_Minv.Data());
  TString NameErecon_Minv = Form("%s/h_Erecon_Minv_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histEreconMinv    = (TH2F*)inputFile->Get(NameErecon_Minv.Data());

  vector<double> allmean;
  vector<double> allsigma;

  if (histEMinv->GetEntries() > 100){
    TCanvas* Minvpi02D = new TCanvas("Minvpi02D","",0,0,1500,800);
    DrawGammaCanvasSettings( Minvpi02D, 0.095, 0.01, 0.01, 0.105);
    Minvpi02D->Divide(2,1);
    Minvpi02D->cd(1);
    histEMinv->GetXaxis()->SetRangeUser(0.,40);
    histEMinv->GetYaxis()->SetRangeUser(0.,0.2);
    SetStyleHistoTH2ForGraphs(histEMinv,"E^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histEMinv->Draw("colz");
    Minvpi02D->cd(2);
    histEreconMinv->GetXaxis()->SetRangeUser(0.,40);
    histEreconMinv->GetYaxis()->SetRangeUser(0.,0.2);
    SetStyleHistoTH2ForGraphs(histEreconMinv,"E^{reco}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histEreconMinv->Draw("colz");
     
    TLatex Tl1;
    Tl1.SetTextSize(0.035);
    Tl1.DrawLatex(0.55*energymax,0.04,perfLabel);
    Tl1.DrawLatex(0.55*energymax,0.03,Form("Calorimeter: %s ", caloName.Data()));
    Tl1.DrawLatex(0.55*energymax,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
    Tl1.DrawLatex(0.55*energymax,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
    Tl1.Draw();
    
    Minvpi02D->Print(Form("%s/Mass_vs_E_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    while(ebin < energymax){
      histEMinvbin[padebin]    = (TH1F*)histEMinv->ProjectionY(Form("projection_Minv_%f", ebin), 
                                             histEMinv->GetXaxis()->FindBin(ebin-1.), histEMinv->GetXaxis()->FindBin(ebin+1.),"e");

      if(histEMinvbin[padebin]->GetEntries() > 10){

        mean = meanErr = sigma = sigmaErr = 0;

        if (doFitting){   

          fithistEMinv[padebin]    = new TF1(Form("fit_Minv_%f", ebin), "crystalball", 0.05, 0.15);
          fithistEMinv[padebin]->SetParameters(0.55*histEMinvbin[padebin]->GetMaximum(),histEMinvbin[padebin]->GetMean(),2.*histEMinvbin[padebin]->GetRMS(),2.*histEMinvbin[padebin]->GetRMS(),2.5*histEMinvbin[padebin]->GetRMS());
          histEMinvbin[padebin]->Fit(fithistEMinv[padebin],"L0RMEQ","",minmassfit, maxmassfit);
          
          fithistEMinvFINAL[padebin] = histEMinvbin[padebin]->GetFunction(Form("fit_Minv_%f", ebin));
               
          mean     = fithistEMinvFINAL[padebin]->GetParameter(1);
          meanErr  = fithistEMinvFINAL[padebin]->GetParError(1);
          sigma    = fithistEMinvFINAL[padebin]->GetParameter(2);
          sigmaErr = fithistEMinvFINAL[padebin]->GetParError(2);
          allmean.push_back(mean);
          allsigma.push_back(sigma);
        }


        Double_t reso     = sigma/mean;
        Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );
          
        histResoEVsMinv->SetBinContent(histResoEVsMinv->FindBin(ebin),  100*reso);
        histResoEVsMinv->SetBinError(histResoEVsMinv->FindBin(ebin),  100*resoErr);  
        histMeanEVsMinv->SetBinContent(histMeanEVsMinv->FindBin(ebin), mean);
        histMeanEVsMinv->SetBinError(histMeanEVsMinv->FindBin(ebin), meanErr);

        histSigmaEVsMinv->SetBinContent(histSigmaEVsMinv->FindBin(ebin), sigma/sigmat);
        histSigmaEVsMinv->SetBinError(histSigmaEVsMinv->FindBin(ebin), sigmaErr/sigmat);
        padebin++;
      }

      ebin = ebin + 2.; 
    }

    TCanvas* mass2D = new TCanvas("mass2D","",0,0,1000,800);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFontSize(0.07); 
    DrawGammaCanvasSettings( mass2D, 0.095, 0.01, 0.01, 0.10);
    mass2D->Divide(4,5);
    int padnumber = 0;
    TLegend* legend1 = new TLegend(0., 0.01, 0.4, 0.25);
    for(int z=0; z< padebin +1 ; z++ ){
    if(z==3){
        mass2D->cd(z+1);
        TLatex Tl;
        Tl.SetTextSize(0.1);
        Tl.DrawLatex(0.03,0.9,perfLabel);
        Tl.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloName.Data()));
        Tl.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
        Tl.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

        legend1->Draw();
        
      } else{
        mass2D->cd(z+1);
        SetStyleHistoTH1ForGraphs( histEMinvbin[padnumber],Form("%0.1d GeV < E <  %0.1d GeV",2*z, 2*z+2), "E^{MC}", "M_{#gamma #gamma}", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
        DrawGammaSetMarker(histEMinvbin[padnumber],34,1, kBlack, kBlack);
        histEMinvbin[padnumber]->GetYaxis()->SetNdivisions(5);
        histEMinvbin[padnumber]->GetXaxis()->SetNdivisions(5);
        histEMinvbin[padnumber]->GetXaxis()->SetRangeUser(0,0.3);
        histEMinvbin[padnumber]->Draw("ep");
        fithistEMinvFINAL[padnumber]->SetLineColor(kBlue);
        fithistEMinvFINAL[padnumber]->SetLineWidth(1);
        fithistEMinvFINAL[padnumber]->Draw("same");
        DrawGammaLines(allmean[padnumber], allmean[padnumber], 0, 0.5*fithistEMinvFINAL[padnumber]->GetMaximum(),2, 2);
        DrawGammaLines(allmean[padnumber] + allsigma[padnumber], allmean[padnumber]+ allsigma[padnumber], 0, fithistEMinvFINAL[padnumber]->GetMaximum(), 1, 12);
        DrawGammaLines(allmean[padnumber] - allsigma[padnumber], allmean[padnumber]- allsigma[padnumber], 0, fithistEMinvFINAL[padnumber]->GetMaximum(), 1, 12);
        if(z==0){
          legend1->AddEntry(histEMinvbin[padnumber], "Simulation", "p");
          legend1->AddEntry(fithistEMinvFINAL[padnumber], "Fit", "l");
        }
        padnumber++;
     }
    }

    mass2D->Print(Form("%s/Mass_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    //============= Fit Results======================================
    TCanvas* fitresults = new TCanvas("fitresults","",0,0,1000,500);
    DrawGammaCanvasSettings( fitresults, 0.15, 0.15, 0.15, 0.15);
   
    fitresults->Divide(2,1);
    fitresults->cd(1);

    histMeanEVsMinv->GetYaxis()->SetRangeUser(0.08,0.16);
    SetStyleHistoTH1ForGraphs( histMeanEVsMinv, "", "(M_{#gamma #gamma}", "", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.,1.5);
    DrawGammaSetMarker(histMeanEVsMinv,34,1, kBlack, kBlack);   
    histMeanEVsMinv->GetYaxis()->SetRangeUser(0,1.3*histMeanEVsMinv->GetMaximum()); 
    histMeanEVsMinv->Draw("e,p");
    DrawGammaLines(0, energymax, mass_pi0PDG, mass_pi0PDG,2, 13);

    TLatex Tl3;
    Tl3.SetTextSize(0.035);
    Tl3.DrawLatex(2,0.18,perfLabel);
    Tl3.DrawLatex(2,0.17,Form("Calorimeter: %s ", caloName.Data()));
    Tl3.DrawLatex(2,0.16,Form("Clusterization type: %s ", clusterizerName.Data()));
    Tl3.DrawLatex(2,0.14,"#pi ^{0} #rightarrow #gamma #gamma");
    Tl3.Draw();
    fitresults->cd(2);

    histSigmaEVsMinv->GetYaxis()->SetRangeUser(0.,0.2);
    SetStyleHistoTH1ForGraphs(histSigmaEVsMinv, "", "E^{MC}", "Width (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.3);  
    DrawGammaSetMarker(histSigmaEVsMinv,34,1, kBlack, kBlack);   
    histSigmaEVsMinv->GetYaxis()->SetRangeUser(0,0.015); 
    histSigmaEVsMinv->Draw("e,p");
    
  
    fitresults->Print(Form("%s/Mean_Width_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    //============= End Fit Results======================================

  }
  
  //============= M02 and M20 Results======================================
  TString Name_M02_Cluster1 = Form("%s/h_E_M02_photon1_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02p1    = (TH2F*)inputFile->Get(Name_M02_Cluster1.Data());
  TString Name_M02_Cluster2 = Form("%s/h_E_M02_photon2_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02p2    = (TH2F*)inputFile->Get(Name_M02_Cluster2.Data());
  TString Name_M02_SC = Form("%s/h_E_M02_singlecluster_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02SC    = (TH2F*)inputFile->Get(Name_M02_SC.Data());

  if (histM02p1->GetEntries() > 100){
  
    TCanvas* M02D = new TCanvas("M02D","",0,0,1500,800);
    DrawGammaCanvasSettings( M02D, 0.095, 0.01, 0.01, 0.105);
    M02D->Divide(3,1);
    M02D->cd(1);
    histM02p1->GetXaxis()->SetRangeUser(0.,energymax);
    histM02p1->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02p1,"E^{MC}", "M02(#gamma 1 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02p1->Draw("colz");
    M02D->cd(2);
    histM02p2->GetXaxis()->SetRangeUser(0.,energymax);
    histM02p2->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02p2,"E^{MC}", "M02(#gamma 2 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02p2->Draw("colz");
    TPaveText *t2 = new TPaveText(0.6*energymax,0.7,0.95*energymax,0.9,"NB");
    t2->SetFillColor(kWhite);
    t2->AddText(perfLabel); ((TText*)t2->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2->AddText(Form("Calorimeter: %s ", caloName.Data()));  ((TText*)t2->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2->AddText(Form("Clusterization type: %s ", clusterizerName.Data()));  ((TText*)t2->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2->AddText("#pi ^{0} #rightarrow #gamma #gamma");  ((TText*)t2->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2->Draw();
    M02D->cd(3);
    histM02SC->GetXaxis()->SetRangeUser(0.,energymax);
    histM02SC->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02SC,"E^{MC}", "M02(Single Cluster)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02SC->Draw("colz");
    
    M02D->Print(Form("%s/M02_vs_E_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));
   }

  TString Name_M20_Cluster1 = Form("%s/h_E_M20_photon1_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20p1    = (TH2F*)inputFile->Get(Name_M20_Cluster1.Data());
  TString Name_M20_Cluster2 = Form("%s/h_E_M20_photon2_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20p2    = (TH2F*)inputFile->Get(Name_M20_Cluster2.Data());
  TString Name_M20_SC = Form("%s/h_E_M20_singlecluster_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20SC    = (TH2F*)inputFile->Get(Name_M20_SC.Data());

  if (histM02p1->GetEntries() > 100){
  
    TCanvas* M20D = new TCanvas("M20D","",0,0,1500,800);
    DrawGammaCanvasSettings( M20D, 0.095, 0.01, 0.01, 0.105);
    M20D->Divide(3,1);
    M20D->cd(1);
    histM20p1->GetXaxis()->SetRangeUser(0.,energymax);
    histM20p1->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM20p1,"E^{MC}", "M20(#gamma 1 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM20p1->Draw("colz");
    M20D->cd(2);
    histM20p2->GetXaxis()->SetRangeUser(0.,energymax);
    histM20p2->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM20p2,"E^{MC}", "M20(#gamma 2 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM20p2->Draw("colz");
    TPaveText *t = new TPaveText(0.03*energymax,0.7,0.6*energymax,0.95,"NB");
    t->SetFillColor(kWhite);
    t->AddText(perfLabel); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
    t->AddText(Form("Calorimeter: %s ", caloName.Data()));  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
    t->AddText(Form("Clusterization type: %s ", clusterizerName.Data()));  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
    t->AddText("#pi ^{0} #rightarrow #gamma #gamma");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlack);
    t->Draw();
    M20D->cd(3);
    histM20SC->GetXaxis()->SetRangeUser(0.,energymax);
    histM20SC->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM20SC,"E^{MC}", "M20(Single Cluster)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM20SC->Draw("colz");
    
    
    
    M20D->Print(Form("%s/M20_vs_E_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));
  }

  //============= End M02 and M20 Results======================================

  //============= Clusters Number ======================================

   TString Name_E_eta_Clusters = Form("%s/h_E_eta_Clusters_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
   hist_E_eta_Clusters    = (TH3F*)inputFile->Get(Name_E_eta_Clusters.Data());

   TCanvas* ProbCanv = new TCanvas("ProbCanv","",0,0,800,800);
   DrawGammaCanvasSettings( ProbCanv, 0.095, 0.01, 0.1, 0.105);

   TLegend* legendprob = new TLegend(0.65, 0.15, 0.98, 0.35);
   legendprob->SetHeader(perfLabel);
   legendprob->AddEntry((TObject*)0, Form("Calorimeter: %s ", caloName.Data()),"");  
   legendprob->AddEntry((TObject*)0,Form("Clusterization type: %s ", clusterizerName.Data()),"");
   legendprob->AddEntry((TObject*)0, "#pi ^{0} #rightarrow #gamma #gamma",""); 

   for (int ieta = 0 ; ieta < etabine-1; ieta++){
   
   histEProb[ieta]  = new TH1F(Form("histEProb_%i",ieta),"",100,0,41);
   ebin = 1;
   int iebin = 0;  

   while(ebin < energymax){
    if(caloNameRead == "FEMC"){
    hist1DResoEtabins3D[ieta][iebin]    = (TH1F*)hist_E_eta_Clusters->ProjectionZ(Form("hist_E_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_E_eta_Clusters->GetXaxis()->FindBin(ebin - 1), hist_E_eta_Clusters->GetXaxis()->FindBin(ebin + 1),
                                                                              hist_E_eta_Clusters->GetYaxis()->FindBin(etabinfwd[ieta]),hist_E_eta_Clusters->GetYaxis()->FindBin(etabinfwd[ieta+1]),"e");

    }
    else if(caloNameRead == "BECAL"){
    hist1DResoEtabins3D[ieta][iebin]    = (TH1F*)hist_E_eta_Clusters->ProjectionZ(Form("hist_E_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_E_eta_Clusters->GetXaxis()->FindBin(ebin - 1), hist_E_eta_Clusters->GetXaxis()->FindBin(ebin + 1),
                                                                              hist_E_eta_Clusters->GetYaxis()->FindBin(etabincentral[ieta]),hist_E_eta_Clusters->GetYaxis()->FindBin(etabincentral[ieta+1]),"e");

    }
    else if(caloNameRead == "EEMC"){
    hist1DResoEtabins3D[ieta][iebin]    = (TH1F*)hist_E_eta_Clusters->ProjectionZ(Form("hist_E_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_E_eta_Clusters->GetXaxis()->FindBin(ebin - 1), hist_E_eta_Clusters->GetXaxis()->FindBin(ebin + 1),
                                                                              hist_E_eta_Clusters->GetYaxis()->FindBin(etabinbck[ieta]),hist_E_eta_Clusters->GetYaxis()->FindBin(etabinbck[ieta+1]),"e");

    }
    else{
      cout << " Not an electromagnetic calorimeter " << endl;
      return 0;
    }
    double cluster1 = hist1DResoEtabins3D[ieta][iebin]->Integral(hist1DResoEtabins3D[ieta][iebin]->GetXaxis()->FindBin(0.5),hist1DResoEtabins3D[ieta][iebin]->GetXaxis()->FindBin(1.5));
    double cluster2 = hist1DResoEtabins3D[ieta][iebin]->Integral(hist1DResoEtabins3D[ieta][iebin]->GetXaxis()->FindBin(1.5),hist1DResoEtabins3D[ieta][iebin]->GetXaxis()->FindBin(2.5));

    double prob = cluster1/(cluster1 + cluster2)*100;
    double err10 = sqrt(cluster1 + cluster2);
    double err1  = pow(err10/(cluster1 + cluster2),2);
    double err2  = cluster1/(cluster1*cluster1);
    double err = sqrt(err1*err1 + err2*err2)*prob*100;

    if(cluster1>0){
      histEProb[ieta]->SetBinContent(histEProb[ieta]->GetXaxis()->FindBin(ebin), prob);
      histEProb[ieta]->SetBinError(histEProb[ieta]->GetXaxis()->FindBin(ebin), err);
      histEProb[ieta]->GetYaxis()->SetRangeUser(0,100);  
    }
    ebin = ebin + 2;
    iebin++;
    }

    SetStyleHistoTH1ForGraphs( histEProb[ieta],"", "E^{MC}", "Probability of Merging (%)", 0.75*textSizeSinglePad,0.85*textSizeSinglePad, 0.75*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.1);
    DrawGammaSetMarker(histEProb[ieta],32+ieta,2, kBlack, kBlack); 

    if(caloNameRead == "FEMC"){
      legendprob->AddEntry(histEProb[ieta], Form("%0.2f < #eta < %0.2f",etabinfwd[ieta], etabinfwd[ieta+1]) , "p");
    }else if(caloNameRead == "BECAL"){
      legendprob->AddEntry(histEProb[ieta], Form("%0.2f < #eta < %0.2f",etabincentral[ieta], etabincentral[ieta+1]) , "p"); 
    }else if(caloNameRead == "EEMC"){
      legendprob->AddEntry(histEProb[ieta], Form("%0.2f < #eta < %0.2f",etabinbck[ieta], etabinbck[ieta+1]) , "p"); 
    }else{
      cout << "check your calorimeter" << endl;

    }
    if(ieta==0) histEProb[ieta]->Draw();
    else histEProb[ieta]->Draw("same");
   }

   legendprob->Draw();
   ProbCanv->Print(Form("%s/Clusters_vs_E_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));


   //============= End Clusters Number ======================================

  TString NamepT_Minv = Form("%s/h_pT_Minv_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histpTMinv    = (TH2F*)inputFile->Get(NamepT_Minv.Data());  
  vector<double> allmeanPt;
  vector<double> allsigmaPt;


  if (histpTMinv->GetEntries() > 100){
    TCanvas* Minvpi02DPt = new TCanvas("Minvpi02DPt","",0,0,800,800);
    DrawGammaCanvasSettings( Minvpi02DPt, 0.15, 0.15, 0.1, 0.1);
    histpTMinv->GetXaxis()->SetRangeUser(0.,20);
    histpTMinv->GetYaxis()->SetRangeUser(0.,0.2);
    SetStyleHistoTH2ForGraphs(histpTMinv,"p_{T}^{MC}", "M_{#gamma #gamma}", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histpTMinv->Draw("colz");      
    TLatex Tl1Pt;
    Tl1Pt.SetTextSize(0.035);
    Tl1Pt.DrawLatex(0.4*energymax/2,0.04,perfLabel);
    Tl1Pt.DrawLatex(0.4*energymax/2,0.03,Form("Calorimeter: %s ", caloName.Data()));
    Tl1Pt.DrawLatex(0.4*energymax/2,0.02,Form("Clusterization type: %s ", clusterizerName.Data()));
    Tl1Pt.DrawLatex(0.4*energymax/2,0.01,"#pi ^{0} #rightarrow #gamma #gamma");
    Tl1Pt.Draw();
    Minvpi02DPt->Print(Form("%s/Mass_vs_Pt_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    while(Ptbin < 20){
      histPtMinvbin[padptbin]    = (TH1F*)histpTMinv->ProjectionY(Form("projection_Minv_pT_%f", Ptbin), 
                                             histpTMinv->GetXaxis()->FindBin(Ptbin-0.5), histpTMinv->GetXaxis()->FindBin(Ptbin+0.5),"e");

      if(histPtMinvbin[padptbin]->GetEntries() > 10){

        mean = meanErr = sigma = sigmaErr = 0;

        if (doFitting){   

          fithistPtMinv[padptbin]    = new TF1(Form("fit_Minv_Pt_%f", Ptbin), "crystalball", 0.05, 0.15);
          fithistPtMinv[padptbin]->SetParameters(0.55*histPtMinvbin[padptbin]->GetMaximum(),histPtMinvbin[padptbin]->GetMean(),
                  2*histPtMinvbin[padptbin]->GetRMS(),2.*histPtMinvbin[padptbin]->GetRMS(),1*histPtMinvbin[padptbin]->GetRMS());
          histPtMinvbin[padptbin]->Fit(fithistPtMinv[padptbin],"L0RMEQ","",0.9*minmassfit, 1.1*maxmassfit);
          fithistPtMinvFINAL[padptbin] = histPtMinvbin[padptbin]->GetFunction(Form("fit_Minv_Pt_%f", Ptbin));
             
          mean     = fithistPtMinvFINAL[padptbin]->GetParameter(1);
          meanErr  = fithistPtMinvFINAL[padptbin]->GetParError(1);
          sigma    = fithistPtMinvFINAL[padptbin]->GetParameter(2);
          sigmaErr = fithistPtMinvFINAL[padptbin]->GetParError(2);
          allmeanPt.push_back(mean);
          allsigmaPt.push_back(sigma);
          
        }

        Double_t reso     = sigma/mean;
        Double_t resoErr  = TMath::Sqrt( TMath::Power(sigmaErr/sigma,2) + TMath::Power(meanErr/mean,2) );
        histResoPtVsMinv->SetBinContent(histResoPtVsMinv->FindBin(Ptbin),  100*reso);
        histResoPtVsMinv->SetBinError(histResoPtVsMinv->FindBin(Ptbin), 100*resoErr);
        histMeanPtVsMinv->SetBinContent(histMeanPtVsMinv->FindBin(Ptbin), mean);
        histMeanPtVsMinv->SetBinError(histMeanPtVsMinv->FindBin(Ptbin), meanErr);
        histSigmaPtVsMinv->SetBinContent(histSigmaPtVsMinv->FindBin(Ptbin), sigma/sigmat);
        histSigmaPtVsMinv->SetBinError(histSigmaPtVsMinv->FindBin(Ptbin), sigmaErr/sigmat);
        padptbin++;
      }
      Ptbin = Ptbin + 1.; 
      
    }
 

    TCanvas* mass2DPt = new TCanvas("mass2DPt","",0,0,1000,800);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFontSize(0.07); 
    DrawGammaCanvasSettings( mass2DPt, 0.095, 0.01, 0.01, 0.10);
    mass2DPt->Divide(4,5);
    int padnumberPt = 0;
    TLegend* legend1Pt = new TLegend(0., 0.01, 0.4, 0.25);
    for(int z=0; z< padptbin +1 ; z++ ){
    if(z==3){
        mass2DPt->cd(z+1);
        TLatex TlPt;
        TlPt.SetTextSize(0.1);
        TlPt.DrawLatex(0.03,0.9,perfLabel);
        TlPt.DrawLatex(0.03,0.7,Form("Calorimeter: %s ", caloName.Data()));
        TlPt.DrawLatex(0.03,0.5,Form("Clusterization type: %s ", clusterizerName.Data()));
        TlPt.DrawLatex(0.03,0.3,"#pi ^{0} #rightarrow #gamma #gamma");

        legend1Pt->Draw();
        
      } else{
          mass2DPt->cd(z+1);
          SetStyleHistoTH1ForGraphs( histPtMinvbin[padnumberPt],Form("%0.1d GeV < p_{T} <  %0.1d GeV",z, z+1), "M_{#gamma #gamma}", "", 1.5*textSizeSinglePad,1.5*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,1.1);
          DrawGammaSetMarker(histPtMinvbin[padnumberPt],34,1, kBlack, kBlack);
          histPtMinvbin[padnumberPt]->GetYaxis()->SetNdivisions(5);
          histPtMinvbin[padnumberPt]->GetXaxis()->SetNdivisions(5);
          histPtMinvbin[padnumberPt]->GetXaxis()->SetRangeUser(0,0.3);
          histPtMinvbin[padnumberPt]->Draw("ep");
          fithistPtMinvFINAL[padnumberPt]->SetLineColor(kBlue);
          fithistPtMinvFINAL[padnumberPt]->SetLineWidth(1);
          fithistPtMinvFINAL[padnumberPt]->Draw("same");
          DrawGammaLines(allmeanPt[padnumberPt], allmeanPt[padnumberPt], 0, 0.5*fithistPtMinvFINAL[padnumberPt]->GetMaximum(),2, 2);
          DrawGammaLines(allmeanPt[padnumberPt] + allsigmaPt[padnumberPt], allmeanPt[padnumberPt]+ allsigmaPt[padnumberPt], 0, fithistPtMinvFINAL[padnumberPt]->GetMaximum(), 1, 12);
          DrawGammaLines(allmeanPt[padnumberPt] - allsigmaPt[padnumberPt], allmeanPt[padnumberPt]- allsigmaPt[padnumberPt], 0, fithistPtMinvFINAL[padnumberPt]->GetMaximum(), 1, 12);
          if(z==0){
           legend1Pt->AddEntry(histPtMinvbin[padnumberPt], "Simulation", "p");
           legend1Pt->AddEntry(fithistPtMinvFINAL[padnumberPt], "Fit", "l");
          }
        padnumberPt++;
     }
    }

    mass2DPt->Print(Form("%s/Mass_Pt_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    //============= Fit Results======================================
    TCanvas* fitresultsPt = new TCanvas("fitresultsPt","",0,0,1000,500);
    DrawGammaCanvasSettings( fitresultsPt, 0.15, 0.15, 0.15, 0.2);
   
    fitresultsPt->Divide(2,1);
    fitresultsPt->cd(1);

    histMeanPtVsMinv->GetYaxis()->SetRangeUser(0.08,0.16);
    SetStyleHistoTH1ForGraphs( histMeanPtVsMinv, "", "p_{T}^{MC}", "Mean (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.5);
    DrawGammaSetMarker(histMeanPtVsMinv,34,1, kBlack, kBlack);   
    histMeanPtVsMinv->GetYaxis()->SetRangeUser(0,1.3*histMeanPtVsMinv->GetMaximum()); 
    histMeanPtVsMinv->Draw("e,p");
    DrawGammaLines(0, energymax, mass_pi0PDG, mass_pi0PDG,2, 13);

    TLatex Tl3;
    Tl3.SetTextSize(0.035);
    Tl3.DrawLatex(2,0.18,perfLabel);
    Tl3.DrawLatex(2,0.17,Form("Calorimeter: %s ", caloName.Data()));
    Tl3.DrawLatex(2,0.16,Form("Clusterization type: %s ", clusterizerName.Data()));
    Tl3.DrawLatex(2,0.14,"#pi ^{0} #rightarrow #gamma #gamma");
    Tl3.Draw();
    fitresultsPt->cd(2);

    histSigmaPtVsMinv->GetYaxis()->SetRangeUser(0.,0.2);
    SetStyleHistoTH1ForGraphs(histSigmaPtVsMinv, "", "p_{T}^{MC}", "Width (M_{#gamma #gamma})", 0.8*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.3);  
    DrawGammaSetMarker(histSigmaPtVsMinv,34,1, kBlack, kBlack);   
    histSigmaPtVsMinv->GetYaxis()->SetRangeUser(0,0.015); 
    histSigmaPtVsMinv->Draw("e,p");
    
  
    fitresultsPt->Print(Form("%s/Mean_Width_Pt_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));

    //============= End Fit Results======================================


  }
  
  //============= M02 and M20 Results======================================
  TString Name_M02_Cluster1Pt = Form("%s/h_pT_M02_photon1_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02p1Pt    = (TH2F*)inputFile->Get(Name_M02_Cluster1Pt.Data());
  TString Name_M02_Cluster2Pt = Form("%s/h_pT_M02_photon2_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02p2Pt    = (TH2F*)inputFile->Get(Name_M02_Cluster2Pt.Data());
  TString Name_M02_SC_Pt = Form("%s/h_pT_M02_singlecluster_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM02SCPt   = (TH2F*)inputFile->Get(Name_M02_SC_Pt.Data());

  if (histM02p1Pt->GetEntries() > 100){

    TCanvas* M02DPt = new TCanvas("M02DPt","",0,0,1500,800);
    DrawGammaCanvasSettings( M02DPt, 0.1, 0.1, 0.01, 0.105);
    M02DPt->Divide(3,1);
    M02DPt->cd(1);
    histM02p1Pt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM02p1Pt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02p1Pt,"p_{T}^{MC}", "M02(#gamma 1 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02p1Pt->Draw("colz");
    M02DPt->cd(2);
    histM02p2Pt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM02p2Pt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02p2Pt,"p_{T}^{MC}", "M02(#gamma 2 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02p2Pt->Draw("colz");
    TPaveText *t2Pt = new TPaveText(0.2*energymax,0.7,0.45*energymax,0.9,"NB");
    t2Pt->SetFillColor(kWhite);
    t2Pt->AddText(perfLabel); ((TText*)t2Pt->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2Pt->AddText(Form("Calorimeter: %s ", caloName.Data()));  ((TText*)t2Pt->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2Pt->AddText(Form("Clusterization type: %s ", clusterizerName.Data()));  ((TText*)t2Pt->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2Pt->AddText("#pi ^{0} #rightarrow #gamma #gamma");  ((TText*)t2Pt->GetListOfLines()->Last())->SetTextColor(kBlack);
    t2Pt->Draw();
    M02DPt->cd(3);
    histM02SCPt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM02SCPt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02SCPt,"p_{T}^{MC}", "M02(Single Cluster)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02SCPt->Draw("colz");
     
    
    
    M02DPt->Print(Form("%s/M02_vs_Pt_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));
   }
   
  TString Name_M20_Cluster1zPt = Form("%s/h_pT_M20_photon1_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20p1Pt    = (TH2F*)inputFile->Get(Name_M20_Cluster1zPt.Data());
  TString Name_M20_Cluster2Pt = Form("%s/h_pT_M20_photon2_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20p2Pt    = (TH2F*)inputFile->Get(Name_M20_Cluster2Pt.Data());
  TString Name_M20_SCPt = Form("%s/h_pT_M20_singlecluster_%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
  histM20SCPt    = (TH2F*)inputFile->Get(Name_M20_SCPt.Data());

  if (histM02p1Pt->GetEntries() > 100){

    TCanvas* M20DPt = new TCanvas("M20DPt","",0,0,1500,800);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.05);
    DrawGammaCanvasSettings( M20DPt, 0.1, 0.1, 0.1, 0.105);
    M20DPt->Divide(3,1);
    M20DPt->cd(1);
    histM02p1Pt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM02p1Pt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM02p1Pt,"p_{T}^{MC}", "M20(#gamma 1 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM02p1Pt->Draw("colz");
    M20DPt->cd(2);
    histM20p2Pt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM20p2Pt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM20p2Pt,"p_{T}^{MC}", "M20(#gamma 2 in 2 clusters)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM20p2Pt->Draw("colz");
    TPaveText *tPt = new TPaveText(0.2*energymax,0.1,0.45*energymax,0.3,"NB");
    tPt->SetFillColor(kWhite);
    tPt->AddText(perfLabel); ((TText*)tPt->GetListOfLines()->Last())->SetTextColor(kBlack);
    tPt->AddText(Form("Calorimeter: %s ", caloName.Data()));  ((TText*)tPt->GetListOfLines()->Last())->SetTextColor(kBlack);
    tPt->AddText(Form("Clusterization type: %s ", clusterizerName.Data()));  ((TText*)tPt->GetListOfLines()->Last())->SetTextColor(kBlack);
    tPt->AddText("#pi ^{0} #rightarrow #gamma #gamma");  ((TText*)tPt->GetListOfLines()->Last())->SetTextColor(kBlack);
    tPt->Draw();
    M20DPt->cd(3);
    histM20SCPt->GetXaxis()->SetRangeUser(0.,0.5*energymax);
    histM20SCPt->GetYaxis()->SetRangeUser(0.,1);
    SetStyleHistoTH2ForGraphs(histM20SCPt,"p_{T}^{MC}", "M20(Single Cluster)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
    histM20SCPt->Draw("colz");
    
    
    
    M20DPt->Print(Form("%s/M20_vs_Pt_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));
  }

  //============= End M02 and M20 Results======================================

  //============= Clusters Number ======================================

   TString Name_Pt_eta_Clusters = Form("%s/hPt_eta_Clusters__%s_%s",caloNameRead.Data(), caloNameRead.Data(), clusterizerName.Data() );
   hist_Pt_eta_Clusters    = (TH3F*)inputFile->Get(Name_Pt_eta_Clusters.Data());

   TCanvas* ProbCanvPt = new TCanvas("ProbCanvPt","",0,0,800,800);
   DrawGammaCanvasSettings( ProbCanvPt, 0.095, 0.01, 0.1, 0.105);

   TLegend* legendprobPt = new TLegend(0.6, 0.15, 0.95, 0.4);
   legendprobPt->SetHeader(perfLabel);
   legendprobPt->AddEntry((TObject*)0, Form("Calorimeter: %s ", caloName.Data()),"");  
   legendprobPt->AddEntry((TObject*)0,Form("Clusterization type: %s ", clusterizerName.Data()),"");
   legendprobPt->AddEntry((TObject*)0, "#pi ^{0} #rightarrow #gamma #gamma",""); 

   for (int ieta = 0 ; ieta < etabine-1; ieta++){
  
   histPtProb[ieta]  = new TH1F(Form("histPtProb%i",ieta),"",100,0,41);
   ebin = 1;
   int iebin = 0;  
   while(ebin < 0.5*energymax){
    if(caloNameRead == "FEMC"){
      hist1DResoPttabins3D[ieta][iebin]    = (TH1F*)hist_Pt_eta_Clusters->ProjectionZ(Form("hist_pT_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin - binpT), hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin + binpT),
                                                                              hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabinfwd[ieta]),hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabinfwd[ieta+1]),"e");
    }else if(caloNameRead == "BECAL"){
      hist1DResoPttabins3D[ieta][iebin]    = (TH1F*)hist_Pt_eta_Clusters->ProjectionZ(Form("hist_pT_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin - binpT), hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin + binpT),
                                                                              hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabincentral[ieta]),hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabincentral[ieta+1]),"e");
    }else if(caloNameRead == "EEMC"){
      hist1DResoPttabins3D[ieta][iebin]    = (TH1F*)hist_Pt_eta_Clusters->ProjectionZ(Form("hist_pT_eta_Clusters_1D%i%i", iebin, iebin), 
                                                                              hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin - binpT), hist_Pt_eta_Clusters->GetXaxis()->FindBin(ebin + binpT),
                                                                              hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabinbck[ieta]),hist_Pt_eta_Clusters->GetYaxis()->FindBin(etabinbck[ieta+1]),"e");
    }else{
      cout << " Not a good calorimeter " << endl;
    }
    double cluster1 = hist1DResoPttabins3D[ieta][iebin]->Integral(hist1DResoPttabins3D[ieta][iebin]->GetXaxis()->FindBin(0.5),hist1DResoPttabins3D[ieta][iebin]->GetXaxis()->FindBin(1.5));
    double cluster2 = hist1DResoPttabins3D[ieta][iebin]->Integral(hist1DResoPttabins3D[ieta][iebin]->GetXaxis()->FindBin(1.5),hist1DResoPttabins3D[ieta][iebin]->GetXaxis()->FindBin(2.5));    
    
    double prob = cluster1/(cluster1 + cluster2)*100;
    double err10 = sqrt(cluster1 + cluster2);
    double err1  = pow(err10/(cluster1 + cluster2),2);
    double err2  = cluster1/(cluster1*cluster1);
    double err = sqrt(err1*err1 + err2*err2)*prob*100;

    if(cluster1>0){
      histPtProb[ieta]->SetBinContent(histPtProb[ieta]->GetXaxis()->FindBin(ebin), prob);
      histPtProb[ieta]->SetBinError(histPtProb[ieta]->GetXaxis()->FindBin(ebin), err);
      histPtProb[ieta]->GetYaxis()->SetRangeUser(0,100); 
      histPtProb[ieta]->GetXaxis()->SetRangeUser(0,0.5*energymax);  
    }
    ebin = ebin +ptskip;
    iebin++;
    }

    SetStyleHistoTH1ForGraphs( histPtProb[ieta],"", "p_{T}^{MC}", "Probability of Merging (%)", 0.75*textSizeSinglePad,0.85*textSizeSinglePad, 0.75*textSizeSinglePad,0.8*textSizeSinglePad,1.1,1.1);
    DrawGammaSetMarker(histPtProb[ieta],32+ieta,2, kBlack, kBlack); 

    if(caloNameRead == "FEMC"){
      legendprobPt->AddEntry(histPtProb[ieta], Form("%0.2f < #eta < %0.2f",etabinfwd[ieta], etabinfwd[ieta+1]) , "p");
    }else if(caloNameRead == "BECAL"){
      legendprobPt->AddEntry(histPtProb[ieta], Form("%0.2f < #eta < %0.2f",etabincentral[ieta], etabincentral[ieta+1]) , "p");
    }else if(caloNameRead == "EEMC"){
      legendprobPt->AddEntry(histPtProb[ieta], Form("%0.2f < #eta < %0.2f",etabinbck[ieta], etabinbck[ieta+1]) , "p");
    } else{
      cout << " Not a good calorimeter " << endl;
    }

    if(ieta==0) histPtProb[ieta]->Draw();
    else  histPtProb[ieta]->Draw("same");


   }

   legendprobPt->Draw();
   ProbCanvPt->Print(Form("%s/Clusters_vs_pT_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));


   //============= End Clusters Number ======================================

   //========Resolution ========================

  TCanvas* ResoCanv = new TCanvas("ResoCanv","",0,0,800,800);
  DrawGammaCanvasSettings( ResoCanv, 0.15, 0.1, 0.1, 0.15);
  SetStyleHistoTH1ForGraphs(histResoEVsMinv,"","E^{MC}", "M_{#gamma #gamma} Resolution (%)", 0.7*textSizeSinglePad,0.8*textSizeSinglePad, 0.7*textSizeSinglePad,0.8*textSizeSinglePad,1,1.5);
  DrawGammaSetMarker(histResoEVsMinv,33,2, kBlack, kBlack); 
  //histResoEVsMinv->GetYaxis()->SetRangeUser(0,60);
  histResoEVsMinv->GetYaxis()->SetRangeUser(0,1.3*histResoEVsMinv->GetMaximum());
  histResoEVsMinv->Draw();
  
  TLegend* legendresoPt = new TLegend(0.5, 0.2, 0.8, 0.4);
  legendresoPt->SetHeader(perfLabel);
  legendresoPt->AddEntry((TObject*)0, Form("Calorimeter: %s ", caloName.Data()),"");  
  legendresoPt->AddEntry((TObject*)0,Form("Clusterization type: %s ", clusterizerName.Data()),"");
  legendresoPt->AddEntry((TObject*)0, "#pi ^{0} #rightarrow #gamma #gamma",""); 
  legendresoPt->Draw();

  ResoCanv->Print(Form("%s/Reso_%s_%s.%s", outputDir.Data(), caloNameRead.Data(), clusterizerName.Data(), suffix.Data()));


  
  
  
  }
