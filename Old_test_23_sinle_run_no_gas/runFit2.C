#include <iostream>
#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TH2D.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TMath.h>
#include <TRandom3.h>

void example() {
    // Create a 2D histogram
    
   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->SetGrid();
   
    
   TFile* file = TFile::Open("recodata.root");
    //creating tree
    TTree * tree = (TTree*) file->Get("recodata");
    //init branches
    float x[60000];
    float y[60000];
    float t[60000];
    UShort_t n;
    // set branch address
    TGraphErrors* gr = new TGraphErrors(6000);
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("y",&y);
    tree->SetBranchAddress("t",&t);
    tree->SetBranchAddress("n",&n);
    //tree entry loop
  
     
     gr->SetMarkerColor(2);
     gr->SetMarkerSize(1);
     gr->SetMarkerStyle(20); 
      
     for(int i =11; i <12/*tree->GetEntries()*/; i++){
         tree->GetEntry(i);
         //per entry hits loop
         for (int k =0 ; k< n ; k++){
           gr->SetPoint(k,x[k],y[k]);
           gr->SetPointError(k,0.3,0.3);  
         }
         
    } 
  
   auto chi2Function = [&](const double *par) {
      //minimisation function computing the sum of squares of residuals
      // looping at the graph points
      int np = gr->GetN();
      double f = 0;
      double *x = gr->GetX();
      double *y = gr->GetY();
      double *ex = gr->GetEX();
      double *ey = gr->GetEY();
      for (int i=0;i<np;i++) {
         double u = x[i] - par[0];
         double v = y[i] - par[1];
         double dr = par[2] - std::sqrt(u*u+v*v);
         double err2 = ((u*u)*(ex[i]*ex[i]) + (v*v)*(ey[i]*ey[i]))/(u*u+v*v);
         f += dr*dr/err2;
      }
      return f;
   };
   // wrap chi2 funciton in a function object for the fit
   // 3 is the number of fit parameters (size of array par)
   ROOT::Math::Functor fcn(chi2Function,3);
   ROOT::Fit::Fitter  fitter;
   double pStart[3] = {0,0,20};
   fitter.SetFCN(fcn, pStart);
   fitter.Config().ParSettings(0).SetName("x0");
   fitter.Config().ParSettings(1).SetName("y0");
   fitter.Config().ParSettings(2).SetName("R");
   // do the fit 
   bool ok = fitter.FitFCN();
   if (!ok) {
      Error("line3Dfit","Line3D Fit failed");
   }   
   const ROOT::Fit::FitResult & result = fitter.Result();
   result.Print(std::cout);
   //Draw the circle on top of the points
   TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
   
   c1->DrawFrame(-100,-100,100,100);
   gr->Draw("p");
  
   arc->SetLineColor(kRed);
   arc->SetLineWidth(4);
   arc->SetFillStyle(0); 
   arc->Draw(); 
   c1->SaveAs("xygraph.png");
      
    
    
  
}
void runFit2() {
    example();
    
}
