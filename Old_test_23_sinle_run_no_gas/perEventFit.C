#include "iostream"
#include "fstream"

using namespace std;

/*
amrit gautam 

agautam@cern.ch


*/
void perEventFit(const char *filename = "recodata.root"){

    TFile* file = TFile::Open(filename);
    
    if(!file){
     cout << "can not find the file" << endl;
     return;   
        
    }
    //creating tree
    TH1D* rad = new TH1D("rad", "", 200, 0,100);
    TH2D* cent = new TH2D("cent", "", 200, -20,20, 200,-20,20); 
    
    TTree * tree = (TTree*) file->Get("recodata");
    
    if(!tree){
       cout << "can not read the tree" << endl;
       return;     
    }
    //init branches
    float x[60000];
    float y[60000];
    float t[60000];
    UShort_t n;
    
    // set branch address
   
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("y",&y);
    tree->SetBranchAddress("t",&t);
    tree->SetBranchAddress("n",&n);
    
    TCanvas*  arccan = new TCanvas("arccan" ,"",600,600);
    arccan->cd();
    TArc *arc2 = new TArc(0,0,70);
     
    //tree entry loop
    for(int i =0; i <tree->GetEntries(); i++){
        //cout<< "getting entries"<< endl;
        tree->GetEntry(i);
        TGraphErrors* gr = new TGraphErrors(n);
        //per entry hits loop
        
         for (int k =0 ; k< n ; k++){
           if(TMath::Abs(t[k]>30)) continue;  
           gr->SetPoint(k,x[k],y[k]);
           gr->SetPointError(k,1.5,1.5);  
          }//end of n hits
         
        
        //defining chisq minimization
        auto chi2Function = [&](const double *par) {
      //minimisation function computing the sum of squares of residuals
      // looping at the graph points
          int np = gr->GetN();
    //  cout<< np<< endl;
          double f = 0;
          double *x = gr->GetX();
          double *y = gr->GetY();
          double *ex = gr->GetEX();
          double *ey = gr->GetEY();
          for (int i=0;i<np;i++) {
             double u = x[i] - par[0];
             double v = y[i] - par[1];
             if( (u ==0 ) || (v==0) )continue;
             double dr = par[2] - std::sqrt(u*u+v*v);
             double err2 = ((u*u)*(ex[i]*ex[i]) + (v*v)*(ey[i]*ey[i]))/(u*u+v*v);
             //cout << err2<< endl;
             f += dr*dr/err2;
             //cout <<f<< endl;
          }
          return f;
       };

       ROOT::Math::Functor fcn(chi2Function,3);
       
       ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100);
       //gROOT->SetPrintLevel(2);
       //ROOT::Math::Minimizer::MaxFunctionCalls(100);
       ROOT::Fit::Fitter  fitter;
  
   double pStart[3] = {0.,0.,40};
  // fitter.F
        
   fitter.SetFCN(fcn, pStart);
   //fitter.FixParameter(0,1);      
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
   rad->Fill(result.Parameter(2));
   cent->Fill(result.Parameter(0),result.Parameter(1));   
   TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
   //arc.SetNPX(1000);
   arc->Draw("SAME");  
        

         
    } //end of entry
    
    arccan->SaveAs("radiusofcircles.png");


    
    TCanvas *can = new TCanvas("can","",1600,800);
    can->cd();
    can->Divide(2,1);
    can->cd(1);
    rad->Draw();
    can->cd(2);
    cent->Draw();
    can->Draw();
    can->SaveAs("circ_cent_and_rad_dist.png");
    
    
    //can->DrawFrame(-100,-100,100,100);

    //arc->SetFillStyle(0);
    //arc->SetLineColor(34);
    //arc->Draw("SAME");
    //gr->Draw("Z");


}