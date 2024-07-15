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




TGraph* ConvertTH2DToTGraph(TH2D* hist) {
    if (!hist) return nullptr;

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();
    int nPoints = nBinsX * nBinsY;

    TGraph* graph = new TGraph(nPoints);
    int pointIndex = 0;

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            double binCenterX = hist->GetXaxis()->GetBinCenter(i);
            double binCenterY = hist->GetYaxis()->GetBinCenter(j);
            double binContent = hist->GetBinContent(i, j);

            graph->SetPoint(pointIndex, binCenterX, binContent);
            ++pointIndex;
        }
    }

    return graph;
}

TGraph* ConvertTH1ToTGraph(TH1* hist) {
    if (!hist) return nullptr;

    int nBins = hist->GetNbinsX();
    TGraph* graph = new TGraph(nBins);

    for (int i = 1; i <= nBins; ++i) {
        double binCenter = hist->GetBinCenter(i);
        double binContent = hist->GetBinContent(i);
        graph->SetPoint(i - 1, binCenter, binContent);
    }

    return graph;
}

// Define the circle function to be used in the fit
Double_t circleFunc(Double_t *v, Double_t *par) {
    
    
    
    Double_t x = v[0];
    Double_t y = v[1];
    Double_t x0 = par[0];
    Double_t y0 = par[1];
    Double_t r = par[2];
    // Residuals from the circle equation (x-x0)^2 + (y-y0)^2 = r^2
    return pow((x - x0), 2) + pow((y - y0), 2) - pow(r, 2);
}



// Function to fit a 2D histogram with a circle


void fitHistogramWithCircle(TH2D* hist,TH2D* hist2, TH1D* hist3) {
    // Define the fitting function
    TF2 *fitFunc = new TF2("fitFunc", circleFunc, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(),
                           hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(), 3);
    fitFunc->SetParNames("x0", "y0", "r");
    

    // Initial parameter guesses (center of the histogram and some reasonable radius)
    fitFunc->SetParameters(hist->GetMean(1), hist->GetMean(2), (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin()) / 4);
    fitFunc->SetParLimits(2,30,70);
    
    fitFunc->SetParLimits(0,-2,2);
    fitFunc->SetParLimits(1,-2,2);
    
    // Fit the histogram
    hist->Fit(fitFunc, "NR");

    // Print the results
    std::cout << "Fitted circle parameters:" << std::endl;
    std::cout << "x0 = " << fitFunc->GetParameter(0) << std::endl;
    std::cout << "y0 = " << fitFunc->GetParameter(1) << std::endl;
    std::cout << "r = " << fitFunc->GetParameter(2) << std::endl;
    TEllipse *circle3 = new TEllipse( fitFunc->GetParameter(0), fitFunc->GetParameter(1),fitFunc->GetParameter(2), fitFunc->GetParameter(2));
    hist2->Fill( fitFunc->GetParameter(0), fitFunc->GetParameter(1)); 
    hist3->Fill( fitFunc->GetParameter(2)); 
    //TArc *arc = new TArc(fitFunc->GetParameter(0),fitFunc->GetParameter(1),fitFunc->GetParameter(2),0,360);
    circle3->Draw("SAME");
}

void example() {
    // Create a 2D histogram
    
    float xmax= 100;
    float xmin =-100;
    float ymax= 80;
    float ymin =-80;

    int nbinx = 50;
    int nbiny = 50;
    TH2D *hist = new TH2D("hist", ";xhits;yhits",nbinx,xmax,xmin,nbiny,ymax,ymin);
    TH1D *histred = new TH1D("histred", ";radius;entry",100,xmax,xmin);
    TH2D *histcent = new TH2D("histcent", ";xcenter;ycenter",nbinx,xmax,xmin,nbiny,ymax,ymin);
    

    //opening the file
    
    TFile* file = TFile::Open("recodata.root");
    //creating tree
    TTree * tree = (TTree*) file->Get("recodata");
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
    //tree entry loop
     for(int i =0; i <tree->GetEntries(); i++){
         tree->GetEntry(i);
         //per entry hits loop
         for (int k =0 ; k< n ; k++){
            TH2D *nhistxy = new TH2D("nhistxy", ";xhits;yhits",nbinx,xmax,xmin,nbiny,ymax,ymin); 
            nhistxy->Fill(x[k],y[k]);
             hist->Fill(x[k],y[k]); 
            //fitHistogramWithCircle(hist);
            //Fit the histogram with a circle
            fitHistogramWithCircle(nhistxy,histcent, histred); 
         }
         
    }
    


    // Example data: fill the histogram with some data points that roughly form a circle
  // TH2D* hist = (TH2D*) nhistxy->Clone("");
    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Circle Fit Example", 800, 600);
    hist->Draw("COLZ");
    
    
    TCanvas *c2 = new TCanvas("c2", "Circle Fit Example", 800, 600);
    histred->Draw("COLZ");
    
    TCanvas *c3 = new TCanvas("c3", "Circle Fit Example", 800, 600);
    histcent->Draw("COLZ");

   

    c1->SaveAs("circle_fit.png");
    c2->SaveAs("rad_fit.png");
    c3->SaveAs("cent_fit.png");
}




void runFit() {
    example();
    
}
