#pragma once
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TString.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include "TH1F.h"
#include "TLine.h"
#include <TLegend.h>
#include "../lib/sensor_scan_lib.C"


using namespace std;


TGraphErrors* get_graph (vector <string> filenames, int par, string scan_axis = "x"){    
        auto *gr = new TGraphErrors(filenames.size());
        for (int i = 0; i< filenames.size(); ++i){
            auto scan = get_sensor_scan(filenames[i], scan_axis);
            auto fit = get_fermi_fit(scan);
            gr->SetPoint(i, i+1 , fit->GetParameter(par));
            gr->SetPointError(i, 0 , fit->GetParError(par) );      
        }
        return gr;
}

TGraphErrors* get_cen_graph (vector <string> filenames, int par, string scan_axis = "x"){    
        auto *gr = new TGraphErrors(filenames.size());
        for (int i = 0; i< filenames.size(); ++i){
            auto scan = get_sensor_scan(filenames[i], scan_axis);
            auto fit = get_scan_fit(scan);
            gr->SetPoint(i, i+1 , fit->GetParameter(par));
            gr->SetPointError(i, 0 , fit->GetParError(par) );      
        }
        return gr;
}


TGraph* PlotPulls(TGraphErrors *graph) {
  // Check if valid TGraphErrors
  if (!graph || graph->GetN() <= 0) {
    std::cerr << "Invalid TGraphErrors object!" << std::endl;
    return nullptr;
  }

  // Get original data
  int npoints = graph->GetN();
  double *x = graph->GetX();
  double *y = graph->GetY();
//   double *ex = graph->GetEX();
  double *ey = graph->GetEY();

  // Fit the graph (replace with your fitting function)
  // Here, we assume a simple linear fit (y = a + b*x)
  TF1 *fit = new TF1("fit", "[0] + [1]*x", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
  graph->Fit(fit);

  // Create a histogram for residuals (pulls)
  TGraph *pulls_hist = new TGraph(npoints);

  // Fill the histogram with pulls (y - fit result) and errors
  for (int i = 0; i < npoints; ++i) {
    double pull = (y[i] - fit->Eval(x[i])) / ey[i];
    pulls_hist->SetPoint(i, x[i], pull);
  }

  // Set plotting options for the histogram (optional)
//   pulls_hist->SetLineColor(kBlue);
//   pulls_hist->SetMarkerColor(kBlue);
//   pulls_hist->SetMarkerStyle(20);
  pulls_hist->GetXaxis()->SetTitle("X");
  pulls_hist->GetYaxis()->SetTitle("Pull [(y - fit) / sigma]");

  // Add a horizontal line at y = 0 (mean pull)
  TLine *zero_line = new TLine(pulls_hist->GetXaxis()->GetXmin(), 0, pulls_hist->GetXaxis()->GetXmax(), 0);
  zero_line->SetLineColor(kRed);
  zero_line->SetLineStyle(2);

  // Draw the pull distribution and zero line
  pulls_hist->Draw();
  zero_line->Draw("same");

    return pulls_hist;
  // Clean up (optional)
  delete fit;
  delete zero_line;
}


