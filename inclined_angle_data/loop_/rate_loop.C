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

TGraphErrors* get_graph (vector <string> filenames, int par){    
        auto *gr = new TGraphErrors(filenames.size());
        for (int i = 0; i< filenames.size(); ++i){
            auto scan = get_sensor_scan(filenames[i]);
            auto fit = get_fermi_fit(scan);
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


void rate_loop() {
    TCanvas *c1 = new TCanvas("c1", "rate analysis", 1200, 800); // Increased canvas height for better visibility


    std::vector<std::string> filenames;
    std::vector<string> measurements;
    vector <float> x;

    // Loop from 0 to 33 (inclusive)
    for (int i = 0; i <= 99; ++i) {
        auto file_end = ".scanx.yoff.chip-0.channel-A1.txt.tree.root";
        std::stringstream filename_stream;

        // Construct the filename with padding
        filename_stream << "loop-";
        filename_stream.width(2);
        filename_stream.fill('0');
        filename_stream << i;

    // Add the filename to the vector
    measurements.push_back(filename_stream.str());
    filenames.push_back(filename_stream.str()+"/"+filename_stream.str()+file_end);
    cout<<filenames[i]<<endl;
    cout<<measurements[i]<<endl;
    x.push_back(i+1);
    
    }

    

    auto *gr_low = get_graph(filenames, 2);
    auto *gr_high = get_graph(filenames, 3);

    auto gr = new TGraphErrors(gr_low->GetN());
    for (int i = 0; i < gr_low->GetN(); i++){
        gr->SetPoint(i, i+1, gr_high->GetPointY(i) - gr_low->GetPointY(i) );
        auto error = propagate_error_add_sub({gr_high->GetErrorY(i), gr_low->GetErrorY(i) });
        gr->SetPointError(i, 0, error );

    }

    //auto gr = PlotPulls(gr_a1);
  //string _data = " Photon rate detected by the sensor";
  string _data = " Sensor width";

  string title = "; loop number ;"+ _data+" [mm]";


  gr->SetTitle(title.c_str());
  gr->Draw("ALP");
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1);
  gr->SetMarkerColor(kRed);



  
    // auto *legend = new TLegend(0.7,0.7, 0.9, 0.9);

    // legend->SetBorderSize(0);
    // legend->SetFillColor(0);
    // legend->Draw();

    // TAxis *axis = gr->GetXaxis();
    // axis->SetNdivisions(measurements.size(), kFALSE);
    // axis->SetLabelSize(0.05); // Adjust axis label size

    // for (size_t i = 0; i < measurements.size(); ++i) {
    //     axis->SetBinLabel(axis->FindBin(x[i]), measurements[i].c_str());
    // }
    c1->SetLeftMargin(0.1362725);
    c1->SetBottomMargin(0.1316309);
    c1->Update();
    c1->Draw();
    auto outfile = "results/"+ _data +".png";
    c1->SaveAs(outfile.c_str());

}
