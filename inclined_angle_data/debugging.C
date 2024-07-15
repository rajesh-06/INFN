// #pragma once
#include <vector>
#include <TCanvas.h>
#include <iostream>
#include "lib/sensor_scan_lib.C"
#include <TString.h>
#include <TAxis.h>

using namespace std;

void 
debugging(){
    TCanvas *c1 = new TCanvas("c1", "All Data", 1200, 400);
    c1->Divide(3, 1);

    vector <string> measurements = {//6 measurements
        "40deg",
        "40deg-a",
        "40deg-b",
        "40deg-c",
        "40deg-d",
        "40deg-e",
    };
    vector<double> x(measurements.size());
    for (size_t i = 0; i < measurements.size(); ++i) {
        x[i] = i + 1; // Assigning numeric x values (1, 2, 3, ...)
    }

    vector <string> file_end = {
        ".scanx.yoff.chip-0.channel-A1.txt.tree.root",//reference sensor
        ".scanx.yoff.chip-1.channel-A1.txt.tree.root",//sensor A1
        ".scanx.yoff.chip-1.channel-A2.txt.tree.root",//sensor A2
        ".scanx.yoff.chip-1.channel-A3.txt.tree.root",//sensor A3
        ".scanx.yoff.chip-1.channel-A4.txt.tree.root",//sensor A4
    };
    // vector < vector <TGraphErrors*> > gr[][];
    vector < vector <TF1*> > fit(file_end.size());
    vector < vector <string> > filenames(file_end.size());//[sensor][measurements]


    for (int sensor = 0; sensor < file_end.size(); ++sensor){
        for (int dir = 0; dir < measurements.size(); ++dir){
            auto file = measurements[dir]+"/"+measurements[dir]+file_end[sensor];
            cout<<file<<"\t";
            auto gr = get_sensor_scan(file);
            auto fit_gr = get_fermi_fit(gr);
            filenames[sensor].push_back(file);
            fit[sensor].push_back(fit_gr);
        }cout<<endl;
    }
    
    vector <TGraphErrors*> gr(file_end.size());
   // auto gr = new TGraphErrors();

    for (int j = 0; j < fit.size(); ++j){
        gr[j] = new TGraphErrors(measurements.size());
        for (int i = 0; i< fit[0].size(); ++i){
            cout<<filenames[j][i]<<endl;
            cout<<fit[j][i]->GetParameter(1)<<endl;
            gr[j]->SetPoint(i, x[i], fit[j][i]->GetParameter(1));

        }
    }
    
    
    c1->cd(1);
    gr[0]->Draw("ap");

    gr[0]->SetTitle("; ;Rate (Hz)");
    gr[0]->SetMarkerStyle(21); // Marker style
    gr[0]->SetMarkerSize(1.2); // Marker size
    gr[0]->SetMarkerColor(kMagenta+2); // Marker color
    //gr_pde->SetLineColor(kCyan); // Line color
    // gr[0]->SetLineWidth(2); // Line width
    for (int i = 1; i < gr.size(); ++i){
        
        gr[i]->SetMarkerSize(1); // Marker size
        gr[i]->SetMarkerColor(20+i); // Marker color
        gr[i]->Draw("P SAME");
    }
    
    TAxis *axis = gr[0]->GetXaxis();
    axis->SetNdivisions(x.size(), kFALSE);
    for (size_t i = 0; i < measurements.size(); ++i) {
        axis->SetBinLabel(axis->FindBin(x[i]), measurements[i].c_str());
    }

    // gr_pde->GetXaxis()->SetTitleSize(0.04); // X axis title size
    // gr_pde->GetYaxis()->SetTitleSize(0.04); // Y axis title size
    // gr_pde->GetXaxis()->SetLabelSize(0.04); // X axis label size
    // gr_pde->GetYaxis()->SetLabelSize(0.04); // Y axis label size
    // gr_pde->GetXaxis()->SetTitleOffset(1.2); // X axis title offset
    // gr_pde->GetYaxis()->SetTitleOffset(1.5); // Y axis title offset
    
    c1->Update();
    // c1->Draw();

}