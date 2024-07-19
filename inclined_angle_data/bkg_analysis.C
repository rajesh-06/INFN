#pragma once
#include <vector>
#include <TCanvas.h>
#include <iostream>
#include <TString.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include "lib/sensor_scan_lib.C"

using namespace std;

void bkg_analysis() {
    TCanvas *c1 = new TCanvas("c1", "All Data", 800, 600); // Increased canvas height for better visibility

    vector<string> measurements = {//   6 measurements
        "40deg",
        "40deg-a",
        "40deg-b",
        "40deg-c",
        "40deg-d",
        "40deg-e",
        "40deg-f",
        
    };

    vector <string> time = {
        "08:56", "14:06", "16:03", "18:07", "20:17", "09:23", "11:34",
    };

    vector <double> x = {
        536, 846, 963, 1087, 1217, 2003, 2134,
    };
    vector<string> file_end = {
        ".scanx.yoff.chip-0.channel-A1.txt.tree.root",//reference sensor
        ".scanx.yoff.chip-1.channel-A1.txt.tree.root",//sensor A1
        ".scanx.yoff.chip-1.channel-A2.txt.tree.root",//sensor A2
        ".scanx.yoff.chip-1.channel-A3.txt.tree.root",//sensor A3
        ".scanx.yoff.chip-1.channel-A4.txt.tree.root",//sensor A4
    };

    vector<string> sensor_names = {
        "Reference Sensor",
        "Sensor A1",
        "Sensor A2",
        "Sensor A3",
        "Sensor A4"
    };

    vector<vector<TF1*>> fit(file_end.size());
    vector<vector<string>> filenames(file_end.size()); //[sensor][measurements]

    // Initialize min and max y-values
    double minY = numeric_limits<double>::max();
    double maxY = -numeric_limits<double>::max();

    for (int sensor = 0; sensor < file_end.size(); ++sensor) {
        for (int dir = 0; dir < measurements.size(); ++dir) {
            auto file = measurements[dir] + "/" + measurements[dir] + file_end[sensor];
            cout << file << "\t";
            auto gr = get_sensor_scan(file);
            auto *fit_gr = get_fermi_fit(gr);
            filenames[sensor].push_back(file);
            fit[sensor].push_back(fit_gr);

            // Find min and max y-values
            // for (int i = 0; i < gr->GetN(); ++i) {
                float y = fit_gr->GetParameter(0);
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            // }
        }
        cout << endl;
    }

    vector<TGraphErrors*> gr(file_end.size());

    for (int j = 0; j < fit.size(); ++j) {
        gr[j] = new TGraphErrors(measurements.size());
        for (int i = 0; i < fit[0].size(); ++i) {
            cout << filenames[j][i] << endl;
            cout << fit[j][i]->GetParameter(0) << endl;
            gr[j]->SetPoint(i, x[i], fit[j][i]->GetParameter(0));
            gr[j]->SetPointError(i, 0, fit[j][i]->GetParError(0));

        }
    }

    // Create a legend
    TLegend *legend = new TLegend(0.8, 0.1, 1., 0.35);
    
    // Define unique marker colors for each graph
    int colors[] = { kRed, kGreen-4, kOrange,  kCyan-4, kViolet+1};

    for (int i = 0; i < gr.size(); ++i) {
        if (i == 0) {
            // c1->cd();
            // Set y-axis range for the first graph
            gr[i]->SetMinimum(minY - 0.1 * (maxY - minY)); // Adjust buffer
            gr[i]->SetMaximum(maxY + 0.1 * (maxY - minY)); // Adjust buffer
            gr[i]->Draw("ALP");
            gr[i]->SetTitle(";;Background Rate (Hz)");
        } else {
            // Set y-axis range for subsequent graphs
            gr[i]->GetYaxis()->SetRangeUser(minY - 0.1 * (maxY - minY), maxY + 0.1 * (maxY - minY)); // Adjust buffer
            gr[i]->Draw("SAMELP");
        }

        gr[i]->SetMarkerStyle(21); // Marker style
        gr[i]->SetMarkerSize(0.8); // Marker size
        gr[i]->SetMarkerColor(colors[i % (sizeof(colors) / sizeof(colors[0]))]); // Marker color
        gr[i]->SetLineColor(colors[i]); // Line color
        gr[i]->SetLineWidth(2); // Line width
        legend->AddEntry(gr[i], sensor_names[i].c_str(), "LPE");
    }

    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->Draw();

    TAxis *axis = gr[0]->GetXaxis();
    axis->SetNdivisions(x.size(), kFALSE);
    axis->SetLabelSize(0.05); // Adjust axis label size

    for (size_t i = 0; i < measurements.size(); ++i) {
        axis->SetBinLabel(axis->FindBin(x[i]), time[i].c_str());
    }
    c1->SetLeftMargin(0.1362725);
    c1->SetBottomMargin(0.1316309);
    c1->Update();
}
