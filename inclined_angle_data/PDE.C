#include "lib/sensor_scan_lib.C"
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

TGraphErrors* pde(string fileinitial = "flat/flat.scanx.yoff", string scan_axis = "x", bool dip = 0) {
    vector<string> filenames = {
        fileinitial + ".chip-0.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A2.txt.tree.root",
        fileinitial + ".chip-1.channel-A3.txt.tree.root",
        fileinitial + ".chip-1.channel-A4.txt.tree.root"
    };

    TGraphErrors* rates = new TGraphErrors(filenames.size()-1);

    auto ref_gr = get_sensor_scan(filenames[0], scan_axis);
    if (!ref_gr) {
        cerr << "Error: Could not get reference graph from file " << filenames[0] << endl;
        return nullptr;
    }

    TF1* ref_res = dip ? get_scan_fit(ref_gr) : get_fermi_fit(ref_gr);
    if (!ref_res) {
        cerr << "Error: Could not get fit for reference graph from file " << filenames[0] << endl;
        return nullptr;
    }

    double ref_rate = ref_res->GetParameter(1);
    double err_ref_rate = ref_res->GetParError(1);

    for (size_t j = 1; j < filenames.size(); ++j) {
        auto graph = get_sensor_scan(filenames[j], scan_axis);
        if (!graph) {
            cerr << "Error: Could not get graph from file " << filenames[j] << endl;
            continue;
        }

        TF1* result = dip ? get_scan_fit(graph) : get_fermi_fit(graph);
        if (!result) {
            cerr << "Error: Could not get fit for graph from file " << filenames[j] << endl;
            continue;
        }
        cout << result->GetParameter(1) << "\t" << result->GetParError(1) << "\t" << ref_res->GetParameter(1) << "\t" << ref_res->GetParError(1) << endl;
        rates->SetPoint(j-1, j, result->GetParameter(1)/ref_res->GetParameter(1));
        double err = propagate_error_division(result->GetParameter(1), result->GetParError(1), ref_res->GetParameter(1), ref_res->GetParError(1));
        cout << err << endl;
        rates->SetPointError(j-1, 0, err);
    }

    vector<string> x = {"A1", "A2", "A3", "A4"};
    TAxis *axis = rates->GetXaxis();
    axis->SetNdivisions(x.size(), kFALSE);
    axis->SetLabelSize(0.05); // Adjust axis label size

    for (size_t k = 0; k < x.size(); ++k) {
        axis->SetBinLabel(axis->FindBin(k+1), x[k].c_str());
    }

    return rates;
}

void PDE(){
    auto can = new TCanvas("can", "Rate", 800, 600);

    auto *flat = pde("flat/flat.scanx.yoff", "x", 0);
    auto *flatv2 = pde("flat-v2.1/flat-v2.1.scanx.yoff", "x", 0);
    auto *flatv2b = pde("flat-v2.1-b/flat-v2.1-b.scanx.yoff", "x", 0);

    auto *flat_ycen = pde("flat/flat.scanx.ycen", "x", 1);
    auto *flatv2_ycen = pde("flat-v2.1/flat-v2.1.scanx.ycen", "x", 1);
    auto *flatv2b_ycen = pde("flat-v2.1-b/flat-v2.1-b.scanx.ycen","x", 1);

    auto *flat_xcen = pde("flat/flat.scany.xcen", "y", 1);
    auto *flatv2_xcen = pde("flat-v2.1/flat-v2.1.scany.xcen", "y", 1);
    auto *flatv2b_xcen = pde("flat-v2.1-b/flat-v2.1-b.scany.xcen", "y", 1);

    flat->SetTitle(" ; Sensor ; PDE"); // scanx and yoff data

    flat->SetMarkerStyle(20);
    flatv2->SetMarkerStyle(21);
    flatv2b->SetMarkerStyle(22);

    flat->SetMarkerColor(kRed);
    flatv2->SetMarkerColor(kBlue);
    flatv2b->SetMarkerColor(kGreen);

    flat->SetLineColor(kRed);
    flatv2->SetLineColor(kBlue);
    flatv2b->SetLineColor(kGreen);

    flat_ycen->SetMarkerStyle(20);
    flatv2_ycen->SetMarkerStyle(21);
    flatv2b_ycen->SetMarkerStyle(22);

    flat_ycen->SetMarkerColor(kRed+2);
    flatv2_ycen->SetMarkerColor(kBlue+2);
    flatv2b_ycen->SetMarkerColor(kGreen+2);

    flat_ycen->SetLineColor(kRed+2);
    flatv2_ycen->SetLineColor(kBlue+2);
    flatv2b_ycen->SetLineColor(kGreen+2);

    flat_xcen->SetMarkerStyle(20);
    flatv2_xcen->SetMarkerStyle(21);
    flatv2b_xcen->SetMarkerStyle(22);

    flat_xcen->SetMarkerColor(kRed-2);
    flatv2_xcen->SetMarkerColor(kBlue-2);
    flatv2b_xcen->SetMarkerColor(kGreen-2);

    flat_xcen->SetLineColor(kRed-2);
    flatv2_xcen->SetLineColor(kBlue-2);
    flatv2b_xcen->SetLineColor(kGreen-2);

    flat->SetMaximum(1.02);
    flat->SetMinimum(0.87);
    can->cd();
    flat->Draw("ALP");
    flatv2->Draw("SAMELP");
    flatv2b->Draw("SAMELP");

    flat_ycen->Draw("SAMELP");
    flatv2_ycen->Draw("SAMELP");
    flatv2b_ycen->Draw("SAMELP");

    flat_xcen->Draw("SAMELP");
    flatv2_xcen->Draw("SAMELP");
    flatv2b_xcen->Draw("SAMELP");

    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(flat, "flat-yoff", "lpfe");
    legend->AddEntry(flat_xcen, "flat-xcen", "lpfe");
    legend->AddEntry(flat_ycen, "flat-ycen", "lpfe");

    legend->AddEntry(flatv2, "flat-v2-yoff", "lpfe");
    legend->AddEntry(flatv2_xcen, "flat-v2.1-xcen", "lpfe");
    legend->AddEntry(flatv2_ycen, "flat-v2.1-ycen", "lpfe");

    legend->AddEntry(flatv2b, "flat-v2.1-b-yoff", "lpfe");
    legend->AddEntry(flatv2b_xcen, "flat-v2.1-b-xcen", "lpfe");
    legend->AddEntry(flatv2b_ycen, "flat-v2.1-b-ycen", "lpfe");

    legend->Draw();

    can->Draw();
    can->SaveAs("results/pde_flat_all.png");
}
