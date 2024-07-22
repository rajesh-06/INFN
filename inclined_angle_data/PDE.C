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


TGraphErrors* background(string fileinitial = "flat/flat.scanx.yoff", string scan_axis = "x", bool dip = 0) {
    vector<string> filenames = {
        fileinitial + ".chip-0.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A2.txt.tree.root",
        fileinitial + ".chip-1.channel-A3.txt.tree.root",
        fileinitial + ".chip-1.channel-A4.txt.tree.root"
    };

    TGraphErrors* bkg = new TGraphErrors(filenames.size()-1);

    for (size_t j = 0; j < filenames.size(); ++j) {
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
        //cout << result->GetParameter(1) << "\t" << result->GetParError(1) << "\t" << ref_res->GetParameter(1) << "\t" << ref_res->GetParError(1) << endl;
        bkg->SetPoint(j, j, result->GetParameter(0));
        bkg->SetPointError(j, 0, result->GetParError(0));
    }

    vector<string> x = {"Reference","A1", "A2", "A3", "A4"};
    TAxis *axis = bkg->GetXaxis();
    axis->SetNdivisions(x.size(), kTRUE);
    axis->SetLabelSize(0.05); // Adjust axis label size
    // axis->LabelsOption("v");
    // axis->RotateTitle();//
    //axis->SetRange(-1, 25);
    for (size_t k = 0; k < x.size(); ++k) {
        axis->SetBinLabel(axis->FindBin(k), x[k].c_str());
    }
    bkg->GetXaxis()->LabelsOption("h");
    bkg->GetXaxis()->SetRangeUser(-1, 5);
    return bkg;
}

void PDE(){
    auto can = new TCanvas("can", "Rate", 1200, 900);

    auto *flat = pde("flat/flat.scanx.yoff", "x", 0);
    auto *flatv2 = pde("flat-v2.1/flat-v2.1.scanx.yoff", "x", 0);
    auto *flatv2b = pde("flat-v2.1-b/flat-v2.1-b.scanx.yoff", "x", 0);
    auto *flatv2c = pde("flat-v2.1-c/flat-v2.1-c.scanx.yoff", "x", 0);


    auto *flat_ycen = pde("flat/flat.scanx.ycen", "x", 1);
    auto *flatv2_ycen = pde("flat-v2.1/flat-v2.1.scanx.ycen", "x", 1);
    auto *flatv2b_ycen = pde("flat-v2.1-b/flat-v2.1-b.scanx.ycen","x", 1);
    auto *flatv2c_ycen = pde("flat-v2.1-c/flat-v2.1-c.scanx.ycen","x", 1);


    auto *flat_xcen = pde("flat/flat.scany.xcen", "y", 1);
    auto *flatv2_xcen = pde("flat-v2.1/flat-v2.1.scany.xcen", "y", 1);
    auto *flatv2b_xcen = pde("flat-v2.1-b/flat-v2.1-b.scany.xcen", "y", 1);
    auto *flatv2c_xcen = pde("flat-v2.1-c/flat-v2.1-c.scany.xcen", "y", 1);

    flat->SetTitle(" ; Sensor ; PDE"); // scanx and yoff data

    flat->SetMarkerStyle(20);
    flatv2->SetMarkerStyle(21);
    flatv2b->SetMarkerStyle(22);
    flatv2c->SetMarkerStyle(23);


    flat->SetMarkerColor(kRed);
    flatv2->SetMarkerColor(kBlue);
    flatv2b->SetMarkerColor(kGreen);
    flatv2c->SetMarkerColor(kOrange);


    flat->SetLineColor(kRed);
    flatv2->SetLineColor(kBlue);
    flatv2b->SetLineColor(kGreen);
    flatv2c->SetLineColor(kOrange);


    flat_ycen->SetMarkerStyle(20);
    flatv2_ycen->SetMarkerStyle(21);
    flatv2b_ycen->SetMarkerStyle(22);
    flatv2c_ycen->SetMarkerStyle(23);

    flat_ycen->SetMarkerColor(kRed+2);
    flatv2_ycen->SetMarkerColor(kBlue+2);
    flatv2b_ycen->SetMarkerColor(kGreen+2);
    flatv2c_ycen->SetMarkerColor(kOrange+2);


    flat_ycen->SetLineColor(kRed+2);
    flatv2_ycen->SetLineColor(kBlue+2);
    flatv2b_ycen->SetLineColor(kGreen+2);
    flatv2c_ycen->SetLineColor(kOrange+2);

    flat_xcen->SetMarkerStyle(20);
    flatv2_xcen->SetMarkerStyle(21);
    flatv2b_xcen->SetMarkerStyle(22);
    flatv2c_xcen->SetMarkerStyle(23);


    flat_xcen->SetMarkerColor(kRed-2);
    flatv2_xcen->SetMarkerColor(kBlue-2);
    flatv2b_xcen->SetMarkerColor(kGreen-2);
    flatv2c_xcen->SetMarkerColor(kOrange-2);


    flat_xcen->SetLineColor(kRed-2);
    flatv2_xcen->SetLineColor(kBlue-2);
    flatv2b_xcen->SetLineColor(kGreen-2);
    flatv2c_xcen->SetLineColor(kOrange-2);


    flat->SetMaximum(1.02);
    flat->SetMinimum(0.87);
    can->cd();
    flat->Draw("ALP");
    flatv2->Draw("SAMELP");
    flatv2b->Draw("SAMELP");
    flatv2c->Draw("SAMELP");


    flat_ycen->Draw("SAMELP");
    flatv2_ycen->Draw("SAMELP");
    flatv2b_ycen->Draw("SAMELP");
    flatv2c_ycen->Draw("SAMELP");

    flat_xcen->Draw("SAMELP");
    flatv2_xcen->Draw("SAMELP");
    flatv2b_xcen->Draw("SAMELP");
    flatv2c_xcen->Draw("SAMELP");


    auto legend = new TLegend(0.85, 0.55, 1. , 0.9);
    legend->AddEntry(flat, "flat-yoff", "lpfe");
    legend->AddEntry(flat_xcen, "flat-xcen", "lpfe");
    legend->AddEntry(flat_ycen, "flat-ycen", "lpfe");

    legend->AddEntry(flatv2, "flat-v2-yoff", "lpfe");
    legend->AddEntry(flatv2_xcen, "flat-v2.1-xcen", "lpfe");
    legend->AddEntry(flatv2_ycen, "flat-v2.1-ycen", "lpfe");

    legend->AddEntry(flatv2b, "flat-v2.1-b-yoff", "lpfe");
    legend->AddEntry(flatv2b_xcen, "flat-v2.1-b-xcen", "lpfe");
    legend->AddEntry(flatv2b_ycen, "flat-v2.1-b-ycen", "lpfe");

    legend->AddEntry(flatv2c, "flat-v2.1-c-yoff", "lpfe");
    legend->AddEntry(flatv2c_xcen, "flat-v2.1-c-xcen", "lpfe");
    legend->AddEntry(flatv2c_ycen, "flat-v2.1-c-ycen", "lpfe");

    legend->Draw();

    can->Draw();
    can->SaveAs("results/pde_flat_all.png");
}

void dcr(){
    auto can = new TCanvas("can", "DCR", 1200, 900);

    auto *flat = background("flat/flat.scanx.yoff", "x", 0);
    auto *flatv2 = background("flat-v2.1/flat-v2.1.scanx.yoff", "x", 0);
    auto *flatv2b = background("flat-v2.1-b/flat-v2.1-b.scanx.yoff", "x", 0);
    auto *flatv2c = background("flat-v2.1-c/flat-v2.1-c.scanx.yoff", "x", 0);


    auto *flat_ycen = background("flat/flat.scanx.ycen", "x", 1);
    auto *flatv2_ycen = background("flat-v2.1/flat-v2.1.scanx.ycen", "x", 1);
    auto *flatv2b_ycen = background("flat-v2.1-b/flat-v2.1-b.scanx.ycen","x", 1);
    auto *flatv2c_ycen = background("flat-v2.1-c/flat-v2.1-c.scanx.ycen","x", 1);


    auto *flat_xcen = background("flat/flat.scany.xcen", "y", 1);
    auto *flatv2_xcen = background("flat-v2.1/flat-v2.1.scany.xcen", "y", 1);
    auto *flatv2b_xcen = background("flat-v2.1-b/flat-v2.1-b.scany.xcen", "y", 1);
    auto *flatv2c_xcen = background("flat-v2.1-c/flat-v2.1-c.scany.xcen", "y", 1);

    flat->SetTitle(" ; Sensor ; DCR [Hz]"); // scanx and yoff data

    flat->SetMarkerStyle(20);
    flatv2->SetMarkerStyle(21);
    flatv2b->SetMarkerStyle(22);
    flatv2c->SetMarkerStyle(23);


    flat->SetMarkerColor(kRed);
    flatv2->SetMarkerColor(kBlue);
    flatv2b->SetMarkerColor(kGreen);
    flatv2c->SetMarkerColor(kOrange);


    flat->SetLineColor(kRed);
    flatv2->SetLineColor(kBlue);
    flatv2b->SetLineColor(kGreen);
    flatv2c->SetLineColor(kOrange);


    flat_ycen->SetMarkerStyle(20);
    flatv2_ycen->SetMarkerStyle(21);
    flatv2b_ycen->SetMarkerStyle(22);
    flatv2c_ycen->SetMarkerStyle(23);

    flat_ycen->SetMarkerColor(kRed+2);
    flatv2_ycen->SetMarkerColor(kBlue+2);
    flatv2b_ycen->SetMarkerColor(kGreen+2);
    flatv2c_ycen->SetMarkerColor(kOrange+2);


    flat_ycen->SetLineColor(kRed+2);
    flatv2_ycen->SetLineColor(kBlue+2);
    flatv2b_ycen->SetLineColor(kGreen+2);
    flatv2c_ycen->SetLineColor(kOrange+2);

    flat_xcen->SetMarkerStyle(20);
    flatv2_xcen->SetMarkerStyle(21);
    flatv2b_xcen->SetMarkerStyle(22);
    flatv2c_xcen->SetMarkerStyle(23);


    flat_xcen->SetMarkerColor(kRed-2);
    flatv2_xcen->SetMarkerColor(kBlue-2);
    flatv2b_xcen->SetMarkerColor(kGreen-2);
    flatv2c_xcen->SetMarkerColor(kOrange-2);


    flat_xcen->SetLineColor(kRed-2);
    flatv2_xcen->SetLineColor(kBlue-2);
    flatv2b_xcen->SetLineColor(kGreen-2);
    flatv2c_xcen->SetLineColor(kOrange-2);


    // flat->SetMaximum(2050);
    // flat->SetMinimum(1250);
    can->cd();
    // can->Range(-1, 1250, 5, 2050);
    flat->SetMaximum(2050);
    flat->SetMinimum(1250);
    flat->Draw("ALP");
    flatv2->Draw("SAMELP");
    flatv2b->Draw("SAMELP");
    flatv2c->Draw("SAMELP");


    flat_ycen->Draw("SAMELP");
    flatv2_ycen->Draw("SAMELP");
    flatv2b_ycen->Draw("SAMELP");
    flatv2c_ycen->Draw("SAMELP");

    flat_xcen->Draw("SAMELP");
    flatv2_xcen->Draw("SAMELP");
    flatv2b_xcen->Draw("SAMELP");
    flatv2c_xcen->Draw("SAMELP");


    auto legend = new TLegend(0.85, 0.55, 1. , 0.9);
    legend->AddEntry(flat, "flat-yoff", "lpfe");
    legend->AddEntry(flat_xcen, "flat-xcen", "lpfe");
    legend->AddEntry(flat_ycen, "flat-ycen", "lpfe");

    legend->AddEntry(flatv2, "flat-v2-yoff", "lpfe");
    legend->AddEntry(flatv2_xcen, "flat-v2.1-xcen", "lpfe");
    legend->AddEntry(flatv2_ycen, "flat-v2.1-ycen", "lpfe");

    legend->AddEntry(flatv2b, "flat-v2.1-b-yoff", "lpfe");
    legend->AddEntry(flatv2b_xcen, "flat-v2.1-b-xcen", "lpfe");
    legend->AddEntry(flatv2b_ycen, "flat-v2.1-b-ycen", "lpfe");

    legend->AddEntry(flatv2c, "flat-v2.1-c-yoff", "lpfe");
    legend->AddEntry(flatv2c_xcen, "flat-v2.1-c-xcen", "lpfe");
    legend->AddEntry(flatv2c_ycen, "flat-v2.1-c-ycen", "lpfe");

    legend->Draw();

    can->SetLeftMargin(0.1362725);
    can->SetBottomMargin(0.1316309);

    can->Draw();
    can->SaveAs("results/dcr_flat_all.png");
}
