// #include "efficiency.C"
#include "lib/sensor_scan_lib.C"
// void test1(){
//     auto *c1 = new TCanvas("chip0", "both fit", 1200, 600);
//     c1->Divide(2,1);
//     c1->cd(1);
//     auto ch1 = get_sensor_scan("flat.scanx.yoff.chip-0.channel-A1.txt.tree.root");
//     ch1->Draw("ap");
//     auto fit1 = get_fit_result(ch1);
//     TLine *l1 = new TLine(-93,fit1->GetParameter(1)+fit1->GetParameter(0) , -87, fit1->GetParameter(1)+fit1->GetParameter(0));
//     l1->Draw("same");
//     c1->Update();
//     c1->cd(2);
//     auto ch2 = get_sensor_scan("flat.scanx.ycen.chip-0.channel-A1.txt.tree.root");
//     ch2->Draw("ap");
//     auto fit2 = get_fit_result(ch2);
//     TLine *l2 = new TLine(-93,fit2->GetParameter(1)+fit2->GetParameter(0) , -87, fit2->GetParameter(1)+fit2->GetParameter(0));
//     l2->Draw("same");
//     c1->Update();




// }
void PDE() {
    vector <string> directory = {
        "flat/flat.scanx.yoff",
        "30deg/30deg.scanx.yoff",
        "40deg/40deg.scanx.yoff",
        "40deg-a/40deg-a.scanx.yoff",
        "40deg-b/40deg-b.scanx.yoff",
        "auto_scan_50deg/test-autoscan.scanx.yoff"
    };

    TCanvas *c1 = new TCanvas("c1", "Relative PDE", 800, 600);

    auto *gr_pde = new TGraphErrors();

    
    auto ref_pde = PDEvsAngle(directory[0]);
    // auto f30_pde = PDEvsAngle(directory[1]);
    // auto f40_pde = PDEvsAngle(directory[2]);
    // auto f50_pde = PDEvsAngle(directory[3]);

    // cout<<"\nref"<<endl;
    // for (auto val: ref_pde){
    //      cout<<val<<"\t";
    // }

    // cout<<"\n30 deg"<<endl;
    // for (auto val: f30_pde){
    //      cout<<val<<"\t";
    // }

    // cout<<"\n40 deg"<<endl;
    // for (auto val: f40_pde){
    //      cout<<val<<"\t";
    // }

    // cout<<"\n50 deg"<<endl;
    // for (auto val: f50_pde){
    //      cout<<val<<"\t";
    // }
    // cout<<endl;



    for (int i = 0; i < directory.size(); ++i) {
        auto target_pde = PDEvsAngle(directory[i]);
        gr_pde->SetPoint(i, target_pde[0], target_pde[2]/ref_pde[2]);
        gr_pde->SetPointError(i, target_pde[1], propagate_error_division(target_pde[2], target_pde[3], ref_pde[2], ref_pde[3]));
        // gr_pde->SetPointError(i, target_pde[], target_pde[3]);

    }


    gr_pde->SetTitle(";Incidence Angle (deg);Relative PDE");
    gr_pde->SetMarkerStyle(21); // Marker style
    gr_pde->SetMarkerSize(1); // Marker size
    gr_pde->SetMarkerColor(kViolet); // Marker color
    //gr_pde->SetLineColor(kCyan); // Line color
    gr_pde->SetLineWidth(2); // Line width

    c1->Range(-5, 0, 60, 1.2);
    c1->SetLeftMargin(0.1362725);
    c1->SetBottomMargin(0.1316309);

    

    gr_pde->GetXaxis()->SetTitleSize(0.04); // X axis title size
    gr_pde->GetYaxis()->SetTitleSize(0.04); // Y axis title size
    gr_pde->GetXaxis()->SetLabelSize(0.04); // X axis label size
    gr_pde->GetYaxis()->SetLabelSize(0.04); // Y axis label size
    gr_pde->GetXaxis()->SetTitleOffset(1.2); // X axis title offset
    gr_pde->GetYaxis()->SetTitleOffset(1.5); // Y axis title offset

    gr_pde->Draw("AP");
    c1->Draw();
    c1->SaveAs("results/relative_photodetction_efficiency.png");

    // delete c1;
    // delete gr_pde;
}


//     // auto pde_40a = PDEvsAngle("40deg-a/40deg-a.scanx.yoff");
//     // auto pde_40b = PDEvsAngle("40deg-b/40deg-b.scanx.yoff");

//     // auto *highlight = new TGraphErrors();
//     // highlight->SetPoint(0, pde_40a[0], pde_40a[2] / ref_pde[2]);
//     // highlight->SetPointError(0, pde_40a[1], propagate_error_division(pde_40a[2], pde_40a[3], ref_pde[2], ref_pde[3]));

//     // highlight->SetPoint(1, pde_40b[0], pde_40b[2] / ref_pde[2]);
//     // highlight->SetPointError(1, pde_40b[1], propagate_error_division(pde_40b[2], pde_40b[3], ref_pde[2], ref_pde[3]));

//     // highlight->SetMarkerColor(kRed); // Change color to red
//     // highlight->SetMarkerStyle(21);
//     // highlight->SetMarkerSize(1.5); // Marker size

//     // auto _40val = gr_pde->GetPointY(2);
    
//     // // Label the coordinates
//     // auto *text = new TLatex(35, _40val - 0.0005, Form("40(%d, %.3f)", 46, _40val));
//     // text->SetTextColor(kMagenta);
//     // text->SetTextSize(0.03);
    
//     // auto *text1 = new TLatex(pde_40a[0], pde_40a[2] / ref_pde[2] - 0.0005, Form("40a(%d, %.3f)", pde_40a[0], pde_40a[2] / ref_pde[2]));
//     // text1->SetTextColor(kRed);
//     // text1->SetTextSize(0.03);

//     // auto *textb = new TLatex(pde_40b[0], pde_40b[2] / ref_pde[2] - 0.0005, Form("40b(%d, %.3f)", pde_40b[0], pde_40b[2] / ref_pde[2]));
//     // textb->SetTextColor(kRed);
//     // textb->SetTextSize(0.03);
// //      highlight->Draw("SAMEP");
// //     text->Draw("SAME");
// //     text1->Draw("SAME");
// //     textb->Draw("SAME");