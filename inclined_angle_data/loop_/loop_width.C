#pragma once
#include "header.C"


using namespace std;



void loop_rate() {
    TCanvas *c1 = new TCanvas("c1", "Sensor width analysis", 1200, 800);


    std::vector<std::string> filenames;
    std::vector<std::string> file_ycen;
    std::vector<std::string> file_xcen;

    std::vector<string> measurements;
    vector <float> x;

    // Loop from 0 to 33 (inclusive)
    for (int i = 0; i <= 99; ++i) {
        auto yoff = ".scanx.yoff.chip-0.channel-A1.txt.tree.root";
        auto ycen = ".scanx.ycen.chip-0.channel-A1.txt.tree.root";
        auto xcen = ".scany.xcen.chip-0.channel-A1.txt.tree.root";
        
        std::stringstream filename_stream;

        // Construct the filename with padding
        filename_stream << "loop-";
        filename_stream.width(2);
        filename_stream.fill('0');
        filename_stream << i;

    // Add the filename to the vector
    measurements.push_back(filename_stream.str());
    filenames.push_back(filename_stream.str()+"/"+filename_stream.str()+yoff);
    file_ycen.push_back(filename_stream.str()+"/"+filename_stream.str()+ycen);
    file_xcen.push_back(filename_stream.str()+"/"+filename_stream.str()+xcen);

    cout<<filenames[i]<<endl;
    cout<<measurements[i]<<endl;
    x.push_back(i+1);
    
    }

    

    auto *gr_yoffl = get_graph(filenames, 2);
    auto *gr_yoffr = get_graph(filenames, 3);

    auto *gr_ycenl = get_cen_graph(file_ycen, 2);
    auto *gr_ycenr = get_cen_graph(file_ycen, 3);

    auto *gr_xcenl = get_cen_graph(file_xcen, 2, "y");
    auto *gr_xcenr = get_cen_graph(file_xcen, 3, "y");

    /*
    auto _ymin = {GetYMinimum(gr), GetYMinimum(gr_ycen), GetYMinimum(gr_xcen)};
    auto _ymax = {GetYMaximum(gr), GetYMaximum(gr_ycen), GetYMaximum(gr_xcen)};
    
    
    auto ymax = std::max_element(_ymax.begin(), _ymax.end());
    if (ymax != _ymax.end()) {
        std::cout << "The maximum element is: " << *ymax << std::endl;
    }

    // Find the minimum element
    auto ymin = std::min_element(_ymin.begin(), _ymin.end());
    if (ymin != _ymin.end()) {
        std::cout << "The minimum element is: " << *ymin << std::endl;
    }
    Double_t ymin_ = *ymin;
    Double_t ymax_ = *ymax;
    */

    gr->SetMinimum(ymin_ - 0.1 * (ymax_ - ymin_)); // Adjust buffer
    gr->SetMaximum(ymax_ + 0.1 * (ymax_ - ymin_));

    // auto *gr_high = get_graph(filenames, 3);

    // auto gr = new TGraphErrors(gr_low->GetN());
    // for (int i = 0; i < gr_low->GetN(); i++){
    //     gr->SetPoint(i, i+1, gr_high->GetPointY(i) - gr_low->GetPointY(i) );
    //     auto error = propagate_error_add_sub({gr_high->GetErrorY(i), gr_low->GetErrorY(i) });
    //     gr->SetPointError(i, 0, error );

    // }

    //auto gr = PlotPulls(gr_a1);
  //string _data = " Photon rate detected by the sensor";
  string _data = "sensor width";

  string title = "; loop number ;"+ _data+" [mm]";


  gr->SetTitle(title.c_str());
  gr->Draw("ALP");
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1);
  gr->SetMarkerColor(kRed);

  gr_ycen->Draw("samelp");
  gr_ycen->SetMarkerStyle(21);
  gr_ycen->SetMarkerSize(1);
  gr_ycen->SetMarkerColor(kBlue);

  gr_xcen->Draw("samelp");
  gr_xcen->SetMarkerStyle(21);
  gr_xcen->SetMarkerSize(1);
  gr_xcen->SetMarkerColor(kViolet);




  
    auto *legend = new TLegend(0.7,0.75, 0.9, 0.87);
    legend->AddEntry(gr, "yoff", "lpfe");
    legend->AddEntry(gr_ycen, "ycen", "lpfe");
    legend->AddEntry(gr_xcen, "xcen", "lpfe");


    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->Draw();

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
    string outfile = "results/loop_sensor_width.png";
    c1->SaveAs(outfile.c_str());

}
