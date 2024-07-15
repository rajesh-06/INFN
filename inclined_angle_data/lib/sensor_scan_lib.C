#pragma once

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <numeric>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <vector>
#include <cmath>   // For std::sqrt, std::acos, and std::abs
#include <TString.h>

using namespace std;

// #pragma once
// #include <iostream>
// #include <vector>
// #include "lib/sensor_scan_lib.C"
// #include <TString.h>
// #include <TMath.h>
// #include <numeric> // For std::accumulate

// using namespace std;

// Function to compute efficiency from sensor scan data


double 
propagate_error_division(double a, double a_err, double b, double b_err) {
    // Check for division by zero in a or b
    if (std::abs(a) < 1e-10 || std::abs(b) < 1e-10) {
        std::cerr << "Error: Division by zero in propagate_error_division." << std::endl;
        return 0.0; // Return 0.0 or handle the error as appropriate for your application
    }

    return std::sqrt((a_err / a) * (a_err / a) + (b_err / b) * (b_err / b)) * (a / b);
}

double 
propagate_error_multiplication(const vector<double>& a, const vector<double>& a_err) {
    if (a.size() != a_err.size()) {
        cout << "Sizes of data vector and error of data vector are not the same" << endl;
        return -1.0; // Indicate an error condition
    }

    double error = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != 0) {
            double relative_error = a_err[i] / a[i];
            error += relative_error * relative_error;
        } else {
            cout << "Division by zero encountered at index " << i << endl;
            return -1.0; // Indicate an error condition
        }
    }

    return std::sqrt(error);
}

double 
propagate_error_add_sub(const vector<double>& err_vec) {
    double error = 0.0;
    for (auto err : err_vec) {
        error += err * err;
    }
    return std::sqrt(error);
}

//To get the minimum Y value of TGraphErrors
Double_t GetYMinimum(TGraphErrors *gr) {
  Int_t n = gr->GetN();
  Double_t *y = gr->GetY();
  return TMath::MinElement(n, y);
}

//To get the maximum Y value of TGraphErrors
Double_t GetYMaximum(TGraphErrors *gr) {
  Int_t n = gr->GetN();
  Double_t *y = gr->GetY();
  return TMath::MaxElement(n, y);
}

//Defining the fit function: 
//Two opposite fermi function has been multiplied and gaussian is subtracted to get the dip at the center

double doubleFermi(double *x, double *par) {
    double xx = x[0];
    double fermi1 = par[1] / (TMath::Exp((xx - par[3]) / par[4]) + 1.0);
    double fermi2 = 1.0 / (TMath::Exp((par[2] - xx) / par[5]) + 1.0);

    //p0 is background
    //p1 is max rate without background
    //p2 is position of the one end of the sensor
    //p3 is position of the other end of the sensor

    return par[0] + (fermi1 * fermi2);// - gaus_dip;

}
double doubleFermiDip(double *x, double *par) {
    double xx = x[0];
    double fermi1 = par[1] / (TMath::Exp((xx - par[3]) / par[4]) + 1.0);
    double fermi2 = 1.0 / (TMath::Exp((par[2] - xx) / par[5]) + 1.0);
    double gaus_dip = par[6] * TMath::Exp(-0.5 * ((xx - par[7]) / par[8]) * ((xx - par[7]) / par[8]));

    //p0 is background
    //p1 is max rate without background
    //p2 is position of the one end of the sensor
    //p3 is position of the other end of the sensor
    //p6 is the dip height
    //p7 is the center of the dip (where anode is there in the sensor)

    return par[0] + (fermi1 * fermi2) - gaus_dip;
}


TF1* get_scan_fit(TGraphErrors *f1) {
    // vector<double> result;
    if (!f1) {
        std::cerr << "Error: TGraphErrors object f1 is null!" << std::endl;
        return nullptr;
    }

    //Defining fitting parameters 
    double xmin, xmax, bkg, max_rate;

    //Optimising the initial parameters 
    xmin = f1->GetPointX(0); //staring point of graph
    xmax = f1->GetPointX(f1->GetN()); //emd point of graph
    bkg = GetYMinimum(f1); //DCR at the begining of the scan taken as background
    max_rate = GetYMaximum(f1)-bkg; //rate of the laser
    double xlow = 1000.;
    double xhigh = -1000.;
    double threshold = (bkg+max_rate)*0.5;

    for(int i = 0; i<f1->GetN(); ++i){
       if (threshold < f1->GetPointY(i)){
            if(f1->GetPointX(i) <xlow){
                xlow = f1->GetPointX(i) ;
            }
            if(f1->GetPointX(i) >xhigh){
                xhigh=f1->GetPointX(i) ;
            }
        }
    }
    cout<<f1->GetTitle()<<endl;

    // Define the fit function
    TF1 *fermidip = new TF1("Double Fermi with Dip", doubleFermiDip, xlow-5, xhigh+5, 9); 
    fermidip->SetParameters(bkg, max_rate, xlow, xhigh, 0.05, 0.05, 5000, (xlow+xhigh)*0.5, 0.2); // Initial guess for the parameters
    fermidip->SetNpx(500);
    //f1->Fit(fermidip, "R"); // "R" stands for fit within the range specified in TF1
    TFitResultPtr fitResult = f1->Fit(fermidip, "R S");
    if (fitResult->Status() == 0) {
        std::cout << "Dip found." << std::endl;
        return fermidip;
    } 
    else{
        cout<<"dip not found"<<endl;
        TF1 *doublefermi = new TF1("Double Fermi without Dip", doubleFermi, xlow-5, xhigh+5, 6);
        doublefermi->SetParameters(bkg, max_rate, xlow, xhigh, 0.05, 0.05); // Initial guess for the parameters
        doublefermi->SetNpx(500);
        f1->Fit(doublefermi, "R"); // "R" stands for fit within the range specified in TF1
        cout<<"fitted with double fermi"<<endl;
        return doublefermi;
    }
    //cout<<fermidip->GetChisquare()<<endl;
    //cout<<fermidip->GetParameter(0)<<endl;

    // for(int i=0; i<4;++i){
    //     result.push_back(fermidip->GetParameter(i));
    // }
    // result.push_back(fermidip->GetParameter(7));

}


TF1* get_fermi_fit(TGraphErrors *f1) {
    // vector<double> result;
    if (!f1) {
        std::cerr << "Error: TGraphErrors object f1 is null!" << std::endl;
        return nullptr;
    }

    //Defining fitting parameters 
    double xmin, xmax, bkg, max_rate;

    //Optimising the initial parameters 
    xmin = f1->GetPointX(0); //staring point of graph
    xmax = f1->GetPointX(f1->GetN()); //emd point of graph
    bkg = GetYMinimum(f1); //DCR at the begining of the scan taken as background
    max_rate = GetYMaximum(f1)-bkg; //rate of the laser
    double xlow = 1000.;
    double xhigh = -1000.;
    double threshold = (bkg+max_rate)*0.5;

    for(int i = 0; i<f1->GetN(); ++i){
       if (threshold < f1->GetPointY(i)){
            if(f1->GetPointX(i) <xlow){
                xlow = f1->GetPointX(i) ;
            }
            if(f1->GetPointX(i) >xhigh){
                xhigh=f1->GetPointX(i) ;
            }
        }
    }
    cout<<f1->GetTitle()<<endl;

    TF1 *doublefermi = new TF1("Double Fermi without Dip", doubleFermi, xmin, xmax, 6);
    doublefermi->SetParameters(bkg, max_rate, xlow, xhigh, 0.05, 0.05); // Initial guess for the parameters
    doublefermi->SetNpx(500);
    f1->Fit(doublefermi, "R"); // "R" stands for fit within the range specified in TF1
    return doublefermi;

}

void FitScan(TGraphErrors *f1) {
    if (!f1) {
        std::cerr << "Error: TGraphErrors object f1 is null!" << std::endl;
        return;
    }


    float xmin, xmax, bkg, max_rate;
    
    xmin = f1->GetPointX(0);
    xmax = f1->GetPointX(f1->GetN());
    bkg = GetYMinimum(f1);
    max_rate = GetYMaximum(f1)-bkg;
    float xlow = 1000.;
    float xhigh = -1000.;
    float threshold = (bkg+max_rate)*0.5;

    for(int i = 0; i<f1->GetN(); ++i){
       if (threshold < f1->GetPointY(i)){
            if(f1->GetPointX(i) <xlow){
                xlow = f1->GetPointX(i) ;
            }
            if(f1->GetPointX(i) >xhigh){
                xhigh=f1->GetPointX(i) ;
            }
        }
    }

    // Define the function to fit
    cout<<f1->GetTitle()<<endl;
    TF1 *fermidip = new TF1("Double Fermi with Dip", doubleFermiDip, xlow-5, xhigh+5, 9); // Adjust the range [-10, 10] as needed
    fermidip->SetParameters(bkg, max_rate, xlow, xhigh, 0.05, 0.05, 5000, (xlow+xhigh)*0.5, 0.2); // Initial guess for the parameters

    // Fit the TGraphErrors
    cout<<f1->GetTitle()<<endl;
    f1->Fit(fermidip, "R"); // "R" stands for fit within the range specified in TF1
    cout<<fermidip->GetChisquare()<<endl;
    cout<<fermidip->GetNDF()<<endl;
    // c1->SaveAs("fit_result.png");
}



TGraphErrors *get_sensor_scan(std::string infilename = "data/scanx.45deg.target.A4.txt.tree.root",  std::string scan_axis = "x") {
    
    TFile *input_file = TFile::Open(infilename.c_str());
    //auto *input_file = new TFile(infilename.c_str());// Open file
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return nullptr;
    }
    auto data_tree = (TTree *)(input_file->Get("tree"));//Get tree

    // Create data structure to read tree and set addresses
    float x, y, rate; // data are in float
    if( scan_axis == "y"){
        data_tree->SetBranchAddress("x", &y);
        data_tree->SetBranchAddress("y", &x);
        data_tree->SetBranchAddress("rate", &rate);
    }
    else{
        data_tree->SetBranchAddress("x", &x);
        data_tree->SetBranchAddress("y", &y);
        data_tree->SetBranchAddress("rate", &rate);
    }
    
    // Vector to store data
    std::vector<float> vx;
    std::vector<float> vrate;

    // Read entries from the TTree
    std::cout << "Entries in tree: " << data_tree->GetEntries() << std::endl;
    for (int i = 0; i < data_tree->GetEntries(); i++) {
        data_tree->GetEntry(i);

        // Filter out rates less than 1.0
        if (rate < 1.0)
            continue;

        vx.push_back(x);
        vrate.push_back(rate);
    }

    // Map to store grouped vrate values by vx
    std::map<float, std::vector<float>> grouped_values;

    // Group vrate values by vx
    for (size_t i = 0; i < vx.size(); ++i) {
        grouped_values[vx[i]].push_back(vrate[i]);
    }

    // Create TGraphErrors to store the results
    TGraphErrors *rate_scan = new TGraphErrors();
    int index = 0;
    for (auto it = grouped_values.begin(); it != grouped_values.end(); ++it) {
        float avg_x = it->first;
        std::vector<float> vvrate = it->second;

        // Calculate average rate and error in average rate
        float avg_rate = std::accumulate(vvrate.begin(), vvrate.end(), 0.0) / vvrate.size();
        float std_rate = TMath::StdDev(vvrate.begin(), vvrate.end());
        float err_rate = std_rate / TMath::Sqrt(vvrate.size());

        rate_scan->SetPoint(index, avg_x, avg_rate);  //Saving the data to the graph
        rate_scan->SetPointError(index, 0, err_rate); // Setting the errors

        index++;
    }

    std::string axislabel = "; Position [mm]; Average Rate [Hz]";

    rate_scan->SetTitle((infilename + axislabel).c_str());
    rate_scan->GetYaxis()->SetTitleOffset(1.5);
    rate_scan->GetYaxis()->SetLabelOffset(0.01);
    rate_scan->SetMarkerStyle(20);
    rate_scan->SetMarkerSize(0.3);
    
    input_file->Close();
    delete input_file;
    return rate_scan;
}
vector<double> PDEvsAngle(string fileinitial = "flat/flat.scanx.yoff", string scan_axis = "x") {
    vector<string> filenames = {
        fileinitial + ".chip-0.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A1.txt.tree.root",
        fileinitial + ".chip-1.channel-A2.txt.tree.root",
        fileinitial + ".chip-1.channel-A3.txt.tree.root",
        fileinitial + ".chip-1.channel-A4.txt.tree.root"
    };

    vector<double> rates;
    vector<double> err_rates;
    vector<double> target_sensor_width, err_target_sensor_width;

    auto ref_gr = get_sensor_scan(filenames[0], scan_axis);
    if (!ref_gr) {
        cerr << "Error: Could not get reference graph from file " << filenames[0] << endl;
        return {};
    }
    auto ref_res = get_scan_fit(ref_gr);

    double ref_rate = ref_res->GetParameter(1);
    double err_ref_rate = ref_res->GetParError(1);
    double ref_sensor_width = ref_res->GetParameter(3) - ref_res->GetParameter(2);
    vector<double> temp_err = {static_cast<double>(ref_res->GetParError(3)), static_cast<double>(ref_res->GetParError(2))};
    double err_ref_sensor_width = propagate_error_add_sub(temp_err);

    for (size_t i = 1; i < filenames.size(); ++i) {
        auto graph = get_sensor_scan(filenames[i], scan_axis);
        if (!graph) {
            cerr << "Error: Could not get graph from file " << filenames[i] << endl;
            continue;
        }
        auto result = get_scan_fit(graph);
        rates.push_back(result->GetParameter(1));
        err_rates.push_back(result->GetParError(1));

        target_sensor_width.push_back(result->GetParameter(3) - result->GetParameter(2));
        vector<double> temp_err1 = {static_cast<double>(result->GetParError(3)), static_cast<double>(result->GetParError(2))};
        err_target_sensor_width.push_back(propagate_error_add_sub(temp_err1));
    }

    double avg_rate = accumulate(rates.begin(), rates.end(), 0.0) / rates.size();
    double err_avg_rate = propagate_error_add_sub(err_rates) / err_rates.size();

    double err_result = propagate_error_division(avg_rate, err_avg_rate, ref_rate, err_ref_rate);

    double avg_width = accumulate(target_sensor_width.begin(), target_sensor_width.end(), 0.0) / target_sensor_width.size();
    double err_width = propagate_error_add_sub(err_target_sensor_width) / err_target_sensor_width.size();

    double angle_rad = std::acos(avg_width / ref_sensor_width);
    double angle_deg = angle_rad * 180.0 / TMath::Pi();

    // Calculate the error in the angle
    double err_angle_rad = std::abs(1 / std::sqrt(1 - std::pow(avg_width / ref_sensor_width, 2))) * propagate_error_division(avg_width, err_width, ref_sensor_width, err_ref_sensor_width);
    double err_angle_deg = err_angle_rad * 180.0 / TMath::Pi();

    // Create and return result vector
    vector<double> result = {angle_deg, avg_rate / ref_rate, err_result, err_angle_deg};
    return result;
}

