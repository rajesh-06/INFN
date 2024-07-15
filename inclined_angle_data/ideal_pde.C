#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <lib/sensor_scan_lib.C>

// Function to compute efficiency from sensor scan data
vector<double> 
PDEvsAngle(string fileinitial = "flat/flat.scanx.yoff", string scan_axis = "x") {
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
    auto ref_res = get_fermi_fit(ref_gr);

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
        auto result = get_fermi_fit(graph);
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
    
    if (avg_width / ref_sensor_width > 1.0){
        double angle_rad = std::acos(ref_sensor_width / avg_width);
        double angle_deg = angle_rad * 180.0 / TMath::Pi();

        // Calculate the error in the angle
        double err_angle_rad = std::abs(1 / std::sqrt(1 - std::pow(ref_sensor_width / avg_width, 2))) * propagate_error_division(ref_sensor_width, err_ref_sensor_width, avg_width, err_width);
        double err_angle_deg = err_angle_rad * 180.0 / TMath::Pi();
        return {angle_deg, err_angle_deg, avg_rate / ref_rate, err_result};
    }
    // else{
    double angle_rad = std::acos(avg_width / ref_sensor_width);
    double angle_deg = angle_rad * 180.0 / TMath::Pi();

    // Calculate the error in the angle
    double err_angle_rad = std::abs(1 / std::sqrt(1 - std::pow(avg_width / ref_sensor_width, 2))) * propagate_error_division(avg_width, err_width, ref_sensor_width, err_ref_sensor_width);
    double err_angle_deg = err_angle_rad * 180.0 / TMath::Pi();
    return {angle_deg, err_angle_deg, avg_rate / ref_rate, err_result};


    // }
    
}

Double_t func(double *theta, double *par) {
    double n = par[0];
    double C = par[1];

    Double_t dell_l = 1.0 / std::sqrt(1.0 - 1.0/(n * n) * std::sin(theta[0]) * std::sin(theta[0])) - 1.0;
    return std::exp(-C * dell_l);
}
void ideal_pde() {
    auto f1 = new TF1("f1", func, 0, TMath::Pi() / 2, 2); // note the 2 at the end for number of parameters
    f1->SetParameter(0, 1.1);
    f1->SetParameter(1, 1.0);

    // Set axis labels
    f1->GetXaxis()->SetTitle("#theta"); // Replace with your x-axis label
    f1->GetYaxis()->SetTitle("f(#theta)"); // Replace with your y-axis label

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    f1->Draw("AL");

    c1->Draw();
    // c1->SaveAs("function_plot.png"); // Saves the plot as a PNG file
}
