#pragma once
#include "lib/sensor_scan_lib.C"
// #include "findcenter_dip.C"
#include <TH1F.h>
#include <filesystem> // Use <experimental/filesystem> if necessary
#include <TString.h>  // ROOT's TString class



namespace fs = std::filesystem;

std::vector<string> getFilesInDirectory(const std::string& directoryPath) {
    std::vector<string> filePaths;
    
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            filePaths.push_back(entry.path().string().c_str());
        }
    }

    return filePaths;
}
void fit(){
    std::string folder = "40deg-b";
    // std::string out_result = folder+"_result.txt";


    std::vector<string> filenames = getFilesInDirectory(folder);

    std::ofstream outFile("results/"+folder+"_result.txt");
    if (!outFile) {
        std::cerr << "Error: Unable to open file " << folder+"_result.txt" << std::endl;
        return 1;
    }
    outFile << "File\t backgroud \t max rate \t xlow \t xhigh \t center \t dip center \t dip height"<<endl; 
    

    for (auto file: filenames){cout<<file<<endl;}
    // if(false){
    TCanvas *c1 = new TCanvas("c1", "Multiple Plots", 2400, 1600);
    c1->Divide(5, 2); // 3 columns and 5 rows

    int run = 0;
    for(auto file: filenames){
        outFile << file;
        if (run<10){
            auto *gr = get_sensor_scan(file);
            if(run<5){ 
                c1->cd(run+1);
                gr->Draw("ap"); 
            }
            else{
                c1->cd(run+1-5);
                gr->Draw("samep");
            }
            
            auto *fit = get_scan_fit(gr);
            outFile << "\t" << fit->GetParameter(0);
            outFile << "\t" << fit->GetParameter(1);
            outFile << "\t" << fit->GetParameter(2);
            outFile << "\t" << fit->GetParameter(3);
            outFile << "\t" << (fit->GetParameter(2)+fit->GetParameter(3))*0.5;
            outFile << "\t" << fit->GetParameter(7);
            outFile << "\t" << fit->GetParameter(6);

           // run++;
        }
        else{
            auto *gr = get_sensor_scan(file, "y");
            c1->cd(run+1-5);
            gr->Draw("ap");
            auto *fit = get_scan_fit(gr);
            outFile << "\t" << fit->GetParameter(0);//bkg
            outFile << "\t" << fit->GetParameter(1);//max rate - bkg
            outFile << "\t" << fit->GetParameter(2);//xlow
            outFile << "\t" << fit->GetParameter(3);//xhigh
            outFile << "\t" << (fit->GetParameter(2)+fit->GetParameter(3))*0.5;//center
            outFile << "\t" << fit->GetParameter(7);//dip center
            outFile << "\t" << fit->GetParameter(6);//dip height

            
        }
        run++;
        outFile<< endl;
        
    }
    // }
    std::string out_graph = "results/"+folder+"_scan.png";
    cout<<out_graph<<endl;
    c1->SaveAs(out_graph.c_str());
    

}
/*
void fit2(){
    std::string folder = "flat";

    std::vector<string> filenames = getFilesInDirectory(folder);

    
    std::string result_data = "result.txt";

    // Create an output file stream
    std::ofstream outFile(result_data);

    // Check if the file is opened successfully
    if (!outFile) {
        std::cerr << "Error: Unable to open file " << result_data << std::endl;
        return 1;
    }
    outFile << "File\t backgroud \t max rate \t xlow \t xhigh \t center \t dip center"<<endl; 
    int run = 0;
    for(auto files : filenames){
        
        if(run < 10){
            // c1->cd(run+1);
            auto gr = get_sensor_scan(files);
            vector<double> results = get_fit_result(gr);
            // run++;
            outFile << files ;
            cout<< results.size()<<endl;
            for (auto result : results){
                outFile<<"\t"<<result;
            }
        outFile<<"\t"<<(results[2]+results[3])*0.5<<endl;
        }
        else{
            // c1->cd(run+1);

            auto gr = get_sensor_scan(files, "y");
            auto results = get_fit_result(gr);
            outFile << files ;
            cout<< results.size()<<endl;
            for (auto result : results){
                outFile<<"\t"<<result;
            }
            // outFile<<endl;
            outFile<<"\t"<<(results[2]+results[3])*0.5<<endl;

        }
        // outFile << files ;
        // cout<< results.size()<<endl;
        // for (auto result : results){
        //     outFile<<"\t"<<result;
        // }
        // outFile<<endl;
        run++;
    }
}
*/