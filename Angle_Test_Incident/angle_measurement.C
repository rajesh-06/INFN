#include "iostream"
#include "fstream"

using namespace std;


//gROOT->gErrorIgnoreLevel = kFatal;// kInfo, kWarning, kError, kBreak, kSysError, kFatal;
const char* image_loc = "images2"; // dir name for output image
const char* input_loc = "input2"; // dir name for input data
const char* output_loc = "output2"; // dirname for output data TGraph
float angle_to_theta  = 180.0/TMath::Pi();// for conversion to angle 
int usefit = 1; // 1 if double fermi , 2 if double erf , 3 if double fermi with dip
float has_flat=false;

float error_in_x = 0.05; // increment of the x position
bool debug_fit = true;

//vector < const char *>  ang_name = {"flat", "15deg", "30deg", "40deg","45deg", "50deg"};
//vector < const char *>  ang_name = {"flat", "30deg", "40deg", "50deg","45deg", "50deg"};
vector < const char *>  ang_name = {"40deg", "40deg-a","40deg-b", "40deg-c", "40deg-d"};
int numangle = ang_name.size();//ang_name.size();// setting number of angle measurement currently we have only 2
    
   // vector < const char *>  ang_name = {"flat", "45deg", "40deg", "15deg","30deg"};
    

void setProfessionalPalette() {
    // Create a TStyle object
    TStyle* myStyle = new TStyle("myStyle", "My ROOT Style");

    // Use color wheel to define colors
    Int_t colors[8] = { TColor::GetColor("#4E79A7"),   // muted blue
                        TColor::GetColor("#F28E2B"),   // safety orange
                        TColor::GetColor("#E15759"),   // brick red
                        TColor::GetColor("#76B7B2"),   // muted cyan
                        TColor::GetColor("#59A14F"),   // muted green
                        TColor::GetColor("#EDC948"),   // muted yellow
                        TColor::GetColor("#B07AA1"),   // muted purple
                        TColor::GetColor("#FF9DA7") }; // muted pink

    myStyle->SetPalette(8, colors); // Set the palette with 8 colors

    // Set the style
    gStyle->SetOptStat(0); // Turn off statistics box
    gStyle->SetTitleAlign(23); // Center-align title
    //gStyle->SetTitleX(0.5); // Center title horizontally
    //gStyle->SetTitleY(0.98); // Adjust title vertically
    gStyle->SetTitleBorderSize(0); // Remove border around title

    gStyle->SetPadTickX(1); // Ticks on the x-axis of the pad
    gStyle->SetPadTickY(1); // Ticks on the y-axis of the pad

    // Apply the style
    gROOT->SetStyle("myStyle");
    //gStyle->SetPalette(kSolar);
    gROOT->ForceStyle(); // Force the style to take effect
}









vector <float> Fitter(TGraphErrors * _graph, float minvalx =-104., float maxvalx =-101. , float bkg = 1000. , float sigtop = 25000){
  if(debug_fit)  cout<< "###  RUNNING THE FITTER ####" << endl; 
    
    //vector to output
    vector <float> outvect;
    float mid = 0.5*(minvalx + maxvalx );
    float xlimhistmin = -200;
    float xlimhistmax = -70;
        
    
    // Multiple Functions to check the fit quality but only double fermi with dip is used 
    //double fermi function
    TF1 *doublefermiFunction = new TF1("doublefermiFunction", [](double *x, double *p) {
            return p[0]+ (p[1]/(TMath::Exp((x[0] - p[3]) /p[4]) + 1.0)/(TMath::Exp((p[2]-x[0] ) /p[5]) + 1.0));
    }, xlimhistmin, xlimhistmax, 6); // Defining the function from 0 to 10, 0 parameters

    doublefermiFunction->SetParameters(bkg,sigtop,minvalx, maxvalx,0.05,0.05); //setting parameters


    //double fermi with dip 
    TF1 *doublefermiFunction_dip = new TF1("doublefermiFunction_dip", [](double *x, double *p) {
                return p[0]+ (p[1]/(TMath::Exp((x[0] - p[3]) /p[4]) + 1.0)/(TMath::Exp((p[2]-x[0] ) /p[5]) + 1.0))-p[6]*exp(-0.5*((x[0]-p[7])/p[8])*((x[0]-p[7])/p[8]));
        }, xlimhistmin, xlimhistmax, 9); // Defining the function from 0 to 10, 0 parameters

        
    // second double fermi with dip  
    TF1 *doublefermiFunction_dip2 = new TF1("doublefermiFunction_dip2", [](double *x, double *p) {
                return p[0]+ (p[1]/(TMath::Exp((x[0] - p[3]) /p[4]) + 1.0)/(TMath::Exp((p[2]-x[0] ) /p[5]) + 1.0))-p[6]*exp(-0.5*((x[0]-p[7])/p[8])*((x[0]-p[7])/p[8]));
        }, xlimhistmin, xlimhistmax, 9); // Defining the function from 0 to 10, 0 parameters

     doublefermiFunction_dip2->SetNpx(4000);   



        //gaussian error function    
    TF1 * double_erf = new TF1("double_erf", [](double *x, double *p) {
            return p[0] + p[1] * (TMath::Erf((x[0] - p[2]) /p[4]) + 1.) * (TMath::Erf((p[3]-x[0] ) /p[5]) + 1.) * 0.25;
    },xlimhistmin, xlimhistmax, 6); // Defining the function from 0 to 10, 0 parameters


    double_erf->SetNpx(1000);
    doublefermiFunction->SetNpx(1000); 
    doublefermiFunction_dip->SetNpx(1000); 
    doublefermiFunction_dip->SetLineColor(kGreen-2);

    double_erf->SetParameters(bkg,sigtop,minvalx,maxvalx,0.1,0.1);
    double_erf->SetLineColor(kGreen+3);
    const char *op;
    if(debug_fit) {
        op = "NR";
    }
    else op = "NQR";
    
    
    if (debug_fit) cout<< "DOUBLE FERMI FUNCTION "<< endl;
    _graph->Fit("doublefermiFunction",op);
    
    
    if (debug_fit) cout<< "DOUBLE ERF FUNCTION "<< endl;
     _graph->Fit("double_erf", op);

 //    TF1 *gaus_dip = new TF1("gaus_dip", "[0]-[1]*exp(-0.5*((x[0]-[2])/[3])*((x[0]-[2])/[3]))",minvalx-2,maxvalx+2);
    
    
    
 /*   if(debug_fit) cout<< "Gaus DIP "<< endl;
    gaus_dip->SetParameters(sigtop,sigtop*0.1, mid, 0.3);    
    _graph->Fit("gaus_dip",op);
   */
    
    doublefermiFunction_dip->SetParameters(doublefermiFunction->GetParameter(0),doublefermiFunction->GetParameter(1),doublefermiFunction->GetParameter(2), doublefermiFunction->GetParameter(3),doublefermiFunction->GetParameter(4),doublefermiFunction->GetParameter(5),sigtop,mid,0.1);
    
    
    /* doublefermiFunction_dip->FixParameter(0,doublefermiFunction->GetParameter(0));
     doublefermiFunction_dip->FixParameter(1,doublefermiFunction->GetParameter(1));
     doublefermiFunction_dip->FixParameter(2,doublefermiFunction->GetParameter(2));
     doublefermiFunction_dip->FixParameter(3,doublefermiFunction->GetParameter(3));
     doublefermiFunction_dip->FixParameter(4,doublefermiFunction->GetParameter(4));
     //doublefermiFunction_dip->FixParameter(7,gaus_dip->GetParameter(2));
     //doublefermiFunction_dip->FixParameter(8,gaus_dip->GetParameter(3));*/
    
    
    //doublefermiFunction_dip->SetParLimits(6,0,10000);
    
    
    if (debug_fit) cout<< "DOUBLE FERMI FUNCTION  WITH DIP "<< endl;
    _graph->Fit("doublefermiFunction_dip",op);
    
    
    doublefermiFunction_dip2->SetParameters(doublefermiFunction_dip->GetParameter(0),doublefermiFunction_dip->GetParameter(1),doublefermiFunction_dip->GetParameter(2), doublefermiFunction_dip->GetParameter(3),doublefermiFunction_dip->GetParameter(4),doublefermiFunction_dip->GetParameter(5),doublefermiFunction_dip->GetParameter(6),doublefermiFunction_dip->GetParameter(7),doublefermiFunction_dip->GetParameter(8));


   //if(debug_fit)  cout<< "DOUBLE FERMI FUNCTION  WITH DIP 2 "<< endl;
   // _graph->Fit("doublefermiFunction_dip2",op);
   
    // doublefermiFunction_dip2->SetLineColor(kRed-2);
    
    
    if(usefit ==1){
        outvect.push_back(doublefermiFunction->GetParameter(3));
        outvect.push_back(doublefermiFunction->GetParameter(2));
        outvect.push_back(doublefermiFunction->GetParError(3));
        outvect.push_back(doublefermiFunction->GetParError(2));
        outvect.push_back(doublefermiFunction->GetParameter(1));
        outvect.push_back(doublefermiFunction->GetParError(1));

    }
     
    else if (usefit==2){
    
        outvect.push_back(doublefermiFunction_dip->GetParameter(3));
        outvect.push_back(doublefermiFunction_dip->GetParameter(2));
        outvect.push_back(doublefermiFunction_dip->GetParError(3));
        outvect.push_back(doublefermiFunction_dip->GetParError(2));
        outvect.push_back(doublefermiFunction_dip->GetParameter(1));
        outvect.push_back(doublefermiFunction_dip->GetParError(1));
    
    
    }
    
   else if (usefit ==3) {
        outvect.push_back(doublefermiFunction_dip->GetParameter(3));
        outvect.push_back(doublefermiFunction_dip->GetParameter(2));
        outvect.push_back(doublefermiFunction_dip->GetParError(3));
        outvect.push_back(doublefermiFunction_dip->GetParError(2));
        outvect.push_back(doublefermiFunction_dip->GetParameter(1));
        outvect.push_back(doublefermiFunction_dip->GetParError(1));
    
   
   
       
   
   }
    
    
   TLine  line(minvalx,0,minvalx,sigtop);
   TLine  line2(maxvalx,0,maxvalx,sigtop); 
   TLine  line3(doublefermiFunction_dip->GetParameter(3),0,doublefermiFunction_dip->GetParameter(3),sigtop ) ;
    
    
  doublefermiFunction->Draw("SAME");
    
   
  //double_erf->Draw("SAME");  
  //doublefermiFunction_dip->Draw("SAME");
    
  // doublefermiFunction_dip2->Draw("SAME");  
  // gaus_dip->Draw("SAME");  
 //  line.Draw("SAME");
  // line2.Draw("SAME");
   //line3.Draw("SAME"); 
 // cout << outvect[0] << "##########"<< endl;  
    return outvect;
}//end fitter








// codes that does the job per file. 
vector <float> findcenter_file(const char * _filename="flat.scanx.yoff.chip-1.channel-A1"){
     //setprecision(2);
    // const char * myfname = _filename+3;
    // cout<<  "my file name "<< myfname << endl;
     const char *filename = Form("%s/%s.txt.tree.root",input_loc, _filename);
     cout<< "reading tree : " << filename<< " creating tree " << _filename << endl;
     int n_mes = 10; // number of measurement done for each position.
    
    // this file is for the second test, we have 10 measurements for each position ,, 
    TFile * File = new TFile(filename); // first sipm
    
    TH1D* ratehist= new TH1D("ratehist","",30,0,3000);
    TFile* outfile_raw_data = new TFile(Form("%s/raw_data.root",output_loc),"UPDATE"); 
    if(!File){
        cout << "can not find the file " << filename<<endl;
        return std::vector<float>();
    }
    // getting treee from each file 
    TTree *tree=(TTree*) File->Get("tree");
    
    // each tree have 3 branches x , y and rate 
    float x, y, rate;
    vector <float> vavx;
    vector <float> vavrate;
    vector <float> err_vavrate;
    
     
    // defining branches for each file
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("rate", &rate);

    

    //setting some initial x and y limits to set graph frame to make sure the 
    float xlim_min =0.;
    float xlim_max =-1000.;  

    //
    float ylim_min = 10000000.;
    float ylim_max  = 0.;   
    TGraphErrors * RatevsX = new TGraphErrors();

    
    
    // loop over measurements , 
    for (int i = 0 ; i<tree->GetEntries();i=i+n_mes){
        
        
        //bool remove = false;
        
        if(rate<ylim_min) ylim_min = rate;
        if(rate> ylim_max) ylim_max = rate;
    
        std::vector <float> vrate; //vector of rate
        std:: vector <float> vx; // vector of x
        

        for (int j = 0 ; j<n_mes ; j++){

            tree->GetEntry(i+j);

            if(rate<0.1) continue; // flag to remove measruement with zero rate
            vrate.push_back(rate);
            vx.push_back(x);

             //if(x<-110 || x>-95){
            if((j+i)<200){
                 ratehist->Fill(rate);
            }
       
        }
        
        //if(remove) continue; //removed the zero rate events
        
        
       
        if(x<xlim_min) xlim_min = x;
      
        if(x > xlim_max) xlim_max = x;

        
        
        
        // cout << vrate.size()<< endl; // checking if the size is correct should be 10 in my case because we have 10 measurements
   
        float average_rate = accumulate( vrate.begin(), vrate.end(), 0.0)/vrate.size(); // computing the averate
       // cout << average_rate<< endl; // checking if average is computed correctly

        float std_dev_rate = TMath::StdDev(vrate.begin(), vrate.end());
        float err_ave_rate = std_dev_rate/TMath::Sqrt(vrate.size());
      
        //  cout << "check calculations "<< "  " << std_dev_rate<< endl; // just to check calculcation if correct
   
        float av_x = accumulate( vx.begin(), vx.end(), 0.0)/vx.size(); // computing the averate 
   
        vavx.push_back(av_x);
        vavrate.push_back(average_rate);
         err_vavrate.push_back(err_ave_rate);
       
       // cout<< ((i+1)/10)<< endl;
        int graph_entry = (i+1)/n_mes; // just to fill one value per 



        //cout << graph_entry<< "    "<<i << endl;
        RatevsX->SetPoint(graph_entry , av_x, average_rate);

        RatevsX->SetPointError(graph_entry , error_in_x, err_ave_rate);
        
        

    
    }// end get Entry A1
    
    
    // computing the bkg initial value
    TCanvas * cani = new TCanvas("cani", "", 1200,1200);
    cani->cd();
    TF1 * mygaus_bkg =new  TF1("mygaus_bkg","gaus" , 500,2500);
    ratehist->Fit("mygaus_bkg","Q");
    cani->SaveAs(Form("%s/%s_bkg_fit.png",image_loc,_filename));
    
    float bkg_val = mygaus_bkg->GetParameter(1);
    
    float bkg_err = mygaus_bkg->GetParameter(2);

    float threshold_val = 10000;
    
        
    float xmin = 1000;
    float xmax = -1000;
    float maxrateval = -1000; 
    float error_max_rate = 0;
    //cout<< vavrate.size()<< "#####"<<endl;    
    for(int k=0; k<vavrate.size();k++){
        // cout<< vavrate[k]<< "###########"<< threshold_val << endl;
        if(vavrate[k]>threshold_val){
         //   cout<< rate<< endl;
          //  cout<< vavx[k]<< "###########"<< xmin << endl;
           if(vavx[k]<xmin){
                xmin = vavx[k] ;
            }
            
            if(vavx[k]>xmax){
                xmax=vavx[k];
            }
            if(vavrate[k]>maxrateval){
                maxrateval = vavrate[k];
                error_max_rate = err_vavrate[k];
            }
    
    
        }
    }
    float dip_val = maxrateval;
    float err_dip_val =0;
     float centx = (xmin +xmax)*0.5;
    float x_min_ = centx-0.4; // slighly more than x min 
    float x_max_ = centx+ 0.4; //slightly less than xmax to get the peak region   
    //for the dip value
    for(int k=0; k<vavrate.size();k++){
        // cout<< vavrate[k]<< "###########"<< threshold_val << endl;
            
            if(vavx[k]<x_min_ || vavx[k]> x_max_) continue;
            if(vavrate[k]<dip_val){
                dip_val = vavrate[k];
                err_dip_val = err_vavrate[k];
            }// to find the dip

    }
    
    
    float ave_max_rate = 0;
    float npeak = 0;
     
    for(int k=0; k<vavrate.size();k++){
        // cout<< value of dip is << "###########"<< threshold_val << endl;
        
        if(TMath::Abs(vavrate[k]-maxrateval)<= error_max_rate){
          ave_max_rate += vavrate[k];
          npeak++;    
              
    
        }
    }   
    ave_max_rate = ave_max_rate/npeak;
    
    
  
        
    
    cout << "############ x max "<< xmax << " x min " << xmin<< "bkg "<< bkg_val  << " maximum rate  "<<maxrateval <<endl;
    
    TCanvas * can = new TCanvas("can", "", 1200,1200);
    
    can->cd();
    RatevsX->SetMarkerStyle(29);
    can->SetLeftMargin(0.2);
    can->DrawFrame(xmin-5,ylim_min,xmax+5,ylim_max);
    //30deg
    TLine * midline = new TLine(centx, 0 , centx,maxrateval);
    midline->SetLineWidth(2);
    midline->SetLineColor(kRed);


    TLine * dipline = new TLine(centx-4, dip_val , centx+4,dip_val);
    dipline->SetLineWidth(2);
    dipline->SetLineColor(kBlue);


    RatevsX->Draw("Z");
    RatevsX->GetXaxis()->SetTitle("X position (mm)");
    
    RatevsX->GetYaxis()->SetTitle("Rate of photon detected ");
    //midline->Draw("SAME");
    //dipline->Draw("SAME");
    
    
    vector <float> _outvect = Fitter(RatevsX,xmin,xmax,bkg_val,maxrateval-bkg_val);// calling the fitfunction here. 
    
    
    //float centx = 0.5 *(xmin+xmax);
  //  _outvect[4] = ave_max_rate-bkg_val; //dip_val -bkg_val; 
   // _outvect[5] = error_max_rate;
    
    
    
    //can->Draw();
    
     outfile_raw_data-> cd();  
     can->SaveAs(Form("%s/%saveragerates_vs_x.png",image_loc,_filename));
     can->SaveAs(Form("%s/%saveragerates_vs_x.pdf",image_loc,_filename));
    
    RatevsX->SetName(Form("%saveragerates_vs_x",_filename));
    RatevsX->Write();
    outfile_raw_data->Close();
    return _outvect;

}// end of findcenter

float compute_err_Angle(float width_flat_A1, float error_width_flat_A1) {
    //gROOT->gErrorIgnoreLevel= kWarning;
    
   
    //gErrorIgnoreLevel = 
    //gErrorIgnoreLevel = kInfo;
    // Compute angle
    float angle = TMath::ACos(width_flat_A1 / 3.0);
    if(isnan(angle)) angle =0;
    // Compute derivative (d(angle) / d(width_flat_A1))
    float derivative = -1.0 / (TMath::Sqrt(1.0 - (width_flat_A1 / 3.0) * (width_flat_A1 / 3.0)) * 3.0);

    // Compute error in angle using error propagation
    float error_angle = TMath::Abs(derivative) * error_width_flat_A1;
    if(isnan(error_angle)) error_angle =0;
    return error_angle;
}


float propagateErrorDivision(float A, float sigma_A, float B, float sigma_B) {
    // Calculate the quotient
    double C = A / B;

    // Calculate the propagated error using the formula
    float sigma_C = std::abs(C) * std::sqrt((sigma_A / A) * (sigma_A / A) + (sigma_B / B) * (sigma_B / B));

    return sigma_C;
}


void angle_measurement(bool debug = true){
    
     
   
     TFile * outfile = new TFile(Form("%s/outfile_2.root",output_loc), "RECREATE");
    
    
    setProfessionalPalette();    
    
    vector <float> av_rate_vect;
    vector <float> err_av_rate_vect;  
    vector <float> av_ang_vect;
    vector <float> err_av_ang_vect;   
    
    // graph for targets 
    TGraphErrors *target_A1 = new TGraphErrors();
    TGraphErrors *target_A2 = new TGraphErrors();
    TGraphErrors *target_A3 = new TGraphErrors();
    TGraphErrors *target_A4 = new TGraphErrors();
    target_A1->SetName("target_A1");
    target_A2->SetName("target_A2");
    target_A3->SetName("target_A3");
    target_A4->SetName("target_A4");
   // for reference 
    TGraphErrors *ref_A1 = new TGraphErrors();
    ref_A1->SetName("ref_A1");
    
    
    
    
    // normalized per sensors
    TGraphErrors *target_A1_N = new TGraphErrors();
    TGraphErrors *target_A2_N = new TGraphErrors();
    TGraphErrors *target_A3_N = new TGraphErrors();
    TGraphErrors *target_A4_N = new TGraphErrors();
    target_A1_N->SetName("target_A1_N");
    target_A2_N->SetName("target_A2_N");
    target_A3_N->SetName("target_A3_N");
    target_A4_N->SetName("target_A4_N");
    
    TGraphErrors *ref_A1_N = new TGraphErrors();
    ref_A1_N->SetName("ref_A1");
    
    
    
    //angle vs rate of each sensors
    TGraphErrors *graph_angle_vs_rate_A1 = new TGraphErrors();
    TGraphErrors *graph_angle_vs_rate_A2 = new TGraphErrors();
    TGraphErrors *graph_angle_vs_rate_A3 = new TGraphErrors();
    TGraphErrors *graph_angle_vs_rate_A4 = new TGraphErrors();
    
    TGraphErrors *graph_avangle_vs_rate = new TGraphErrors();
    graph_avangle_vs_rate->SetName("Average_of_all_4_sensors");
    
    
    
    float flat_rate_A1 =0;
    float flat_rate_A2 =0;
    float flat_rate_A3 =0;
    float flat_rate_A4 =0;
    float err_flat_rate_A1 =0;
    float err_flat_rate_A2 =0;
    float err_flat_rate_A3 =0;
    float err_flat_rate_A4 =0;
    
    
    
   
for (int nangle = 0 ; nangle<numangle ; nangle++){
    
    float init_width = 3.0; 
    
    /*vect 0 = max 
      vect 1 = min 
      vect 2 = max error
      vect 3 = min error
      vect 4 = rate
      vect 5 = rate error
      // are the elements of vector float returns
    flat.scanx.yoff.chip-1.channel-A1
    */ 
    
    // calling the file to compute rate , minimum and maximum of peak from each code 
    vector <float> outputvect_A1_ref = findcenter_file(Form("%s/%s.scanx.yoff.%s.channel-A1",ang_name[nangle],ang_name[nangle],"chip-0"));
    // this is uncalibrated data 
    vector <float> outputvect_A1_ = findcenter_file(Form("%s/%s.scanx.yoff.%s.channel-A1",ang_name[nangle],ang_name[nangle],"chip-1"));

    vector <float> outputvect_A2_ = findcenter_file(Form("%s/%s.scanx.yoff.%s.channel-A2",ang_name[nangle],ang_name[nangle],"chip-1"));
    
    vector <float> outputvect_A3_ = findcenter_file(Form("%s/%s.scanx.yoff.%s.channel-A3",ang_name[nangle],ang_name[nangle],"chip-1"));
    
    vector <float> outputvect_A4_ = findcenter_file(Form("%s/%s.scanx.yoff.%s.channel-A4",ang_name[nangle],ang_name[nangle],"chip-1")); 
    
    
    
    
    
    // calibrating the data here
    float calibrated_rate_A1 = outputvect_A1_[4]/ outputvect_A1_ref[4];
    float err_calibrated_rate_A1 =  propagateErrorDivision(outputvect_A1_[4],outputvect_A1_[5],outputvect_A1_ref[4],outputvect_A1_ref[5]);
    
    
    
    float calibrated_rate_A2 = outputvect_A2_[4]/ outputvect_A1_ref[4];
    float err_calibrated_rate_A2 =  propagateErrorDivision(outputvect_A2_[4],outputvect_A2_[5],outputvect_A1_ref[4],outputvect_A1_ref[5]);
    
    float calibrated_rate_A3 = outputvect_A3_[4]/ outputvect_A1_ref[4];
    float err_calibrated_rate_A3 =  propagateErrorDivision(outputvect_A3_[4],outputvect_A3_[5],outputvect_A1_ref[4],outputvect_A1_ref[5]);
    
    
    float calibrated_rate_A4 = outputvect_A4_[4]/ outputvect_A1_ref[4];
    float err_calibrated_rate_A4 =  propagateErrorDivision(outputvect_A4_[4],outputvect_A4_[5],outputvect_A1_ref[4],outputvect_A1_ref[5]);
    
    
    // new vector with calibrated rate

    vector <float> outputvect_A1 ={outputvect_A1_[0],outputvect_A1_[1],outputvect_A1_[2],outputvect_A1_[3],calibrated_rate_A1,err_calibrated_rate_A1};
    
    
    vector <float> outputvect_A2 ={outputvect_A2_[0],outputvect_A2_[1],outputvect_A2_[2],outputvect_A2_[3],calibrated_rate_A2,err_calibrated_rate_A2};
    
    
    vector <float> outputvect_A3 ={outputvect_A3_[0],outputvect_A3_[1],outputvect_A3_[2],outputvect_A3_[3],calibrated_rate_A3,err_calibrated_rate_A3};
    
    vector <float> outputvect_A4 ={outputvect_A4_[0],outputvect_A4_[1],outputvect_A4_[2],outputvect_A4_[3],calibrated_rate_A4,err_calibrated_rate_A4};
    
    

    
   if(debug){
     cout << Form(" \n \n \n For A1 Rate is %f +- %f before calibration \n and rate of reference %f +- %f \n and Calibrated rate is %f +- %f \n \n \n \n",outputvect_A1_[4],outputvect_A1_[5],outputvect_A1_ref[4],outputvect_A1_ref[5],outputvect_A1[4],outputvect_A1[5] ) << endl;
    
    cout << Form(" \n \n \n For A2 Rate is %f +- %f before calibration \n and rate of reference %f +- %f \n and Calibrated rate is %f +- %f \n \n \n \n",outputvect_A2_[4],outputvect_A2_[5],outputvect_A1_ref[4],outputvect_A1_ref[5],outputvect_A2[4],outputvect_A2[5] ) << endl;
    
    cout << Form(" \n \n \n For A3 Rate is %f +- %f before calibration \n and rate of reference %f +- %f \n and Calibrated rate is %f +- %f \n \n \n \n",outputvect_A3_[4],outputvect_A3_[5],outputvect_A1_ref[4],outputvect_A1_ref[5],outputvect_A3[4],outputvect_A3[5] ) << endl;
    
    cout << Form(" \n \n \n For A4 Rate is %f +- %f before calibration \n and rate of reference %f +- %f \n and Calibrated rate is %f +- %f \n \n \n \n",outputvect_A4_[4],outputvect_A4_[5],outputvect_A1_ref[4],outputvect_A1_ref[5],outputvect_A4[4],outputvect_A4[5] ) << endl;
    
    
    }
    
     
    
    
    
    // computing width and angle 
    // For A1
    float width_A1 =  outputvect_A1[0]- outputvect_A1[1];    
    float error_width_A1 = TMath::Sqrt(outputvect_A1[2]*outputvect_A1[2] +outputvect_A1[3]*outputvect_A1[3]);         

    float angle_A1 =    TMath::ACos(width_A1/init_width);
        if(isnan(angle_A1)) angle_A1 = 0;
    float error_angle_A1 =compute_err_Angle(width_A1,error_width_A1);    

    //for A2   
    float width_A2 =  outputvect_A2[0]- outputvect_A2[1];    
    float error_width_A2 = TMath::Sqrt(outputvect_A2[2]*outputvect_A2[2] +outputvect_A2[3]*outputvect_A2[3]);         

    float angle_A2 =    TMath::ACos(width_A2/init_width);
        if(isnan(angle_A2)) angle_A2 = 0;
    float error_angle_A2 =compute_err_Angle(width_A2,error_width_A2);    


    
    
    // For A3
    float width_A3 =  outputvect_A3[0]- outputvect_A3[1];    
    float error_width_A3 = TMath::Sqrt(outputvect_A3[2]*outputvect_A3[2] +outputvect_A3[3]*outputvect_A3[3]);         
    
    float angle_A3 =    TMath::ACos(width_A3/init_width);
    if(isnan(angle_A3)) angle_A3 = 0;
    float error_angle_A3 =compute_err_Angle(width_A3,error_width_A3);    



           

    //  For A4

    float width_A4 =  outputvect_A4[0]- outputvect_A4[1];    
    float error_width_A4 = TMath::Sqrt(outputvect_A4[2]*outputvect_A4[2] +outputvect_A4[3]*outputvect_A4[3]);         

    float angle_A4 =    TMath::ACos(width_A4/init_width);
    if(isnan(angle_A4)) angle_A4 = 0;
    float error_angle_A4 =compute_err_Angle(width_A4,error_width_A4);    

    
    
    
    
    target_A1->SetPoint(nangle,angle_A1*angle_to_theta,outputvect_A1_[4]);
    target_A2->SetPoint(nangle,angle_A2*angle_to_theta,outputvect_A2_[4]);
    target_A3->SetPoint(nangle,angle_A3*angle_to_theta,outputvect_A3_[4]);
    target_A4->SetPoint(nangle,angle_A4*angle_to_theta,outputvect_A4_[4]);
    ref_A1->SetPoint(nangle,angle_A1*angle_to_theta,outputvect_A1_ref[4]);
    
    
    target_A1->SetPointError(nangle,error_angle_A1*angle_to_theta,outputvect_A1_[5]);
    target_A2->SetPointError(nangle,error_angle_A2*angle_to_theta,outputvect_A2_[5]);
    target_A3->SetPointError(nangle,error_angle_A3*angle_to_theta,outputvect_A3_[5]);
    target_A4->SetPointError(nangle,error_angle_A4*angle_to_theta,outputvect_A4_[5]);
    //ref_A1->SetPointError(nangle,error_angle_A1_ref,outputvect_A1_ref[5]);
    
    
    target_A1_N->SetPoint(nangle,angle_A1*angle_to_theta,outputvect_A1[4]);
    target_A2_N->SetPoint(nangle,angle_A2*angle_to_theta,outputvect_A2[4]);
    target_A3_N->SetPoint(nangle,angle_A3*angle_to_theta,outputvect_A3[4]);
    target_A4_N->SetPoint(nangle,angle_A4*angle_to_theta,outputvect_A4[4]);
    ref_A1_N->SetPoint(nangle,angle_A1*angle_to_theta,outputvect_A1_ref[4]/outputvect_A1_ref[4]);
    
    
    target_A1_N->SetPointError(nangle,error_angle_A1*angle_to_theta,outputvect_A1[5]);
    target_A2_N->SetPointError(nangle,error_angle_A2*angle_to_theta,outputvect_A2[5]);
    target_A3_N->SetPointError(nangle,error_angle_A3*angle_to_theta,outputvect_A3[5]);
    target_A4_N->SetPointError(nangle,error_angle_A4*angle_to_theta,outputvect_A4[5]); 
    
    


   
    if (debug) {
         
         cout << "\n \n \n for "<< ang_name[nangle]<< endl; 
        cout << " width " << " sensor name " << "errro in width " << " angle " << "  error in angle  " <<endl;
         cout<< width_A1 << "  A1  " << error_width_A1<< "   "<<angle_A1*angle_to_theta<< "  "<< error_angle_A1*angle_to_theta << endl;
         cout<< width_A2 << "  A2  " << error_width_A2<< "   "<<angle_A2*angle_to_theta<< "  "<< error_angle_A2*angle_to_theta << endl;
         cout<< width_A3 << "  A3  " << error_width_A3<< "   "<<angle_A3*angle_to_theta<< "  "<< error_angle_A3*angle_to_theta << endl; 
         cout<< width_A4 << "  A4  " << error_width_A4<< "   "<<angle_A4*angle_to_theta<< " "<< error_angle_A4*angle_to_theta << endl;
         
         cout << "\n \n \n "<< endl;
         
         
    
     }
      
    


    float av_angle =      (angle_A1 + angle_A2+ angle_A3 +angle_A4)/4.0;


     float err_av_angle = TMath::Sqrt(error_angle_A1* error_angle_A1+error_angle_A2* error_angle_A2+error_angle_A3* error_angle_A3 +error_angle_A4* error_angle_A4)*0.25 ;  

     //normalizing the rate per sensors
    
    //this is to save the rate of flat to normalize
    if(std::strcmp(ang_name[nangle],"flat")==0) {
      flat_rate_A1 = outputvect_A1[4] ;
      flat_rate_A2 = outputvect_A2[4] ;
      flat_rate_A3 = outputvect_A3[4] ;
      flat_rate_A4 = outputvect_A4[4] ;
        
      err_flat_rate_A1 = outputvect_A1[5] ;
      err_flat_rate_A2 = outputvect_A2[5] ;
      err_flat_rate_A3 = outputvect_A3[5] ;
      err_flat_rate_A4 = outputvect_A4[5] ;  
        
        
      has_flat = true;
    }
    
        float Normalized_A1;
        float Normalized_A2;
        float Normalized_A3;
        float Normalized_A4; 
    
   
    
        float err_Normalized_A1;

        float err_Normalized_A2;


        float err_Normalized_A3;

        float err_Normalized_A4;

    
    if(has_flat){
    
        Normalized_A1 = outputvect_A1[4]/flat_rate_A1;
        Normalized_A2 = outputvect_A2[4]/flat_rate_A2;
        Normalized_A3 = outputvect_A3[4]/flat_rate_A3;
        Normalized_A4 = outputvect_A4[4]/flat_rate_A4; 
    
   
    
        err_Normalized_A1 = propagateErrorDivision(outputvect_A1[4],outputvect_A1[5],flat_rate_A1,err_flat_rate_A1);

        err_Normalized_A2 = propagateErrorDivision(outputvect_A2[4],outputvect_A2[5],flat_rate_A2,err_flat_rate_A2);


        err_Normalized_A3 = propagateErrorDivision(outputvect_A3[4],outputvect_A3[5],flat_rate_A3,err_flat_rate_A3);

        err_Normalized_A4 = propagateErrorDivision(outputvect_A4[4],outputvect_A4[5],flat_rate_A4,err_flat_rate_A4);

    
    
    }
    
    else{
    
        Normalized_A1 = outputvect_A1[4];
        Normalized_A2 = outputvect_A2[4];
        Normalized_A3 = outputvect_A3[4];
        Normalized_A4 = outputvect_A4[4]; 
    
   
    
        err_Normalized_A1 = outputvect_A1[5];
        err_Normalized_A2 = outputvect_A2[5];


        err_Normalized_A3 = outputvect_A3[5];
        err_Normalized_A4 = outputvect_A4[5];
    
    
    }
    

    float av_rate = (Normalized_A1 + Normalized_A2 + Normalized_A3 + Normalized_A4)/4.0;
    float err_av_rate = 0.25 * TMath::Sqrt(err_Normalized_A1*err_Normalized_A1 + err_Normalized_A2*err_Normalized_A2 + err_Normalized_A3*err_Normalized_A3 +err_Normalized_A4*err_Normalized_A4);
    
    av_rate_vect.push_back(av_rate);
    err_av_rate_vect.push_back(err_av_rate);
    
    av_ang_vect.push_back(av_angle);
    err_av_ang_vect.push_back(err_av_angle);
    
    
    
       
     if(debug){
         
              cout << " Checking  normalization per sensor with respect to initial flat sensor " << endl;
     
             cout << Form("\n \n \n %f / %f =  Normalized_A1 %f +- %f \n %f / %f =  Normalized_A2 %f +- %f \n %f / %f =  Normalized_A3 %f+- %f \n %f / %f =  Normalized_A4 %f+-%f \n\n\n " ,  outputvect_A1[4],flat_rate_A1,Normalized_A1,err_Normalized_A1,outputvect_A2[4],flat_rate_A2,Normalized_A2,err_Normalized_A2,outputvect_A3[4],flat_rate_A3,Normalized_A3,err_Normalized_A3,outputvect_A4[4],flat_rate_A4,Normalized_A4,err_Normalized_A4)<< endl;


          cout << Form("\n \n \n average normalized rate %f +-%f \n",av_rate , err_av_rate) <<endl;

     
     
     } 
    
    
     graph_angle_vs_rate_A1->SetPoint(nangle,angle_A1,Normalized_A1);
     graph_angle_vs_rate_A1->SetPointError(nangle,error_angle_A1,err_Normalized_A1);
     
     graph_angle_vs_rate_A2->SetPoint(nangle,angle_A2,Normalized_A2);
     graph_angle_vs_rate_A2->SetPointError(nangle,error_angle_A2,err_Normalized_A2);
      
     graph_angle_vs_rate_A3->SetPoint(nangle,angle_A3,Normalized_A3);
     graph_angle_vs_rate_A3->SetPointError(nangle,error_angle_A3,err_Normalized_A3);
    
    
     graph_angle_vs_rate_A4->SetPoint(nangle,angle_A4,Normalized_A4);
     graph_angle_vs_rate_A4->SetPointError(nangle,error_angle_A4,err_Normalized_A4);
     
     
    
     
    
    
     graph_avangle_vs_rate->SetPoint(nangle,av_angle, av_rate);
     graph_avangle_vs_rate->SetPointError(nangle,err_av_angle, err_av_rate);
   
    


}// end n angle loop     

    
    //cout<< " Average Rate 40deg:  " << av_rate_40deg << " errror average rate 40deg:  "<< err_av_rate_40deg << endl;    
    
 
    
    
    
    
   
    
    
    TGraphErrors *graph_av_angle_vs_normrate = new TGraphErrors();
    float av_rate_flat=0;
    float err_av_rate_flat=0;
    for (int i =0; i< numangle ;i++){
        float av_rate_ = av_rate_vect[i];
        float err_av_rate_ = err_av_rate_vect[i];
        float av_angle_ = av_ang_vect[i];
        float err_av_angle_ = err_av_ang_vect[0];
        
        if(std::strcmp(ang_name[i],"flat")==0) {
            av_rate_flat = av_rate_ ;
            err_av_rate_flat = av_rate_;    
        }
        
        float Norm_rate_ = av_rate_;//av_rate_flat;
        float err_Norm_rate_ =err_av_rate_;//0; //propagateErrorDivision(av_rate_flat,err_av_rate_flat,av_rate_flat,err_av_rate_flat);
            
        graph_av_angle_vs_normrate->SetPoint(i,av_angle_*angle_to_theta, Norm_rate_); 
        graph_av_angle_vs_normrate->SetPointError(0,err_av_angle_*angle_to_theta, err_Norm_rate_);
        
        
         if(debug){
             cout<< Form("\n \n \naverage_angle of %s is :  %f +- %f ",ang_name[i], av_angle_,err_av_angle_)<< endl; 
         cout<< Form("average_rate of %s is :  %f +- %f ",ang_name[i], av_rate_,err_av_rate_)<< endl; 
         }//END DEBUG
         
    }
    
     TCanvas* canvas_av_ang_vs_normrate = new TCanvas("canvas_av_ang_vs_normrate","",1600,1200);
     canvas_av_ang_vs_normrate->cd();
     canvas_av_ang_vs_normrate->DrawFrame(0,0,1,1);
     graph_av_angle_vs_normrate->SetLineColor(46);
     graph_av_angle_vs_normrate->SetMarkerStyle(24);
     graph_av_angle_vs_normrate->SetLineWidth(3);
     
     graph_av_angle_vs_normrate->GetXaxis()->SetTitle("Incidence angle (deg)");
     graph_av_angle_vs_normrate->GetYaxis()->SetTitle(" Relative photodetection efficiency");
     graph_av_angle_vs_normrate->SetName("Normalized_rate ");
     
         
      graph_av_angle_vs_normrate->Draw("AP");
     canvas_av_ang_vs_normrate->SetLeftMargin(0.2); 
     //canvas_av_ang_vs_normrate->BuildLegend();   
     canvas_av_ang_vs_normrate->SaveAs(Form("%s/average_angle_vs_normrate.png",image_loc)); 
     canvas_av_ang_vs_normrate->SaveAs(Form("%s/average_angle_vs_normrate.pdf",image_loc));   


    
     TCanvas* canvas_av_ang_vs_rate = new TCanvas("canvas_av_ang_vs_rate","",1600,1200);
     canvas_av_ang_vs_rate->cd();
     canvas_av_ang_vs_rate->DrawFrame(0,10000,1,25000);
     graph_avangle_vs_rate->SetLineColor(46);
     graph_avangle_vs_rate->SetMarkerStyle(24);
     graph_avangle_vs_rate->SetLineWidth(3);

     graph_avangle_vs_rate->GetXaxis()->SetTitle("Average angle");
     graph_avangle_vs_rate->GetYaxis()->SetTitle("Average rate");   

     
     graph_avangle_vs_rate->Draw("AP");
     canvas_av_ang_vs_rate->SetLeftMargin(0.2); 
     canvas_av_ang_vs_rate->BuildLegend();   
     canvas_av_ang_vs_rate->SaveAs(Form("%s/average_angle_vs_rate.png",image_loc)); 
     canvas_av_ang_vs_rate->SaveAs(Form("%s/average_angle_vs_rate.pdf",image_loc));   
    
    
    
    TMultiGraph *mg3 = new TMultiGraph();
    
    TMultiGraph *mg4 = new TMultiGraph();
    mg4->SetName("calibrated_rate");
     
     target_A1->SetMarkerStyle(25);
     target_A1->SetLineColor(42);
     target_A1->SetLineWidth(3);   


     target_A2->SetMarkerStyle(24);
     target_A2->SetLineColor(36);
     target_A2->SetLineWidth(3);   


     target_A3->SetMarkerStyle(27);
     target_A3->SetLineColor(30);
     target_A3->SetLineWidth(3);   


     target_A4->SetMarkerStyle(28);
     target_A4->SetLineColor(kRed);
     target_A4->SetLineWidth(3);   

    
     ref_A1->SetMarkerStyle(30);
     ref_A1->SetLineColor(52);
     ref_A1->SetLineWidth(3); 
    
    
    
    target_A1_N->SetMarkerStyle(25);
     target_A1_N->SetLineColor(42);
     target_A1_N->SetLineWidth(3);   


     target_A2_N->SetMarkerStyle(24);
     target_A2_N->SetLineColor(36);
     target_A2_N->SetLineWidth(3);   


     target_A3_N->SetMarkerStyle(27);
     target_A3_N->SetLineColor(30);
     target_A3_N->SetLineWidth(3);   


     target_A4_N->SetMarkerStyle(28);
     target_A4_N->SetLineColor(kRed);
     target_A4_N->SetLineWidth(3);   

    
     ref_A1_N->SetMarkerStyle(30);
     ref_A1_N->SetLineColor(52);
     ref_A1_N->SetLineWidth(3); 
    
    mg3->Add(target_A1);
    mg3->Add(target_A2);
    mg3->Add(target_A3);
    mg3->Add(target_A4);
    mg3->Add(ref_A1);
    
    
    mg4->Add(target_A1_N);
    mg4->Add(target_A2_N);
    mg4->Add(target_A3_N);
    mg4->Add(target_A4_N);
    mg4->Add(ref_A1_N);
    
    
     
    
    TCanvas* can_org = new TCanvas("can_org","",1600,1200);
    can_org->cd();
    mg3->Draw("AP");
    can_org->BuildLegend();
    can_org->SaveAs(Form("%s/original_rate.png",image_loc));
    
    
    TCanvas* can_org_cal = new TCanvas("can_org_cal","",1600,1200);
    
    can_org_cal->cd();
    can_org_cal->DrawFrame(0,0.8,60,1.2);
    mg4->Draw("AP");
    can_org_cal->BuildLegend();
    can_org_cal->SaveAs(Form("%s/calibrated_rate.png",image_loc));
    
    
      
    
 
    
     TCanvas* canvas_ang_vs_rate = new TCanvas("canvas_ang_vs_rate","",1600,1200);
     TMultiGraph *mg = new TMultiGraph();
     //mg->SetTitle("");   

     canvas_ang_vs_rate->cd();
     canvas_ang_vs_rate->DrawFrame(0,10000,1,25000);   

     graph_angle_vs_rate_A1->SetMarkerStyle(25);
     graph_angle_vs_rate_A1->SetLineColor(42);
     graph_angle_vs_rate_A1->SetLineWidth(3);   


     graph_angle_vs_rate_A2->SetMarkerStyle(24);
     graph_angle_vs_rate_A2->SetLineColor(36);
     graph_angle_vs_rate_A2->SetLineWidth(3);   


     graph_angle_vs_rate_A3->SetMarkerStyle(27);
     graph_angle_vs_rate_A3->SetLineColor(30);
     graph_angle_vs_rate_A3->SetLineWidth(3);   


     graph_angle_vs_rate_A4->SetMarkerStyle(28);
     graph_angle_vs_rate_A4->SetLineColor(41);
     graph_angle_vs_rate_A4->SetLineWidth(3);   


      mg->Add(graph_angle_vs_rate_A1);
      mg->Add(graph_angle_vs_rate_A2);
      mg->Add(graph_angle_vs_rate_A3);
      mg->Add(graph_angle_vs_rate_A4);


      canvas_ang_vs_rate->SetLeftMargin(0.2);  
      //mg->Add(graph_angle_vs_rate_A5);  
     //graph_angle_vs_rate_A1->Draw("ALP"); 
     //graph_angle_vs_rate_A2->Draw("ALP SAME");

     //graph_angle_vs_rate_A3->Draw("ALP SAME"); 
     //graph_angle_vs_rate_A4->Draw("ALP SAME");
     graph_angle_vs_rate_A1->SetName("Sensor_A1");
     graph_angle_vs_rate_A2->SetName("Sensor_A2");
     graph_angle_vs_rate_A3->SetName("Sensor_A3");
     graph_angle_vs_rate_A4->SetName("Sensor_A4");   

     mg->Draw("AP");
     mg->GetYaxis()->SetTitle("Rate of photon detection");
     mg->GetXaxis()->SetTitle("Angle of incidence");      
     //mg->GetHistogram()->SetTitle("Global title");   
     canvas_ang_vs_rate->BuildLegend(); 
     canvas_ang_vs_rate->SaveAs(Form("%s/angle_rate_check_final.png",image_loc));
     canvas_ang_vs_rate->SaveAs(Form("%s/angle_rate_check_final.pdf",image_loc));  
    
     
    TMultiGraph *mg2 = new TMultiGraph();
    
    
   // TCanvas *canvas_fit = new TCanvas("canvas_fit","",800,800);
    //canvas_fit->SetLeftMargin(0.2);
   // mg2->Draw("ALP");
    //mg2->GetYaxis()->SetTitle("Normalized rate of photon detection");
    //mg2->GetXaxis()->SetTitle("Cosine of Angle of Incidence");      
    //canvas_fit->SaveAs(Form("%s/cos_check.png",image_loc));
    //outfile->Add()
    //mg2->SetName("graph");
    
    mg3->SetName("graph_of_rate");
    mg4->SetName("graph_of_calibrated_rate");
    outfile->cd();
    mg2->Write();
    graph_av_angle_vs_normrate->Write();
    graph_avangle_vs_rate->Write();
    graph_angle_vs_rate_A1->Write();
    graph_angle_vs_rate_A2->Write();
    graph_angle_vs_rate_A3->Write();
    graph_angle_vs_rate_A4->Write();
    mg3->Write();
    mg4->Write();
    outfile->Close();
    //outfile->Write();
}


