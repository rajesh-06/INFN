#include "iostream"
#include "fstream"

using namespace std;


//gROOT->gErrorIgnoreLevel = kFatal;// kInfo, kWarning, kError, kBreak, kSysError, kFatal;
const char* image_loc = "image"; // dir name for output image
const char* input_loc = "Data"; // dir name for input data
const char* output_loc = "output"; // dirname for output data TGraph
float angle_to_theta  = 180.0/TMath::Pi();// for conversion to angle 
int usefit = 1; // 1 if double fermi , 2 if double erf , 3 if double fermi with dip
float has_flat=false;

float error_in_x = 0.05; // increment of the x position
bool debug_fit = true;
//vector <float> peakrate = {}



int numangle = 99;//ang_name.size();//ang_name.size();// setting number of angle measurement currently we have only 2
    
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
        outvect.push_back(doublefermiFunction->GetParameter(0));
        outvect.push_back(doublefermiFunction->GetParError(0));
        if (debug_fit) cout<< "DOUBLE FERMI FUNCTION "<< endl;
        _graph->Fit("doublefermiFunction",op);
    

    }
     
    else if (usefit==2){
    
        outvect.push_back(doublefermiFunction_dip->GetParameter(3));
        outvect.push_back(doublefermiFunction_dip->GetParameter(2));
        outvect.push_back(doublefermiFunction_dip->GetParError(3));
        outvect.push_back(doublefermiFunction_dip->GetParError(2));
        outvect.push_back(doublefermiFunction_dip->GetParameter(1));
        outvect.push_back(doublefermiFunction_dip->GetParError(1));
        if (debug_fit) cout<< "DOUBLE ERF FUNCTION "<< endl;
         _graph->Fit("double_erf", op);

    
    
    }
    
   else if (usefit ==3) {
        outvect.push_back(doublefermiFunction_dip->GetParameter(3));
        outvect.push_back(doublefermiFunction_dip->GetParameter(2));
        outvect.push_back(doublefermiFunction_dip->GetParError(3));
        outvect.push_back(doublefermiFunction_dip->GetParError(2));
        outvect.push_back(doublefermiFunction_dip->GetParameter(1));
        outvect.push_back(doublefermiFunction_dip->GetParError(1));
        if (debug_fit) cout<< "DOUBLE FERMI FUNCTION  WITH DIP "<< endl;
        _graph->Fit("doublefermiFunction_dip",op);
   
   
       
   
   }
    
   else{
       if (debug_fit) cout<< "DOUBLE FERMI FUNCTION "<< endl;
       _graph->Fit("doublefermiFunction",op);
    
       if (debug_fit) cout<< "DOUBLE ERF FUNCTION "<< endl;
       _graph->Fit("double_erf", op);

       if (debug_fit) cout<< "DOUBLE FERMI FUNCTION  WITH DIP "<< endl;
       _graph->Fit("doublefermiFunction_dip",op);
       
   
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
vector <float> findcenter_file(const char * _filename="flat.scanx.yoff.chip-1.channel-A1", TGraphErrors * pos_rate_gr = new TGraphErrors()){
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
    _outvect[4] = ave_max_rate;//-bkg_val; //dip_val -bkg_val; 
    _outvect[5] = error_max_rate;
    
    
    
    
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
    
     
   
     TFile * outfile = new TFile(Form("%s/outfile_.root",output_loc), "RECREATE");
    
    
    setProfessionalPalette();    
    
    
   // for reference 
    TGraphErrors *ref_A1 = new TGraphErrors();
    ref_A1->SetName("ref_A1");
    
    
    TGraphErrors *ref_A1_corr = new TGraphErrors();
    ref_A1_corr->SetName("ref_A1_corr");
    
    
    TGraphErrors *ref_A1_dc = new TGraphErrors();
    ref_A1_dc->SetName("ref_A1_dc");
    
    TH2D* corr_ref = new TH2D("corr_ref",";DCR; Rate of Sensor",10,1300,1400,10,25600,26200);
    
    
    
    
    
    
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
     vector <float> outputvect_A1_ref;
    // calling the file to compute rate , minimum and maximum of peak from each code 
    if(nangle<10){
        outputvect_A1_ref = findcenter_file(Form("loop-0%i.scanx.yoff.%s.channel-A1",nangle,"chip-0"));
    }
    else {
    
        outputvect_A1_ref = findcenter_file(Form("loop-%i.scanx.yoff.%s.channel-A1",nangle,"chip-0"));
    }
    
    
    
    
    
    //loop-26.scanx.ycen.chip-0.channel-A1.txt.tree.root
    
    
    
   if(debug){
     cout << Form(" \n \n \n loop-%i the rate of reference %f +- %f \n \n \n \n \n",nangle,outputvect_A1_ref[4],outputvect_A1_ref[5]) << endl;
    
    
    }
    
     
    
    
    
    ref_A1->SetPoint(nangle,nangle,outputvect_A1_ref[4]);
    
    ref_A1->SetPointError(nangle,0,outputvect_A1_ref[5]);
    
    ref_A1_corr->SetPoint(nangle,outputvect_A1_ref[6],outputvect_A1_ref[4]);
    
    ref_A1_corr->SetPointError(nangle,outputvect_A1_ref[7],outputvect_A1_ref[5]);
    corr_ref->Fill(outputvect_A1_ref[6],outputvect_A1_ref[4]);
    
    ref_A1_dc->SetPoint(nangle,nangle,outputvect_A1_ref[6]);
    
    ref_A1_dc->SetPointError(nangle,0,outputvect_A1_ref[7]);
    
    
    
    
    
    

}// end n angle loop     

    
    //cout<< " Average Rate 40deg:  " << av_rate_40deg << " errror average rate 40deg:  "<< err_av_rate_40deg << endl;    
    TMultiGraph *mg3 = new TMultiGraph();
    
    
    
     ref_A1->SetMarkerStyle(30);
     ref_A1->SetLineColor(52);
     ref_A1->SetLineWidth(3); 
    
    
     ref_A1_dc->SetMarkerStyle(20);
     ref_A1_dc->SetLineColor(22);
     ref_A1_dc->SetLineWidth(3); 
    
    cout << "Integral is " << ref_A1->Integral(0,1000)<< endl;
    cout << "Integral is dc " << ref_A1_dc->Integral(0,1000)<< endl;
    
    TF1 f1("f",[&](double *x, double *){ return ref_A1->Eval(x[0]); },1,ref_A1->GetN(),0);
    TF1 f2("f",[&](double *x, double *){ return ref_A1_dc->Eval(x[0]); },1,ref_A1_dc->GetN(),0);
        
    //cout << "Integral new  " << f1.Integral(1,ref_A1->GetN(),1e-10)<< endl;
    float integ = f1.Integral(1,ref_A1->GetN(),1e-10);
    float integ2 = f2.Integral(1,ref_A1_dc->GetN(),1e-10);
   // ref_A1->Scale((1.0/integ),"y");
    
   // ref_A1_dc->Scale((1.0/integ2),"y");
    
    mg3->Add(ref_A1);
    
    //mg3->Add(ref_A1_dc);
    
    TCanvas* can_org = new TCanvas("can_org","",1600,1200);
    can_org->cd();
    mg3->Draw("ALP");
    mg3->GetXaxis()->SetTitle("Loop Number");
    
    mg3->GetYaxis()->SetTitle("Rate");
    
        
   // can_org->BuildLegend();
    can_org->SaveAs(Form("%s/original_rate.png",image_loc));
    
    
    
    outfile->cd();
      
    ref_A1->Write();
    
    ref_A1_dc->Write();
    
    mg3->Write();
    
    TCanvas* can_org_C = new TCanvas("can_org_C","",1600,1200);
    can_org_C->SetLeftMargin(0.2);
    gStyle->SetOptStat(0);
    can_org_C->cd();
    ref_A1_corr->Draw("AP");
    //corr_ref->Draw("COL");
    
    can_org_C->SaveAs("correlataion.png");
    corr_ref->Write();
    ref_A1_corr->Write(); 
    outfile->Close();
    //outfile->Write();
}


