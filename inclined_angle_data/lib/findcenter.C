#include "iostream"
#include "fstream"

using namespace std;


vector <float> Fitter(TGraphErrors * _graph, float minvalx =-104., float maxvalx =-101. , float bkg = 1000. , float sigtop = 25000){
    vector <float> fitresult;
    
    //double fermi function

    TF1 *doublefermiFunction = new TF1("doublefermiFunction", [](double *x, double *p) {
            return p[0]+ p[1]/(TMath::Exp((x[0] - p[3]) /p[4]) + 1.0)/(TMath::Exp((p[2]-x[0] ) /p[5]) + 1.0);
    }, -120, -95, 6); // Defining the function from 0 to 10, 0 parameters

    doublefermiFunction->SetParameters(bkg,sigtop,minvalx, maxvalx,0.05,0.05);




        //gaussian error function    
    TF1 * double_erf = new TF1("double_erf", [](double *x, double *p) {
            return p[0] + p[1] * (TMath::Erf((x[0] - p[2]) /p[4]) + 1.) * (TMath::Erf((p[3]-x[0] ) /p[5]) + 1.) * 0.25;
    }, -120, -95, 6); // Defining the function from 0 to 10, 0 parameters


    double_erf->SetNpx(1000);
    doublefermiFunction->SetNpx(1000);    

    double_erf->SetParameters(bkg,sigtop,minvalx,maxvalx,0.1,0.1);
    double_erf->SetLineColor(kGreen+3);
    cout<< "DOUBLE FERMI FUNCTION "<< endl;
    _graph->Fit("doublefermiFunction","R");
    
    
     cout<< "DOUBLE ERF FUNCTION "<< endl;
     _graph->Fit("double_erf","R");


    //fitresult.push_back(doublefermiFunction->GetParameter(0));
    
    //fitresult.push_back(doublefermiFunction->GetParameter(1));
   
    
    fitresult.push_back(doublefermiFunction->GetParameter(2));
    fitresult.push_back(doublefermiFunction->GetParError(2));
    fitresult.push_back(doublefermiFunction->GetParameter(3));
    fitresult.push_back(doublefermiFunction->GetParError(3));
    
    //fitresult.push_back(doublefermiFunction->GetParameter(4));
    //fitresult.push_back(doublefermiFunction->GetParameter(5));
    //fitresult.push_back(doublefermiFunction->GetParameter(6));
    doublefermiFunction->Draw("SAME");
    double_erf->Draw("SAME");  


    return fitresult;
}//end fitter









void findcenter(TString filename = "scanx.200mrad.A1.txt.tree.root"){

     int n_mes = 10; // number of measurement done for each position.
    // this file is for the second test, we have 10 measurements for each position ,, 
    TFile * File = new TFile(filename.Data()); // first sipm
    
    TH1D* ratehist= new TH1D("ratehist","",30,0,3000);
    
    if(!File){
        cout << "can not find the file " << filename.Data()<<endl;
        return;
    }
    // getting treee from each file 
    TTree *tree=(TTree*) File->Get("tree");
    
    // each tree have 3 branches x , y and rate 
    float x, y, rate;
    vector <float> vavx;
    vector <float> vavrate;
     
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

    for (int i = 0 ; i<tree->GetEntries();i=i+n_mes){
        
        
        //bool remove = false;
        
        if(rate<ylim_min) ylim_min = rate;
        if(rate> ylim_max) ylim_max = rate;
    
        std::vector <float> vrate; //vector of rate
        std::vector <float> vx; // vector of x
        

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

        
        
        
        //cout << vrate.size()<< endl; // checking if the size is correct should be 10 in my case because we have 10 measurements
   
        float average_rate = accumulate( vrate.begin(), vrate.end(), 0.0)/vrate.size(); // computing the averate
       //cout << average_rate<< endl; // checking if average is computed correctly

        float std_dev_rate = TMath::StdDev(vrate.begin(), vrate.end());
        float err_ave_rate = std_dev_rate/TMath::Sqrt(vrate.size());
   //  cout << "check calculations "<< "  " << std_dev_rate<< endl; // just to check calculcation if correct
   
        float av_x = accumulate( vx.begin(), vx.end(), 0.0)/vx.size(); // computing the averate 
        vavx.push_back(av_x);
        vavrate.push_back(average_rate);
       
       // cout<< ((i+1)/10)<< endl;
        int graph_entry = (i+1)/n_mes; // just to fill one value per 

        RatevsX->SetPoint(graph_entry , av_x, average_rate);
        RatevsX->SetPointError(graph_entry , 0., err_ave_rate);
 
    }// end get Entry A1
    
    // computing the bkg initial value
    TCanvas * cani = new TCanvas("cani", "", 1200,1200);
    cani->cd();
    TF1 * mygaus_bkg =new  TF1("mygaus_bkg","gaus" , 500,2500);
    ratehist->Fit("mygaus_bkg");
    cani->SaveAs("bkg_fit.png");
    
    float bkg_val = mygaus_bkg->GetParameter(1);
    float bkg_err = mygaus_bkg->GetParameter(2);
    float threshold_val = 10000;
    
        
    float xmin = 1000;
    float xmax = -1000;
    float maxval = -1000;    
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
            if(vavrate[k]>maxval){
                maxval = vavrate[k];
            }
    
    
        }
    }   
        
        
    
    
    
    TCanvas * can = new TCanvas("can", "", 1200,1200);
    can->cd();
    can->DrawFrame(xlim_min,ylim_min,xlim_max,ylim_max);

    RatevsX->Draw("Z");
    // Fitter(TGraphErrors * _graph, float minvalx =-104., float maxvalx =-101. , float bkg = 1000. , float sigtop = 25000)
        
    //vector <float> myresult = Fitter(RatevsX,-104,-101,1000,25000);
    
    vector <float> myresult = Fitter(RatevsX,xmin,xmax,bkg_val,maxval);
    
    for(int i = 0; i< myresult.size(); i++){
        cout<<myresult[i]<<endl;
    }
   // cout<< xmax << "  "<< xmin << endl;
    
    
    TLine line(myresult[0],0,myresult[0],maxval);

    TLine line2(myresult[2],0,myresult[2],maxval);

    line.SetLineColor(12);
    line.SetLineWidth(4);
    line.SetLineStyle(10);

    line2.SetLineColor(32);
    line2.SetLineStyle(10);
    line2.SetLineWidth(4);


    line.Draw("SAME");
    line2.Draw("SAME");
    
    can->Draw();
    
    
    // can->SaveAs("averagerates_vs_x.png");
    
    //cout<< "low x val " <<   myresult[0]<< " +- " <<  myresult[1]<<endl;
    //cout<< "high x val " <<   myresult[2]<< " +- " <<  myresult[3]<<endl;
    

}// end of findcenter


