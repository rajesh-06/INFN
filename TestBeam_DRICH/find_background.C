#include "iostream"
using namespace std;
float center_x = -999;
float center_y =-999;
//float myoutrad_ = -999;
//float myinrad_ = -999;
float myrad_ = -999;
float err_myrad_ = -999;



void fitcircle(TGraphErrors* Graph_XY_10, float est_rad =10){
    cout << "fitting with radius"<< est_rad <<endl;
    auto chi2Function = [&](const double *par) {
      //minimisation function computing the sum of squares of residuals
      // looping at the graph points
      int np = Graph_XY_10->GetN();
    //  cout<< np<< endl;
      double f = 0;
      double *x = Graph_XY_10->GetX();
      double *y = Graph_XY_10->GetY();
      double *ex = Graph_XY_10->GetEX();
      double *ey = Graph_XY_10->GetEY();
      for (int i=0;i<np;i++) {
         double u = x[i] - par[0];
         double v = y[i] - par[1];
         if( (u ==0 ) || (v==0) )continue;
         double dr = par[2] - std::sqrt(u*u+v*v);
         if(dr>40 && dr<20) continue; 
         double err2 = ((u*u)*(ex[i]*ex[i]) + (v*v)*(ey[i]*ey[i]))/(u*u+v*v);
         //cout << err2<< endl;
         f += dr*dr/err2;
         //cout <<f<< endl;
      }
      return f;
   };



ROOT::Math::Functor fcn(chi2Function,3);
    ROOT::Fit::Fitter  fitter;
    
   double pStart[3] = {center_x,center_y,est_rad};
   
   fitter.SetFCN(fcn, pStart);
   fitter.Config().ParSettings(0).Fix();
   fitter.Config().ParSettings(1).Fix(); 
    
   fitter.Config().ParSettings(0).SetName("x0");
   fitter.Config().ParSettings(1).SetName("y0");
   fitter.Config().ParSettings(2).SetName("R");
   
   
   // do the fit 
   bool ok = fitter.FitFCN();
   

   if (!ok) {
      Error("line3Dfit","Line3D Fit failed");
      return; 
   }  

     
   const ROOT::Fit::FitResult & result = fitter.Result();
   
   result.Print(std::cout);
   
  myrad_ = result.Parameters()[2];
    
  err_myrad_ = result.Errors()[2];  
    
  
   
   //Draw the circle on top of the points

  //return fitter->;

}//end

std::pair<float, float> computeCircumcenter(float x1, float y1, float x2, float y2, float x3, float y3) {
    float D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    
    float x_c = ((x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) + (x3*x3 + y3*y3) * (y1 - y2)) / D;
    float y_c = ((x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3) + (x3*x3 + y3*y3) * (x2 - x1)) / D;

    return std::make_pair(x_c, y_c);
}

void find_center(const char* filename = "20231017-194703/recodata.root"){
    
   
    float est_centx = 0.2;
    float est_centy  =0.2;
        



    // file to analyze 
    TFile *File;

    bool debug = false;
    File = TFile::Open(filename);
    TTree * _Tree= (TTree*) File->Get("recodata");

     //init branches
    float x_[60000];
    float y_[60000];
    float t_[60000];
    UShort_t n_;

        // set branch address
    cout << _Tree->GetEntries() << endl;
       // TGraphErrors* gr = new TGraphErrors();
    _Tree->SetBranchAddress("x",&x_);
    _Tree->SetBranchAddress("y",&y_);
    _Tree->SetBranchAddress("t",&t_);
    _Tree->SetBranchAddress("n",&n_);
    //tree entry loop
    //10 GeV tree Loop
    TH2D * hist_cent = new TH2D("hist_cent", ";Est_Cent; Entries", 200,-30,30,200,-30,30 );
    TH1D * rad_cent = new TH1D ("rad_cent",";Radius; Entries", 30, 0,100);
    TH1D * rad_cent2 = new TH1D ("rad_cent2",";Radius; Entries", 30, 0,100);
    TGraph * gr_cent =  new TGraph();

    for(int i=0; i<_Tree->GetEntries();i++){
        _Tree->GetEntry(i);
        //finding 3 points with in some distance
        int ncomb = 0;
        //cout<< n_<< endl;
       // if(n_ < 20) continue;
        for (int j =0; j<n_-2;j++){
           
           // cout<< j+1 << endl;

            for(int k =j+1 ; k< n_-1 ; k++){
            
                for(int l = k+1; l< n_ ;l++){
                    if(TMath::Abs(t_[j]) >10) continue;
                    if(TMath::Abs(t_[k]) >10) continue;
                    if(TMath::Abs(t_[l]) >10) continue;
                   
                    
                    
                  // float rad2 = TMath::Sqrt(); 
                   float x1 = x_[j];
                   float y1 = y_[j];
                   float x2 = x_[k];
                   float y2 = y_[k]; 
                   float x3 = x_[l];
                   float y3 = y_[l];
                    
                   float rad1 = TMath::Sqrt(x1*x1+y1*y1);
                   float rad2 = TMath::Sqrt(x2*x2+y2*y2);
                   float rad3 = TMath::Sqrt(x3*x3+y3*y3);  
                   if(TMath::Abs(rad1-rad2)>3 || TMath::Abs(rad2-rad3)>3 || TMath::Abs(rad1-rad3)>3 ) continue;  // making sure the points are part of possible circle
                   
                   float midx = 0.5 *(x1+x2);
                   float midy = 0.5 *(y1+y2);
                   float area =0.5 * TMath::Sqrt((x3-midx)*(x3-midx) +(y3-midy)*(y3-midy)) * TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)); 
                    if( area <0.01) continue;
                    
                   auto [centx, centy] = computeCircumcenter(x1, y1, x2, y2, x3, y3);
                   if(TMath::Abs(centx)>10||TMath::Abs(centy)>10) continue;
                       

                    
                 // float  centx = (x1+x2+x3)/3.0;
                 // float  centy = (y1+y2+y3)/3.0;
                  
                  float final_rad = TMath::Sqrt( (centx-x1)*(centx-x1) +(centy-y1)*(centy-y1));
                  float final_rad2 = TMath::Sqrt( (centx-x2)*(centx-x2) +(centy-y2)*(centy-y2));  
                                                
                                                
                                          
                    
                  hist_cent->Fill(centx,centy);  
                  rad_cent->Fill(final_rad);  
                  rad_cent2->Fill(final_rad2);    
                  //cout << rad1 <<"   " << rad2<< "    "<< rad3<< endl;  
                  //cout << "\n \n \n " << "center is " << Form("%f, %f",centx,centy)<< endl;  

                 gr_cent->SetPoint(ncomb, centx,centy);
                 ncomb++;    
                // cout<< ncomb<< endl;   
                }


            }


        }




    }//end entries

    TCanvas * can = new TCanvas("can","",1200,1200);
    can ->Divide(2,2);
    can->cd(1);
   // can->DrawFrame(-5,-5,5,5);
    //gr_cent->Draw("P");
    hist_cent->Draw("COLZ");
    can->Draw();
    can->cd(2);
    TH1D*  hist_cent_x = hist_cent->ProjectionX("hist_centx");
    hist_cent_x->Draw("e");
    TF1 * xgaus = new TF1("xgaus", "gaus" , -2,2);
    xgaus->SetParameters(10,0,1);
    xgaus->SetParLimits(2,0.1,1);
    hist_cent_x->Fit("xgaus","QR");
    can->cd(3);
    TH1D*  hist_cent_y = hist_cent->ProjectionY("hist_centy");
    hist_cent_y->Draw("e");
     TF1 * ygaus = new TF1("ygaus", "gaus" , -2,2);
    ygaus->SetParameters(10,0,1);
    hist_cent_y->Fit("ygaus","QR");
    
    can->cd(4);
    gPad->SetLogy();
    rad_cent->Draw();
    rad_cent2->SetLineColor(kRed);
    rad_cent2->SetLineStyle(10);
    rad_cent2->Draw("eSAME");
    
   // can->Draw();
    can->SaveAs("center.png");
   // can->cd(5);
   //gPad->DrawFrame(-5,-5,5,5);
    //gr_cent->Draw("");
    
center_x = xgaus-> GetParameter(1);
center_y = ygaus->GetParameter(1);
}


void find_radius(){
     const char *_filename = "20231017-190808/recodata.root";
     //20231017-194703 10 GeV  polarity Neg
    //20231017-190808 8GeV  polarity pos
    find_center(_filename);
    
    
    cout<< Form("%f , %f \n \n \n", center_x, center_y) <<endl;   
    


    // file to analyze 
    TFile *File = TFile::Open(_filename);
    if(!File){
     cout << "no file by name " << _filename << endl;
    } 
    
        
    bool debug = false;
   //
    TH2D* nvsr = new TH2D("nvsr",";r;number of hits per event",33,0,100,100,0,100);
    
    TH2D* tvsr = new TH2D("tvsr",";r;time",33,0,100,20,-30,30);
    
   //  TH2D* nvsr = new TH2D("nvsr",";r;number of hits per event",33,0,100,5,0,5);
    
    TTree * _Tree= (TTree*) File->Get("recodata");

    //TTree * Eight_Gev_Tree= (TTree*) En_8_File->Get("recodata");

     //init branches
    float x_[60000];
    float y_[60000];
    float t_[60000];
    UShort_t n_;

        // set branch address
    // cout << _Tree->GetEntries() << endl;
       // TGraphErrors* gr = new TGraphErrors();
    _Tree->SetBranchAddress("x",&x_);
    _Tree->SetBranchAddress("y",&y_);
    _Tree->SetBranchAddress("t",&t_);
    _Tree->SetBranchAddress("n",&n_);
    
    
    for(int i=0; i<_Tree->GetEntries(); i++){
        _Tree->GetEntry(i);
        //float  ev_rad = 0;
        int num_hits1 = 0;
        int num_hits2 = 0;
        int num_hits3 = 0;
        vector <float> inrad ;
        vector <float> midrad;
        vector <float> outrad;
        
        
        if(n_>=18){
            for(int j =0; j< n_;j++){
            if(t_[j]<-6 || t_[j]>20) continue;
            //float prev_rad =  ev_rad;
             
            float rad = TMath::Sqrt((x_[j]-center_x)*(x_[j]-center_x)+(y_[j]-center_y)*(y_[j]-center_y));
             
            if(rad<=52){
                num_hits1++;
                inrad.push_back(rad);
                tvsr->Fill(rad,t_[j]);
                
                
            }
            if(rad>52){
                if(t_[j]>6) continue;
                num_hits2++;
                outrad.push_back(rad);
                tvsr->Fill(rad,t_[j]);
                
                
            
            }     
                
          //  if((TMath::Abs(rad)-prev_rad)>3){
            
            //}
            //v_rad = rad;
            
         //   cout<< "radius"<< rad<< endl;
           // nvsr->Fill(rad,n_);
            
        }// vector loop
        }//n>18    
        
         if(n_<18){
            //cout << "####################"<< n_<< endl; 
            for(int j =0; j< n_;j++){
             if(t_[j]<-6 || t_[j]>20) continue;   
            //float prev_rad =  ev_rad;
            if(t_[j]>6) continue; 
            float rad = TMath::Sqrt((x_[j]-center_x)*(x_[j]-center_x)+(y_[j]-center_y)*(y_[j]-center_y));
            
            //if(rad>50 && rad<80) {
                
                
                num_hits3++;
                midrad.push_back(rad);
                tvsr->Fill(rad,t_[j]);
            
            //}     
                
          //  if((TMath::Abs(rad)-prev_rad)>3){
            
            //}
            //v_rad = rad;
            
         //   cout<< "radius"<< rad<< endl;
           // nvsr->Fill(rad,n_);
            
        }// vector loop
        }  //n<18
       
      for(auto v: inrad){
           //cout << "********" << v << endl;
          nvsr->Fill(v,num_hits1);
      } 
        
       for(auto v: outrad){
         // cout << "********" << v << endl;
          nvsr->Fill(v,num_hits2);
      }   

      for(auto v: midrad){
         // cout << "********" << v << endl;
          nvsr->Fill(v,num_hits3);
      }   
  
        
    }// entries
    TH1D* hr = nvsr->ProjectionX("hr");
    TCanvas* mcan = new TCanvas("mcan","",1800,600);
    mcan->cd();
    hr->GetXaxis()->SetRangeUser(65,80);
    
    hr->SetLineWidth(6);
    hr->SetLineStyle(8);
    hr->SetLineColor(32);
    hr->Draw("e");
   
    mcan->SaveAs("only_radius.png");
    TCanvas* radcan = new TCanvas("radcan","",1800,600);
    radcan->Divide(2,1);
    radcan->cd(1);
    nvsr->Draw("COLZ");
    radcan->cd(2);
    
    tvsr->Draw("COLZ");
    
    radcan->SaveAs("radius_nhits.png");  
    
    
    
}


void find_background(){
 const char *_filename = "20231017-190808/recodata.root";
  find_center(_filename);
   // file to analyze 
    TFile *File = TFile::Open(_filename);
    if(!File){
     cout << "no file by name " << _filename << endl;
     return;   
    } 
   // nbin = 160/2.5;
    TH2D* tvsr = new TH2D("tvsr",";r;time",75,0,150,60,-80,80);
    TH2D* tvsr_nogas = new TH2D("tvsr_nogas",";r;time",75,0,150,60,-80,80);
    TH2D* tvsr_gas = new TH2D("tvsr_nogas",";r;time",75,0,150,60,-80,80);
   // TH2D* xy = new TH2D("xy",";r;time",75,0,150,60,-80,80);
    TH1D * calc_rad = new TH1D("calc_rad", ";radius;", 30,20,100);
    
    TH1D * est_rad = new TH1D("est_rad", ";radius;", 30,20,100);
    
    TTree * _Tree= (TTree*) File->Get("recodata");
    
    //TTree * Eight_Gev_Tree= (TTree*) En_8_File->Get("recodata");

     //init branches
    float x_[60000];
    float y_[60000];
    float t_[60000];
    UShort_t n_;

        // set branch address
    // cout << _Tree->GetEntries() << endl;
       // TGraphErrors* gr = new TGraphErrors();
    _Tree->SetBranchAddress("x",&x_);
    _Tree->SetBranchAddress("y",&y_);
    _Tree->SetBranchAddress("t",&t_);
    _Tree->SetBranchAddress("n",&n_);
    
    for(int i=0; i<_Tree->GetEntries(); i++){
        _Tree->GetEntry(i);
        //float  ev_rad = 0;
        int num_hits_in = 0;
        int num_hits_mid = 0;
        int num_hits_out = 0;
        vector <float> inrad ;
        vector <float> midrad;
        vector <float> outrad;
        
          
          TH1D* hinrad = new TH1D("hinrad", ";radius;entries", 20 ,0,60);
          TH1D* hmidrad = new TH1D("hmidrad", ";radius;entries", 20 ,30,80);
          TH1D* houtrad = new TH1D("houtrad", ";radius;entries", 20 ,30,100); 
          TGraphErrors *gr_outrad = new TGraphErrors();
          TGraphErrors *gr_inrad = new TGraphErrors();
          TGraphErrors *gr_midrad = new TGraphErrors(); 
          for(int j =0; j< n_;j++){
            if(TMath::Abs(t_[j])>20) continue;
            float rad = TMath::Sqrt((x_[j]-center_x)*(x_[j]-center_x)+(y_[j]-center_y)*(y_[j]-center_y));
            
            tvsr->Fill(rad,t_[j]);
            est_rad->Fill(rad);  
              
           if(rad<60 && rad>30){
           num_hits_in++;
           hinrad->Fill(rad);    
           gr_inrad->SetPoint(num_hits_in,x_[j],y_[j]);
           gr_inrad->SetPointError(j,1.5,1.5);    
           } //if inrad 
              
          if(n_<20 && rad>55 && rad<70) {
             num_hits_mid++;  
            hmidrad->Fill(rad);
            gr_midrad->SetPoint(num_hits_mid,x_[j],y_[j]);
             gr_midrad->SetPointError(j,1.5,1.5);     
        }//if midrad
        
        if( rad<80&&rad>60 && n_>20) {
            cout << "found out radius " << endl;
             num_hits_out++;
            houtrad->Fill(rad);
            gr_outrad->SetPoint(num_hits_mid,x_[j],y_[j]);
             gr_outrad->SetPointError(j,1.5,1.5);   
            
        }//if outrad
            
              
        
        }// j hits
        
        
        
     TF1 * mgaus = new TF1("mgaus","gaus", 0, 100);
     mgaus->SetParameters(1,40,3);   
     //hinrad->Fit("mgaus","R");
     if(num_hits_in>=3){
     fitcircle(gr_inrad,40);  
     float inrad_ = myrad_;
     float err_inrad_ = err_myrad_;
     
     calc_rad->Fill(inrad_);    
     
     }
        
     if(num_hits_mid>=3){    
     fitcircle(gr_midrad,65);  
     //mgaus->SetParameters(1,65,3);
     //hmidrad->Fit("mgaus","R");   
     float midrad_ = myrad_;
     float err_midrad_ =err_myrad_;
     calc_rad->Fill(midrad_);    
         
     } 
        
     if(num_hits_out>=3){
     fitcircle(gr_outrad,70);     
     //mgaus->SetParameters(1,77,3);
     //houtrad->Fit("mgaus","R");       
     
     float outrad_ = myrad_;
     float err_outrad_ = err_myrad_;    
       
          calc_rad->Fill(outrad_);    
     }   
    
    
        
        
        
    
}// entry loop
TCanvas * can_1 = new TCanvas("can_1", "",600,600);
TCanvas * can_2 = new TCanvas("can_2", "",600,600); 
TCanvas * can_3 = new TCanvas("can_3", "",600,600);
    
can_1->cd();
calc_rad->SetLineWidth(4);    
est_rad->SetLineWidth(4); 
est_rad->SetLineColor(kRed);    


est_rad->Draw("");     
calc_rad->Draw("SAME"); 
can_1->SaveAs("radius_calculated_and_Est.png")  ;  

can_2->cd();

est_rad->SetLineWidth(4);    
est_rad->Draw(""); 
can_2->SaveAs("radius_estimated.png");    
    
} //find_background    