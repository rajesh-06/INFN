#include "iostream"
using namespace std;
float center_x = -999;
float center_y =-999;
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
    hist_cent_x->Fit("xgaus","R");
    can->cd(3);
    TH1D*  hist_cent_y = hist_cent->ProjectionY("hist_centy");
    hist_cent_y->Draw("e");
     TF1 * ygaus = new TF1("ygaus", "gaus" , -2,2);
    ygaus->SetParameters(10,0,1);
    hist_cent_y->Fit("ygaus","R");
    
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
