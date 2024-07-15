#pragma once


class radius{

public:
      
    
    radius(){
    //constructor
    }//const
    
    ~radius(
      //  rad_ delete;
       // time delete;
       // nhits delete;
    ){
    
    //destructor
    
    }//dest
    
    
    
   /* stuct Data{
        float  center_x;
        float  center_y;
        float  ring_id;
        
        vector float event; 
        
       
    
    
    }*/
    
     vector <float>rad_;
     vector <float> time;
     vector <float> nhits;
     vector <pair<float,float> > cent;
     vector < float> cent_y;
     vector < float> cent_x;
     vector<int> ring_tag; // if 1 then in 2 if mid and 3 if out
     int entrynum = 25342;
    
     
 
     
    
     float center_x = -999;
     float center_y =-999;
     TGraph * gr_cent =  new TGraph();
     virtual void find_center(const char* filename); 
     virtual void find_radius(const char* filename);
     
    void plot_circle(std::array<float, 3> parameters, int line_color = kBlack, int line_style = kSolid, int line_width = 1)
        {
          auto result = new TEllipse(parameters[0], parameters[1], parameters[2]);
          result->SetFillStyle(0);
          result->SetLineColor(line_color);
          result->SetLineStyle(line_style);
          result->SetLineWidth(line_width);
          result->DrawEllipse(parameters[0], parameters[1], parameters[2], 0, 0, 360, 0, "same");
        }
    
    
    virtual void beautify(TH1D* object,int col, int style){    
    object->SetLineWidth(3);
    object->SetLineColor(col);
    //object->SetOptStat(0);
    object->SetMarkerStyle(style);     
    
    }
    
    
    virtual void beautify(TLine* object,int col, int style){    
    object->SetLineWidth(3);
    object->SetLineColor(col);
    //object->SetOptStat(0);
    object->SetLineStyle(style);     
    
    }
    
    virtual void beautify(TArc* object,int col, int style){    
    object->SetLineWidth(3);
    object->SetLineColor(col);
    //object->SetOptStat(0);
    object->SetLineStyle(style);     
    
    }
    
    virtual void SetBinWidth(float bwidth){
        bin_width = bwidth;
    }
    // to divide the canvas
    virtual void Set_Divide_Canvas(bool dodivide, int nrow, int ncol){
         divide_canvas = dodivide;
         can_row = nrow;
         can_column = ncol;
     }
    
    virtual void Set_Time_Cut(float Time_Cut_Min,float Time_Cut_Max){
         tcut_min = Time_Cut_Min;
         tcut_max = Time_Cut_Max;
      }
    
    {
        outfilename = _outfilename; 
    }
    
    
   /* void WriteFile(auto file){
        
        new_file->cd();
         file->Write();
       
        new_file->Close();
    }*/
    
    
    
    TH2D *plot_histo2d(const char* histname, vector <float> x , vector <float> y){
      //TCanvas * can = new TCanvas("can", "", 1200,1200);
      
      if(x.empty() || y.empty()) 
      {  cout << "size of x  is " << x.size()<< endl;
         cout << "size of y  is " << y.size()<< endl;
          
          return 0x0;}
      // Use std::max_element to find the maximum element
      auto xmax_ = std::max_element(x.begin(), x.end());
      auto  xmin_ = std::min_element(x.begin(), x.end());
      
      
      auto ymax_ = std::max_element(y.begin(), y.end());
      auto  ymin_ = std::min_element(y.begin(), y.end());
     
      
      float xmax= *xmax_ ;
      float xmin =*xmin_;
      float ymax =*ymax_;
      float ymin =*ymin_;  
      int nbin_x = (int) TMath::Abs(xmax-xmin)/bin_width;  
      int nbin_y  = (int) TMath::Abs(ymax-ymin)/bin_width;  
        
        
      TH2D* hist_2 = new TH2D(histname,"",nbin_x,xmin,xmax,nbin_y,ymin,ymax);
          
      for(int i=0; i < x.size();i++){
      
      hist_2->Fill(x[i],y[i]);
      }//end loop
      //can->cd();
      //hist_2->Draw(option);
      //can->SaveAs(Form("%s.png",histname));   
      return hist_2;
     
       
    } // 
    
    
    
     TH1D*  plot_histo1d(const char* histname, vector <float> x ){
     
         TCanvas * can = new TCanvas(histname, "", 1200,1200);

          if(x.empty()) 
          {  cout << "size of x  is " << x.size()<< endl;


              return 0x0;}
          // Use std::max_element to find the maximum element
          auto xmax_ = std::max_element(x.begin(), x.end());
          auto  xmin_ = std::min_element(x.begin(), x.end());



          float xmax= *xmax_ ;
          float xmin =*xmin_;
          int nbin_x = (int) TMath::Abs(xmax-xmin)/3.0;  


          TH1D* hist_1 = new TH1D(histname,"",nbin_x,xmin,xmax);

          for(int i=0; i < x.size();i++){

              hist_1->Fill(x[i]);
          }//end loop
         
         return hist_1;    
    
    } 
    
    TCanvas * BuildCanvas(const char* name)
    {
        
    TCanvas * can = new TCanvas(name, "", 1200, 1200);
    if(divide_canvas){
    can->Divide(can_row,can_column);
   
    }    
        
    
    
    return can; 
    }
    

    
private:
    
    float x_[60000];
    float y_[60000];
    float t_[60000];
    UShort_t n_;
    int canvas_x =1200;
    int canvas_y =1200;
    bool divide_canvas = false;
    int can_row = 2;
    int can_column = 2;
    float bin_width =3.0;
    
    const char* outfilename ="default.root";
    
   // TFile * new_file =  TFile::Open(outfilename,"UPDATE");
    
    float tcut_max = 40;// selection on timing this is +20 to -20
    float tcut_min = -40;
    
    
    
    //compute circumcenter
    std::pair<float, float> computeCircumcenter(float x1, float y1, float x2, float y2, float x3, float y3) {
        float D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
        float x_c = ((x1*x1 + y1*y1) * (y2 - y3) + (x2*x2 + y2*y2) * (y3 - y1) + (x3*x3 + y3*y3) * (y1 - y2)) / D;
        float y_c = ((x1*x1 + y1*y1) * (x3 - x2) + (x2*x2 + y2*y2) * (x1 - x3) + (x3*x3 + y3*y3) * (x2 - x1)) / D;
        return std::make_pair(x_c, y_c);
            }//end compute circumcenter
    //compute area
   float compute_area (float x1, float y1,float x2,float y2,float x3,float y3){
          float midx = 0.5 *(x1+x2);
          float midy = 0.5 *(y1+y2);
          float area =0.5 * TMath::Sqrt((x3-midx)*(x3-midx) +(y3-midy)*(y3-midy)) * TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
          return area;
        }// to compute area
    
                
  // class to read the tree
   TTree *read_recodata(const char* filename){
     TFile  *File;
     File = TFile::Open(filename);
     TTree * _Tree= (TTree*) File->Get("recodata");
     _Tree->SetBranchAddress("x",&x_);
     _Tree->SetBranchAddress("y",&y_);
     _Tree->SetBranchAddress("t",&t_);
     _Tree->SetBranchAddress("n",&n_);
       
    return _Tree;   
       

   
   }


        












};





void radius::find_center(const char* filename){

   
    auto _Tree =(TTree*) read_recodata(filename);
     
    
    // set branch address
    cout << _Tree->GetEntries() << endl;
    
    
    TH2D * hist_cent = new TH2D("hist_cent", ";Est_Cent; Entries", 60,-100,100,60,-100,100 );
    TH1D * rad_cent = new TH1D ("rad_cent",";Radius; Entries", 50, 0,150);
    TH1D * rad_cent2 = new TH1D ("rad_cent2",";Radius; Entries", 50, 0,150);
    
    TH2D * hist_xy_single_ =  new TH2D("hist_xy_single_", ";; Entries", 70,-100,100,70,-100,100 );
    
    TH2D * hist_xy_single_ring =  new TH2D("hist_xy_single_ring", ";; Entries", 70,-100,100,70,-100,100);
    
    
    
    
    
    
    
    for(int i=0; i<_Tree->GetEntries();i++){
        
        bool has_bkg_ring =false; 
        _Tree->GetEntry(i);
        int ncomb = 0;
        vector<float> possible_center;
        //float centx_ =0.;
        //float centy_ =0.
        
        //if(n_<=4 ||n_>=10) continue;
        //if(n_<=4 ) continue;
       
        for (int j =0; j<n_-2;j++){
           
           // cout<< j+1 << endl;
           
         
            
            
          if( i == entrynum){
           
               hist_xy_single_->Fill(x_[j],y_[j]);
           } 

            for(int k =j+1 ; k< n_-1 ; k++){
            
                for(int l = k+1; l< n_ ;l++){
                    if(t_[j] >tcut_max) continue;
                    if(t_[k] >tcut_max) continue;
                    if(t_[l] >tcut_max) continue;
                    
                    if(t_[j] <tcut_min) continue;
                    if(t_[k] <tcut_min) continue;
                    if(t_[l] <tcut_min) continue;
                    
                   
                    
                    
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
                   
                    
                    
                   
                    // making sure the points are part of possible circle
                    
                   // if(TMath::Abs(rad1-rad2)>3 || TMath::Abs(rad2-rad3)>3 || TMath::Abs(rad1-rad3)>3 ) continue; 
                    
                   //   if(TMath::Abs(rad1-rad2)<3 && TMath::Abs(rad2-rad3)<3 && TMath::Abs(rad1-rad3)<3 ) continue; 
                    
                   //computing area so that we can put some threshold on area
                    auto area = compute_area(x1,y1,x2,y2,x3,y3);
                    if( area <0.1) continue;
                    
                   float dist1 = TMath::Sqrt((x1-x2)*(x1-x2)+ (y1-y2)*(y1-y2));
                   float dist2=  TMath::Sqrt((x2-x3)*(x2-x3)+ (y2-y3)*(y2-y3));
                   float dist3 = TMath::Sqrt((x1-x3)*(x1-x3)+ (y1-y3)*(y1-y3));
                    
                  //  if(dist1< ||dist2<20 ||dist3<20)continue;
                    
                    
                   auto [centx, centy] = computeCircumcenter(x1, y1, x2, y2, x3, y3);
                  if( TMath::Abs(centx)>100 ||  TMath::Abs(centy)>100 ) continue; 
                  if(TMath::Abs(centx)<10&&TMath::Abs(centy)<10) continue;
                 // if(TMath::Abs(centx)>10||TMath::Abs(centy)>10) continue;  
                       
                  float final_rad = TMath::Sqrt( (centx-x1)*(centx-x1) +(centy-y1)*(centy-y1));
                  float final_rad2 = rad1;//TMath::Sqrt( (centx-x2)*(centx-x2) +(centy-y2)*(centy-y2));
                  //float weight = 0;
                  float n_rad = 0; 
                    
                    
                //  single_fill = true;  
                 /* for(int m =0; m< n_; m++){
                      if(m!=j && m!=k && m!=l){
                          
                          float rad_all = TMath::Sqrt((centx- x_[m])*(centx- x_[m])+(centy- y_[m])*(centy- y_[m]));
                          rad_cent->Fill(rad_all,TMath::Gaus(rad_all,final_rad,3.));
                          float weight = TMath::Gaus(rad_all,final_rad,3.); 
                          
                          if(weight>0.8 && rad_all <25  && rad_all>15){
                                n_rad++;
                               if(n_<10) {
                              
                                  cout << "entry number with ring is "<< i << endl;
                          
                              }// end if n<10   
                             // cout << "entry number with ring is "<< i << endl;
                          
                          } 
                              
                          
                          else continue;    
                      
                         // hist_newring->Fill(x_[m],x_[m]); 
                          //rad_cent2->Fill(final_rad2,TMath::Gaus(rad_all,final_rad2,3.) );  
                          
                      }// compute center of all 3 points except 3 points used for the center finding 
                      
                      */
                    
                    
                    
                      //if(n_rad>3){
                          has_bkg_ring = true;
                          
                          
                          if( i == entrynum ){
                          
                               
                             // single_fill = true;
                            // cout << "filling"<< endl;
                              float rad_all_ = TMath::Sqrt((centx- x_[m])*(centx- x_[m])+(centy- y_[m])*(centy- y_[m]));
                              //if(TMath::Abs(centx-17)<3 && TMath::Abs(centy+46)<3 && weight>0.7){
                              cout<< n_ << " number of hits  " << n_rad << "  centx "<< centx << "   cent y" << centy << "   "<< rad_all_ <<endl;
                             // cout << "entry number with ring is "<< i << endl;
                              hist_xy_single_ring->Fill(x_[m],y_[m]);
                              //}
                              
                              
                              }  //only for fix entry
                  
                      //}
                      
                      
                  }//computing radius for all the points   
                    
                  
                  
                  std::pair centpair = std::make_pair(centx,centy) ;                             
                  
                  
                  cent_x.push_back(centx);
                  cent_y.push_back(centy);  
                  cent.push_back(centpair);                        
                    
                  hist_cent->Fill(centx,centy);  
                    
                  rad_cent2->Fill(final_rad2);    
                  //cout << rad1 <<"   " << rad2<< "    "<< rad3<< endl;  
                  //cout << "\n \n \n " << "center is " << Form("%f, %f",centx,centy)<< endl;  

                 gr_cent->SetPoint(ncomb, centx,centy);
                 ncomb++;    
                
                } // end final loop


            } //end second loop


        }// end third loop


     
                     

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
    hist_cent_x->Fit("xgaus","NQR");
    can->cd(3);
    TH1D*  hist_cent_y = hist_cent->ProjectionY("hist_centy");
    hist_cent_y->Draw("e");
     TF1 * ygaus = new TF1("ygaus", "gaus" , -2,2);
    ygaus->SetParameters(10,0,1);
    hist_cent_y->Fit("ygaus","NQR");
    
    can->cd(4);
    gPad->SetLogy();
    rad_cent->Draw();
    rad_cent2->SetLineColor(kRed);
    rad_cent2->SetLineStyle(10);
    rad_cent2->Draw("eSAME");
    
    TCanvas * canrad = new TCanvas("canrad","",1200,1200);
    canrad->cd();
   
    hist_xy_single_->SetMarkerSize(3);
    hist_xy_single_->SetMarkerStyle(20);
    hist_xy_single_ring->SetMarkerColor(kBlack);
    hist_xy_single_->Draw("SCAT");
    
    
    plot_circle({ 64.7297,-31.4095,17.7613});
    
    hist_xy_single_ring->SetMarkerSize(2);
    hist_xy_single_ring->SetMarkerStyle(2);
    
    hist_xy_single_ring->SetMarkerColor(kRed);
    
    //hist_xy_single_ring->Draw("SCAT SAME");
    
    canrad->SaveAs(Form("%s_%i-xy.png",filename,entrynum));
    
    
   // can->Draw();
    can->SaveAs(Form("%s-center.png",filename));
   // can->cd(5);
   //gPad->DrawFrame(-5,-5,5,5);
    //gr_cent->Draw("");
    
    
    
center_x = xgaus-> GetParameter(1);
center_y = ygaus->GetParameter(1);



}
//}//find_center

void radius::find_radius(const char *_filename){

    
    find_center(_filename);
    
    auto _Tree =(TTree*) read_recodata(_filename);
    
    
    cout<< Form(" ( center x, center y ) = ( %f , %f ) \n \n \n", center_x, center_y) <<endl;  
    
     for(int i=0; i<_Tree->GetEntries();i++){
        _Tree->GetEntry(i);
         for(int j =0; j< n_;j++){
         
         if(t_[j] >tcut_max) continue;
       
         if(t_[j] <tcut_min) continue;
        
         
        float rad = TMath::Sqrt((x_[j]-center_x)*(x_[j]-center_x)+(y_[j]-center_y)*(y_[j]-center_y));
        rad_.push_back(rad);
        time.push_back(t_[j]);
         }//j     
         
     }//tree entries
    
    
    
    

}//end find_radius