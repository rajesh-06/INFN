#include "cherenkov.h"
void calculate(const char *particle_="electron" ,const char *medium_ = "aerogel", float momentum_ = 10.){

    float length_aerogel = 379.5 ;
    float length_gas = 1312. ; // to do correct this
    cout << "making some montecarlo"<< endl;
    float sigma_radius = 2.0;
    const char *particle = particle_;
    const char *medium  = medium_;
    float give_mass = 0; 
    float give_index = 0; 
    float give_momentum = momentum_;
    float give_length = 0;
    float give_med_length = 0;

    if(std::strcmp(particle,"electron")==0) {
      give_mass =   myest.kElectronMass;

    }

    else if(std::strcmp(particle,"pion")==0) {
      give_mass =   myest.kPionMass;

    }

    else if(std::strcmp(particle,"proton")==0) {
      give_mass =   myest.kProtonMass;

    }


    else if(std::strcmp(particle,"muon")==0) {
      give_mass =   myest.kMuonMass;

    }

    else if(std::strcmp(particle,"kaon")==0) {
      give_mass =   myest.kKaonMass;

    }

    else {

        cout << "no valid particle defined "<< endl;
        cout << "use " << "kaon , muon , pion , proton , electron"<< endl;
        //TH2D * emptyhist = TH2D();
        //return ;
    }


    if(std::strcmp(medium,"aerogel")==0) {
      give_index =   1.021;
      give_length = length_aerogel; //mm
      give_med_length = 4.0; // mm  

    }

    else if(std::strcmp(medium,"gas")==0) {
      give_index =   1.0008;
      give_length = length_gas; //mm
      give_med_length = 0.5* 1273;  // mm // to do fix this 



    }
    else {
        cout << "wrong  medium"<< endl;
        //return ;//emptyhist;

    }

    Cherenkov_Est myest; // intializing the class    
    
    
    
    
    TH2D* hr_vs_photon = new TH2D("hr_vs_photon"," ;x;y",70,-100,100, 70, -100,100);
    TH2D* hxy = new TH2D("hxy"," ;x;y",70,-100,100, 70, -100,100);
    hxy->Reset();
    //hxy->SetName();
    float radius_med_ = myest.get_expected_radius(give_mass, give_momentum ,give_index, give_length);    
    float exp_photons_ = myest.get_expected_photons(give_index, give_mass, give_momentum);
    float exp_photons = exp_photons_ *give_med_length *0.3 *2;// 0.3 is photon detection efficiency of SiPM
    
    cout<<" for momentum of "<< give_momentum <<  " GeV/c and for "<< particle <<" expected photons in "<< medium<< " is "  << exp_photons << endl;
    cout<<" for momentum of "<< give_momentum <<  " GeV/c and for "<< particle <<" expected radius in  "<< medium << " is " << radius_med_ << endl;


    // Seed the random number generator
        std::srand(static_cast<unsigned int>(std::time(nullptr)));

        // Create a random number generator
        std::mt19937 generator(std::rand());

        // Define the distribution for random numbers between -100 and 100
       // std::uniform_int_distribution<int> distribution(-100, 100);
       std::uniform_real_distribution<float> distribution(-100.0f, 100.0f);

    std::uniform_real_distribution<> dis(0.0, 1.0);    
    for(int j =0;j < 10000; j++ ){

           for(int i=0; i< int(exp_photons);){

             //  double randomValue = dis(generator);

               float random_X = distribution(generator);
               float random_Y = distribution(generator);
               float    rand_rad = TMath::Sqrt(random_X*random_X + random_Y * random_Y);

               if(TMath::Abs(rand_rad- radius_med_)>3.0)continue;   
               //cout<< "filling histograms"<< endl ;
               hxy->Fill(random_X , random_Y);

               i++;   
           }//photon loop

       }//event loop

     
} 