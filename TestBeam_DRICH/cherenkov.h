#ifndef CHERENKOV_H
#define CHERENKOV_H

#include <random>
#include <iostream>


using namespace std;

class Cherenkov_Est {
private:
    
    // Private utility functions
    float get_relativistic_beta(float gamma) {
        return TMath::Sqrt(1 - gamma * gamma);
    }

    float get_relativistic_beta(float mass, float momentum) {
        return 1. / sqrt(1 + (mass * mass) / (momentum * momentum));
    }

    float get_relativistic_gamma(float beta) {
        return 1. / (TMath::Sqrt(1 - beta * beta));
    }

    float get_relativistic_gamma(float mass, float momentum) {
        return get_relativistic_gamma(get_relativistic_beta(mass, momentum));
    }

    float get_theta_cherenkov(float beta, float ref_index) {
        return TMath::ACos(1. / (beta * ref_index));
    }

    float get_theta_atan(float radius, float arm_length) {
        return TMath::ATan(radius / arm_length);
    }

    

public:
    // Constructor
    Cherenkov_Est() {}

    // Destructor
    ~Cherenkov_Est() {}
    const float alhpa_EM = 1. / 137.055;   //[n.u.]
    const float kElectronMass = 0.0005110; // [GeV]
    const float kMuonMass = 0.1057;        // [GeV]
    const float kPionMass = 0.1396;        // [GeV]
    const float kKaonMass = 0.4937;        // [GeV]
    const float kProtonMass = 0.9383;      // [GeV]

    
    float get_expected_radius(float mass, float momentum, float ref_index, float arm_length) {
        
        
        return TMath::Tan(get_theta_cherenkov(get_relativistic_beta(mass, momentum), ref_index)) * arm_length;
    }

    float get_relativistic_beta_cherenkov_threshold(float ref_index) {
        return 1. / ref_index;
    }

    float get_momentum_cherenkov_threshold(float ref_index, float mass) {
        return mass * get_relativistic_beta_cherenkov_threshold(ref_index) * get_relativistic_gamma(get_relativistic_beta_cherenkov_threshold(ref_index));
    }

    float get_mass_cherenkov_threshold(float ref_index, float momentum) {
        return momentum / (get_relativistic_beta_cherenkov_threshold(ref_index) * get_relativistic_gamma(get_relativistic_beta_cherenkov_threshold(ref_index)));
    }

    float get_expected_photons(float ref_index, float mass, float momentum, int z_charge = 1, float lambda_low = 350.e-6, float lambda_high = 650.e-6) {
        return ((2 * TMath::Pi() * z_charge * z_charge * alhpa_EM) * ((1. / lambda_low) - (1. / lambda_high)) * (1. - 1. / (get_relativistic_beta(mass, momentum) * get_relativistic_beta(mass, momentum) * ref_index * ref_index)));
    }

    // Public member functions
    void show_expected_photons(float ref_index, float momentum) {
        cout << "electron: " << get_expected_photons(ref_index, kElectronMass, momentum) << endl;
        cout << "muon: " << get_expected_photons(ref_index, kMuonMass, momentum) << endl;
        cout << "pion: " << get_expected_photons(ref_index, kPionMass, momentum) << endl;
        cout << "kaon: " << get_expected_photons(ref_index, kKaonMass, momentum) << endl;
        cout << "proton: " << get_expected_photons(ref_index, kProtonMass, momentum) << endl;
    }

    void show_thresholds_per_momentum(float ref_index, float momentum) {
        cout << "electron: " << (get_momentum_cherenkov_threshold(ref_index, kElectronMass) < momentum) << endl;
        cout << "muon: " << (get_momentum_cherenkov_threshold(ref_index, kMuonMass) < momentum) << endl;
        cout << "pion: " << (get_momentum_cherenkov_threshold(ref_index, kPionMass) < momentum) << endl;
        cout << "kaon: " << (get_momentum_cherenkov_threshold(ref_index, kKaonMass) < momentum) << endl;
        cout << "proton: " << (get_momentum_cherenkov_threshold(ref_index, kProtonMass) < momentum) << endl;
    }

    void show_expected_radii(float ref_index, float momentum, float arm_length) {
        cout << "electron: " << get_expected_radius(kElectronMass, momentum, ref_index, arm_length) << endl;
        cout << "muon: " << get_expected_radius(kMuonMass, momentum, ref_index, arm_length) << endl;
        cout << "pion: " << get_expected_radius(kPionMass, momentum, ref_index, arm_length) << endl;
        cout << "kaon: " << get_expected_radius(kKaonMass, momentum, ref_index, arm_length) << endl;
        cout << "proton: " << get_expected_radius(kProtonMass, momentum, ref_index, arm_length) << endl;
    }

    void plot_prediction(float ref_index, float momentum, float mass, float distance, std::array<float, 2> xy = {0., 0.}) {
        // Placeholder function, you can implement the actual plotting logic here
        // plot_circle({xy[0], xy[1], get_expected_radius(mass, momentum, ref_index, distance)}, kRed, kDashed, 2);
        cout << "Plotting prediction..." << endl;
    }
};


#endif // CHERENKOV_H
/*int main() {
    // Example usage of Cherenkov_Est class
    Cherenkov_Est cherenkov_estimator;

    cherenkov_estimator.show_expected_photons(1.5, 100); // Example parameters
    cherenkov_estimator.show_thresholds_per_momentum(1.5, 100); // Example parameters
    cherenkov_estimator.show_expected_radii(1.5, 100, 500); // Example parameters

    // Example of plotting prediction (not implemented here)
    cherenkov_estimator.plot_prediction(1.5, 100, 0.1057, 500, {10., 20.}); // Example parameters

    return 0;
}*/

