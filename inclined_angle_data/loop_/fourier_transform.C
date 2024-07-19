#include <TGraphErrors.h>
#include <TVirtualFFT.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>

// Function to perform Fourier Transform on TGraphErrors and plot the result
void FourierTransform(TGraphErrors* graph) {
    // Get the number of points in the graph
    int nPoints = graph->GetN();

    // Get the y data points (assuming x data points are uniformly spaced by 28 minutes)
    double* y = graph->GetY();
    
    // FFT parameters
    double samplingInterval = 28 * 60; // 28 minutes in seconds
    int nFFT = nPoints;

    // Create TVirtualFFT object for R2C transform
    TVirtualFFT* fft = TVirtualFFT::FFT(1, &nFFT, "R2C EX K");

    // Set the input for the FFT
    fft->SetPoints(y);
    fft->Transform();

    // Get the transformed (Fourier) data
    double* re = new double[nFFT];
    double* im = new double[nFFT];
    fft->GetPointsComplex(re, im);

    // Calculate magnitudes of the Fourier coefficients
    double* magnitudes = new double[nFFT];
    double* frequencies = new double[nFFT];
    for (int i = 0; i < nFFT; ++i) {
        magnitudes[i] = sqrt(re[i] * re[i] + im[i] * im[i]);
        frequencies[i] = i / (samplingInterval * nFFT); // Frequency in Hz
    }

    // Create a graph for the magnitudes
    TGraph* magnitudeGraph = new TGraph(nFFT/2, frequencies, magnitudes);

    // Create a canvas to draw the result
    TCanvas* c1 = new TCanvas("c1", "Fourier Transform", 800, 600);

    // Draw the graph
    magnitudeGraph->SetTitle("Fourier Transform;Frequency (Hz);Magnitude");
    magnitudeGraph->Draw("AL");
    c1->SetLogy();
    c1->SetLogx();


    // Clean up
    delete[] re;
    delete[] im;
    delete[] magnitudes;
    delete[] frequencies;
}
