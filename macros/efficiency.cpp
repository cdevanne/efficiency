/*
In the code, N is for the number of real particles in the detector, k1 is the number of tracks in the 1st detector (REF) and k2 in the second one (TEST)

This code is not optimal and take some approximation and 'shortcut method' to write it faster. I just want to do a quick analysis, not a qualitative one. So it may be incorrect for MAX ~1 or 2% for the case k<5. And more for the case k=5
Anyway, even if I want to upgrade this, it's mandatory to do some approximation because we miss some information P(N), P(k=0), P(k>5)
*/
#include <numeric>
#include <TROOT.h> //gROOT
#include <TTimeStamp.h>
#include <TEfficiency.h>
#include <iomanip>
#include <TH1D.h>
#include <TFile.h>

#include "TGraphAsymmErrors.h"
#include <TGraphAsymmErrors.h>

#include "readModulesEfficiency.h"
#include "eventGenerator.h"

using namespace std;

void efficiency()
{
    // Seeding the random generator
    TTimeStamp currentTime; 
    long seed = currentTime.GetSec() * 1000 + currentTime.GetNanoSec() / 1000000;
    TRandom3 randomGenerator;
    randomGenerator.SetSeed(seed);
    cout << fixed << setprecision(2);

    //options
    bool switchStation = false;
    int nbEvents_for = 1e6;
    double MergingErrorRate = -0.5;
    
    TH1D* hMultiplicity = new TH1D("hMultiplicity", "Multiplicity Andrea&Katie", Nmax + 1, 0., Nmax + 1);
    TH1D* hMatch = new TH1D("hMatch", "Match Andrea&Katie", Nmax + 1, 0., Nmax + 1);

    TH1D* hMCMultiplicity = new TH1D("hMCMultiplicity", "MC Multiplicity", Nmax + 1, 0., Nmax + 1);
    TH1D* hMCMatch = new TH1D("hMCMatch", "MC Match", Nmax + 1, 0., Nmax + 1);

    TH1D* hMC2Multiplicity = new TH1D("hMC2Multiplicity", "MC2 Multiplicity Upgrade", Nmax + 1, 0., Nmax + 1);
    TH1D* hMC2Match = new TH1D("hMC2Match", "MC2 Match Upgrade", Nmax + 1, 0., Nmax + 1);

    TGraphAsymmErrors *graphRatioMatchMultiplicity = new TGraphAsymmErrors();
    TGraphAsymmErrors *graphRatioMCMatchMCMultiplicity = new TGraphAsymmErrors();
    TGraphAsymmErrors *graphRatioMC2MatchMC2Multiplicity = new TGraphAsymmErrors();



    // Opening root file
    vector<TGraphAsymmErrors> efficiencyModules;
    const char* fileName = "../datas/efficiency_allModules_CICbyCIC_window130um_run3232_3233_merged_1_yOffset.root";
    accessAndSaveEfficiencies(fileName, efficiencyModules);

    vector<double> proba_k = probabilities_k();
    vector<double> proba_N = probabilities_N(meanTrackEfficiency, true);

    int nbEvents_while = 0;
    
    double effMod = 0.;         //to have the mean efficiency of all modules
    double effMod_count = 0.;

    int nextk1 = 0;             //to simulate desynch stations
    int nextk2 = 0;

    // use 'while' or 'for' loop according to what you want to do

    // while (nbEvents_while < (total_trakcs_incoming)){ 
    for (int event = 0; event < nbEvents_for; ++event) {
        double random_N = randomGenerator.Rndm();
        int N = -1;
        while (random_N > 0) {
            N++;
            random_N -= proba_N[N];
        }
        if(N > 5) continue; //should not happend ?

        // Generating Tracks and position at targer
        vector<Tracks> tracks;
        int count_k1 = nextk1;
        int count_k2 = nextk2;
        int count_k1_upgrade = 0;
        int count_k2_upgrade = 0;

        nextk1 = 0;
        nextk2 = 0;
       
        generateTracks(N, tracks, randomGenerator, efficiencyModules, count_k1, count_k2, count_k1_upgrade, count_k2_upgrade, nextk1, nextk2, MergingErrorRate, effMod, effMod_count, switchStation);

        int nbMatch = 0;
        int nbMatch_upgrade = 0;

        matchTracks(N, tracks, nbMatch, nbMatch_upgrade);

        nbEvents_while += count_k1;

        //fill histo
        hMCMatch->Fill(count_k1, nbMatch);
        hMC2Match->Fill(count_k1_upgrade, nbMatch_upgrade);

        hMCMultiplicity->Fill(count_k1, count_k1);
        hMC2Multiplicity->Fill(count_k1_upgrade, count_k1_upgrade);

        if (nbEvents_while > 2e6) break;
    }

    cout << "Mean efficiency Modules: " <<  effMod/effMod_count << endl << endl;

    // Fill graoh
    for (int bin = 1; bin < hMCMatch->GetNbinsX(); ++bin) 
    {
        //value from Andre&Katie
        graphRatioMatchMultiplicity->SetPoint(bin-1, bin, matchEffNew[bin]);
        //value sim
        if(hMCMultiplicity->GetBinContent(bin+1) == 0) graphRatioMCMatchMCMultiplicity->SetPoint(bin-1, bin, 0);
            else graphRatioMCMatchMCMultiplicity->SetPoint(bin-1, bin, hMCMatch->GetBinContent(bin+1)/hMCMultiplicity->GetBinContent(bin+1));

        if(hMC2Multiplicity->GetBinContent(bin+1) == 0) graphRatioMC2MatchMC2Multiplicity->SetPoint(bin-1, bin, 0.); 
            else graphRatioMC2MatchMC2Multiplicity->SetPoint(bin-1, bin, hMC2Match->GetBinContent(bin+1)/hMC2Multiplicity->GetBinContent(bin+1)); 
        
    }


    // Sauvegarde dans le fichier ROOT
    TFile *outputFile = new TFile("output.root", "RECREATE");

    graphRatioMatchMultiplicity->Write("graphRatioMatchMultiplicity");
    graphRatioMCMatchMCMultiplicity->Write("graphRatioMCMatchMCMultiplicity");
    graphRatioMC2MatchMC2Multiplicity->Write("graphRatioMC2MatchMC2Multiplicity");
    
    // hMatch->Write();
    hMCMatch->Write();
    hMC2Match->Write();

    // hMultiplicity->Write();
    hMCMultiplicity->Write();
    hMC2Multiplicity->Write();

    outputFile->Close();
    
}
