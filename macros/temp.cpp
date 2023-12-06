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

#include "Utils.h"
#include "andreKatieResults.h"
#include "probaTools.h"
#include "tracks.h"

using namespace std;

void efficiency()
{
    cout << fixed << setprecision(2);
    bool plot = true; //if wanna plot results
    bool plot3 = true;
    bool plotNewAndrea = true;
    bool switchStation = false;
    bool connectAnything = false;
    int nbEvents_for = 1e6;

    TH1D* hMulitplicity = new TH1D("hMulitplicity", "Mulitplicity Andrea&Katie", Nmax+1, 0., Nmax+1);
    TH1D* hMCMulitplicity = new TH1D("hMCMulitplicity", "MC Mulitplicity", Nmax+1, 0., Nmax+1);
    TH1D* hMC2Mulitplicity = new TH1D("hMC2Mulitplicity", "MC Mulitplicity Upgrade", Nmax+1, 0., Nmax+1);

    TH1D* hMatch = new TH1D("hMatch", "Match Andrea&Katie", Nmax+1, 0., Nmax+1);
    TH1D* hMCMatch = new TH1D("hMCMatch", "MC Match", Nmax+1, 0., Nmax+1);
    TH1D* hMC2Match = new TH1D("hMC2Match", "MC_Match_reco_Anything", Nmax+1, 0., Nmax+1);

    TH1D* hEfficiency = new TH1D("hEfficiency", "hMatch / hMulitplicity", Nmax+1, 0., Nmax+1);
    TH1D* hMCEfficiency = new TH1D("hMCEfficiency", "hMCMatch / hMCMulitplicity", Nmax+1, 0., Nmax+1);
    TH1D* hMC2Efficiency = new TH1D("hMC2Efficiency", "hMC2Match / hMC2Mulitplicity", Nmax+1, 0., Nmax+1);

    TH1D* hMCnbMatchM1 = new TH1D("temphMCnbMatchM1", "Mulitplicity 1", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM2 = new TH1D("temphMCnbMatchM2", "Mulitplicity 2", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM3 = new TH1D("temphMCnbMatchM3", "Mulitplicity 3", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM4 = new TH1D("temphMCnbMatchM4", "Mulitplicity 4", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM5 = new TH1D("temphMCnbMatchM5", "Mulitplicity 5", Nmax+1, 0., Nmax+1);

    TH1D* hMCnbMatchM1Norm = new TH1D("hMCnbMatchM1", "Mulitplicity 1", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM2Norm = new TH1D("hMCnbMatchM2", "Mulitplicity 2", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM3Norm = new TH1D("hMCnbMatchM3", "Mulitplicity 3", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM4Norm = new TH1D("hMCnbMatchM4", "Mulitplicity 4", Nmax+1, 0., Nmax+1);
    TH1D* hMCnbMatchM5Norm = new TH1D("hMCnbMatchM5", "Mulitplicity 5", Nmax+1, 0., Nmax+1);
    
    

    double MergingErrorRate;

    const double zijMax2 = zijMax;
    // const double zijMax2 = 1.e6; // accept anything

    // Seeding the random generator
    TTimeStamp currentTime; 
    long seed = currentTime.GetSec() * 1000 + currentTime.GetNanoSec() / 1000000;
    TRandom3 randomGenerator;
    randomGenerator.SetSeed(seed);

    vector<TGraphAsymmErrors> efficiencyModules;

    // Opening root file
    const char* fileName = "../datas/efficiency_allModules_CICbyCIC_window130um_run3232_3233_merged_1_yOffset.root";


    accessAndSaveEfficiencies(fileName, efficiencyModules);
    cout << endl << "Number of Modules: " << efficiencyModules.size()/2 << endl; //1module=2CIC

    vector<double> proba_k = probabilities_k();
    double meanTrackEfficiency = 0.866;
    vector<double> proba_N = probabilities_N(meanTrackEfficiency, true);

    vector<int> MC_k1(Nmax+1, 0);
    vector<int> MC_k2(Nmax+1, 0);
    vector<int> MC_k1_upgrade(Nmax+1, 0);
    vector<int> MC_N(Nmax+1, 0);
    vector<int> MC_match(Nmax+1, 0);
    vector<int> MC_match_upgrade(Nmax+1, 0);

    // use 'while' or 'for' loop according to what you want to do, while is to try compare datas with Andrea and Katie
    double beamXmean = 9.2*0.01;
    double beamYmean = 72.9*0.01;
    double beamXsigma = 150.1*0.01; 
    double beamYsigma = 94.5*0.01; //strip * stripLegnght (cm) (approx) -> Neglectible anyway
    double beamXmeanError = 0.02;   
    double beamYmeanError = 0.02;
    double beamXsigmaError = .0;
    double beamYsigmaError = .0;
    // correction for specific case on boundaries (can just ruin an event). It was a test to see something but I don't use it
    vector<double> boundaryLeft = {-3.5, -2., -2., -2., -3., -1.9, -4., -1.7, -1.5, -1.5, -4., -1.9};
    vector<double> boundaryRight = {3.5, 2.8, 3., 3.5, 3., 2.6, 3.4, 2.5, 3., 2.7, 3.5, 2.4};
    vector<double> boundaryValue = {0.985, 0.975, 0.982, 0.98, 0.965, 0.982, 0.96, 0.965, 0.985, 0.985, 0.985, 0.975};
    //if don't want to use 
    for (int i = 0; i < 12; ++i)
    {
        boundaryLeft[i] = -100;
        boundaryRight[i] = 100;
    }
    

    int nbEvents_while = 0;
    

    int debug = 0;
    double effMod = 0.;
    double effMod_count = 0.;
    int nextk1 = 0;
    int nextk2 = 0;

    // while (nbEvents_while - MC_k1[0] < (total_trakcs_incoming)){ 
    for (int event = 0; event < nbEvents_for; ++event) {
        // cout << event << endl;
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
       
        for (int i = 0; i < N; ++i)
        {
            double xPos = randomGenerator.Gaus(beamXmean, beamXsigma);
            double yPos = randomGenerator.Gaus(beamYmean, beamYsigma);

            if (abs(xPos) > 4.5 || abs(yPos) > 4.5) {
                Tracks newTrack(xPos, yPos, 0, 0, i, false, false,false, false);
                tracks.push_back(newTrack);
                continue;
            }

            double uPos = xPos*sqrt(2.)/2. - yPos*sqrt(2.)/2.;
            double vPos = xPos*sqrt(2.)/2. + yPos*sqrt(2.)/2.;
            double xPosErr = randomGenerator.Gaus(beamXmeanError, beamXsigmaError);
            double yPosErr = randomGenerator.Gaus(beamYmeanError, beamYsigmaError);
            double random_k1 = randomGenerator.Rndm();
            double random_k2 = randomGenerator.Rndm();
            bool k1_bool = false;
            bool k2_bool = false;
            bool k1_bool_upgrade = false;
            bool k2_bool_upgrade = false;

            
            // Use real modulse efficiency for the associated position

            // vector <double> modulesEfficiencyLocal(12, 0.966);
            vector <double> modulesEfficiencyLocal(12, 0.);
            
            for (int i = 0; i < 12; ++i) //12 Modules
            {
                double Xrange;
                set<int> xModule = {0, 4, 6, 10};
                set<int> yModule = {1, 5, 7, 11};
                set<int> uModule = {2, 8};
                set<int> vModule = {3, 9};

                if (xModule.count(i) > 0) {
                    Xrange = xPos;
                    if (xPos < boundaryLeft[i] || xPos  > boundaryRight[i] ) modulesEfficiencyLocal[i] = boundaryValue[i];
                    else if (yPos > 0) modulesEfficiencyLocal[i] = efficiencyModules[2*i].Eval(Xrange);
                    else modulesEfficiencyLocal[i] = efficiencyModules[2*i+1].Eval(Xrange);
                }
                else if (yModule.count(i) > 0) {
                    Xrange = yPos;
                    if (yPos < boundaryLeft[i] || yPos  > boundaryRight[i]) modulesEfficiencyLocal[i] = boundaryValue[i];
                    else if (xPos > 0) modulesEfficiencyLocal[i] = efficiencyModules[2*i].Eval(Xrange);
                    else modulesEfficiencyLocal[i] = efficiencyModules[2*i+1].Eval(Xrange);
                }
                else if (uModule.count(i) > 0) {
                    Xrange = uPos;
                    if (uPos < boundaryLeft[i] || uPos  > boundaryRight[i]) modulesEfficiencyLocal[i] = boundaryValue[i];
                    else if (vPos > 0) modulesEfficiencyLocal[i] = efficiencyModules[2*i].Eval(Xrange);
                    else modulesEfficiencyLocal[i] = efficiencyModules[2*i+1].Eval(Xrange);
                }
                else if (vModule.count(i) > 0) {
                    Xrange = vPos;
                    if (vPos < boundaryLeft[i] || vPos > boundaryRight[i]) modulesEfficiencyLocal[i] = boundaryValue[i];
                    else if (uPos > 0) modulesEfficiencyLocal[i] = efficiencyModules[2*i].Eval(Xrange);
                    else modulesEfficiencyLocal[i] = efficiencyModules[2*i+1].Eval(Xrange);
                }
          
                effMod += modulesEfficiencyLocal[i];
                effMod_count +=1.;
            }
            int stationA = 1;
            int stationB = 0;
            if(switchStation == true)
            {
                stationA = 0;
                stationB = 1;
            }
            
            double efficiency_k1 = efficiencyTrackReco(modulesEfficiencyLocal, stationA);
            double efficiency_k2 = efficiencyTrackReco(modulesEfficiencyLocal, stationB);

            // reco anything (change z condition later)
            double efficiency_k1_upgrade = efficiencyTrackReco(modulesEfficiencyLocal, stationA);
            double efficiency_k2_upgrade = efficiencyTrackReco(modulesEfficiencyLocal, stationB);

            //Algo upgrade
            // double efficiency_k1_upgrade = efficiencyTrackRecoUpgrade(modulesEfficiencyLocal, stationA);
            // double efficiency_k2_upgrade = efficiencyTrackRecoUpgrade(modulesEfficiencyLocal, stationB);

            //1S Modules
            // double efficiency_k1_upgrade = efficiencyTrackReco1S(modulesEfficiencyLocal, stationA);
            // double efficiency_k2_upgrade = efficiencyTrackReco1S(modulesEfficiencyLocal, stationB);

            // cout << efficiency_k1 << endl << efficiency_k2 << endl << efficiency_k1_upgrade << endl << efficiency_k2_upgrade << endl << endl;
 
            
            double errorMerging = randomGenerator.Rndm();
            double errorMerging2 = randomGenerator.Rndm();
            MergingErrorRate = -0.14;


            if (random_k1 < efficiency_k1) {
                if (errorMerging < MergingErrorRate)
                {
                    nextk1++;
                }else{
                    k1_bool = true;
                    ++count_k1;    
                }
            } 

            if (random_k2 < efficiency_k2) { 
                if (errorMerging2 < MergingErrorRate)
                {
                    nextk2++;
                }else{
                    k2_bool = true;
                    ++count_k2;    
                }
            }
            
            if (random_k1 < efficiency_k1_upgrade) {
                k1_bool_upgrade = true;
                ++count_k1_upgrade;     
            } 
            if (random_k2 < efficiency_k2_upgrade) { 
                k2_bool_upgrade = true;
                ++count_k2_upgrade;
            } 


            

            Tracks newTrack(xPos, yPos, xPosErr, yPosErr, i, k1_bool, k2_bool, k1_bool_upgrade, k2_bool_upgrade);
            tracks.push_back(newTrack);
        }

        // Comparing Tracks and try match them. The stats is by event: P(all tracks match for a k1 value (REF))
        



        // Algo fairmuone
        int nbMatch = 0;
        for (int i = 0; i < N; ++i)
        {
            if (tracks[i].Getk1() && tracks[i].Getk2()) {
                nbMatch++;
                tracks[i].SetMatched();
            } 
        }
        // create possible wrong connexions for 1st algo
        for (int i = 0; i < N; ++i) {
            if (!tracks[i].GetMatched()) {
                for (int j = i+1; j < N; ++j) {
                    double zij = tracks[i].z(tracks[j]);
                    if (zij <= zijMax && !tracks[j].GetMatched()){
                        if (((tracks[i].Getk1() && !tracks[i].Getk2()) && (!tracks[j].Getk1() && tracks[j].Getk2())) || ((!tracks[i].Getk1() && tracks[i].Getk2()) && (tracks[j].Getk1() && !tracks[j].Getk2()))) {
                            nbMatch++;
                            // MC_match[count_k1] += 1;
                            tracks[i].SetMatched();
                            tracks[j].SetMatched();
                        }
                    }
                }
            }
        }
      
        // Upgrade Algo
        int nbMatch_upgrade = 0;
        for (int i = 0; i < N; ++i)
        {
            if (tracks[i].Getk1_upgrade() && tracks[i].Getk2_upgrade()) {
                nbMatch_upgrade++;
                tracks[i].SetMatched_upgrade();
            } 
        }
        // create possible wrong connexions for 2nd algo
        for (int i = 0; i < N; ++i) {
            if (!tracks[i].GetMatched_upgrade()) {
                for (int j = i+1; j < N; ++j) {
                    double zij = tracks[i].z(tracks[j]);
                    if (zij <= zijMax2 && !tracks[j].GetMatched_upgrade()){
                        if (((tracks[i].Getk1_upgrade() && !tracks[i].Getk2_upgrade()) && (!tracks[j].Getk1_upgrade() && tracks[j].Getk2_upgrade())) || ((!tracks[i].Getk1_upgrade() && tracks[i].Getk2_upgrade()) && (tracks[j].Getk1_upgrade() && !tracks[j].Getk2_upgrade()))) {
                            nbMatch_upgrade++;
                            // MC_match_upgrade[count_k1] += 1;
                            tracks[i].SetMatched_upgrade();
                            tracks[j].SetMatched_upgrade();
                        }
                    }
                }
            }
        }

        //fix 
        if (count_k1 > 5) count_k1 = 5;
        if (count_k2 > 5) count_k2 = 5;

        if (connectAnything == true && count_k1 <= count_k2)  
        {
            MC_match[count_k1] += count_k1;
        } else if(0 < nbMatch){
            MC_match[count_k1] += nbMatch;
        }
        if (0 < nbMatch_upgrade)  
        {
            MC_match_upgrade[count_k1_upgrade] += nbMatch_upgrade;
        }

        nbEvents_while += count_k1;
        MC_k1[count_k1] += count_k1;
        MC_k2[count_k2] += count_k2;
        MC_k1_upgrade[count_k1_upgrade] += count_k1_upgrade;

        //fill histo
        hMCMatch->Fill(count_k1, nbMatch);
        hMC2Match->Fill(count_k1_upgrade, nbMatch_upgrade);

        hMCMulitplicity->Fill(count_k1, count_k1);
        hMC2Mulitplicity->Fill(count_k1_upgrade, count_k1_upgrade);

        if (count_k1 == 1) {
            hMCnbMatchM1->Fill(nbMatch);
        }else if (count_k1 == 2) {
            hMCnbMatchM2->Fill(nbMatch);
        }else if (count_k1 == 3) {
            hMCnbMatchM3->Fill(nbMatch);
        }else if (count_k1 == 4) {
            hMCnbMatchM4->Fill(nbMatch);
        }else if (count_k1 == 5) {
            hMCnbMatchM5->Fill(nbMatch);
        }

        if (nbEvents_while > 2e6) break;
    }

    cout << "Mean efficiency Modules: " <<  effMod/effMod_count << endl << endl;

    cout << "k" << setw(9) << " MC_k1" << setw(9) << "k1" << setw(9) << "MC_k2" << setw(9) << "k2" << setw(9) << "MC_mat. match" << setw(9) << "k1_u" << endl;
    for (int k = 1; k <= Nmax; ++k) {
        cout << k << setw(9) << MC_k1[k] << setw(9) << incomingTracks[k] << setw(9) << MC_k2[k] << setw(9) <<  outgoingTracks[k] << setw(9) << MC_match[k]  << setw(9) << match[k] << setw(9) << MC_k1_upgrade[k] << endl;
    }
    cout << endl;

    //void plotVectors(const vector<double>& match, const vector<double>& MC_match) {

    // Créer un TGraph avec les données des vecteurs
    int size = match.size() - 1;
    double* x = new double[size];
    double* y = new double[size];
    double* y2 = new double[size];
    double* y3 = new double[size];
    double* y4 = new double[size];

    for (int i = 0; i < size; ++i) {
        x[i] = i+1; // Utiliser les indices comme valeurs x
        y[i] = (double)match[i+1]/(double)incomingTracks[i+1];
        y2[i] = (double)MC_match[i+1]/(double)MC_k1[i+1];
        y3[i] = (double)MC_match_upgrade[i+1]/(double)MC_k1_upgrade[i+1];
        y4[i] = matchEffNew[i+1];
    }

    if (plot)
    {

        TGraph* graphMatch = new TGraph(size, x, y);
        TGraph* graphMatchNewAndrea = new TGraph(size, x, y4);
        TGraph* graphMCMatch = new TGraph(size, x, y2);
        TGraph* graphMCMatch_upgrade = new TGraph(size, x, y3);

        TCanvas* canvas = new TCanvas("canvas", "Comparison", 1800, 1200);
        graphMatch->SetMarkerStyle(21);
        graphMatch->SetMarkerSize(2.0);
        graphMatch->SetMarkerColor(kBlue);
        graphMatch->SetLineColor(kBlue);
        graphMatch->SetLineWidth(2);

        graphMatchNewAndrea->SetMarkerStyle(21);
        graphMatchNewAndrea->SetMarkerSize(2.0);
        graphMatchNewAndrea->SetMarkerColor(kBlue);
        graphMatchNewAndrea->SetLineColor(kBlue);
        graphMatchNewAndrea->SetLineWidth(2);

        graphMCMatch->SetMarkerStyle(21); 
        graphMCMatch->SetMarkerSize(2.0);
        graphMCMatch->SetMarkerColor(kRed);
        graphMCMatch->SetLineColor(kRed);
        graphMCMatch->SetLineWidth(2);

        graphMCMatch_upgrade->SetMarkerStyle(21); 
        graphMCMatch_upgrade->SetMarkerSize(2.0);
        graphMCMatch_upgrade->SetMarkerColor(kGreen);
        graphMCMatch_upgrade->SetLineColor(kGreen);
        graphMCMatch_upgrade->SetLineWidth(2);

        TMultiGraph *mg = new TMultiGraph();
        if (plotNewAndrea == true) {
        mg->Add(graphMatchNewAndrea, "pl");
        } else mg->Add(graphMatch, "pl");
        
        mg->Add(graphMCMatch, "pl");
        if(plot3 == true) mg->Add(graphMCMatch_upgrade, "pl");
        mg->Draw("a");

        // legend
        TLegend* legend = new TLegend(0.1, 0.1, 0.3, 0.3);
        legend->AddEntry(graphMatch, "Match Andrea&Katie", "l");
        legend->AddEntry(graphMCMatch, "MC Match", "l");
        if (plot3 == true) legend->AddEntry(graphMCMatch_upgrade, "MC_Match_reco_Anything", "l");
        legend->Draw();
    }

    cout << "k  match   MC_mat. MC_match_upgrade" << endl;
    for (int i = 0; i < size; ++i) {
        cout << i+1 << "  " << y[i] << "    " << y2[i] << "     " << y3[i] << endl;
        /* 
        The case k=5 is better than k=4 just because of my approximation that N>5 = 0, so k=5 is not polluted by N>5 events
        */
    }
    cout << endl;

    for (int i = 1; i <= Nmax + 1; ++i)
    {
        hMulitplicity->Fill(i, incomingTracks[i]); 
        hMatch->Fill(i, match[i]);  
        //fill histo efficiency
        if (hMulitplicity->GetBinContent(i) != 0) {
            double efficiency = hMatch->GetBinContent(i) / hMulitplicity->GetBinContent(i);
            hEfficiency->SetBinContent(i, efficiency);
        } else hEfficiency->SetBinContent(i, 0);

        if (hMCMulitplicity->GetBinContent(i) != 0) {
            double efficiency = hMCMatch->GetBinContent(i) / hMCMulitplicity->GetBinContent(i);
            hMCEfficiency->SetBinContent(i, efficiency);
        } else hEfficiency->SetBinContent(i, 0);

        if (hMC2Mulitplicity->GetBinContent(i) != 0) {
            double efficiency = hMC2Match->GetBinContent(i) / hMC2Mulitplicity->GetBinContent(i);
            hMC2Efficiency->SetBinContent(i, efficiency);
        } else hEfficiency->SetBinContent(i, 0);        
    }
    
    for (int i = 0; i <= Nmax; ++i)
    {
        double total;

        total = hMCMatch->GetBinContent(1);
        if (total != 0) {
            double value = hMCnbMatchM1->GetBinContent(i) / total;
            hMCnbMatchM1Norm->SetBinContent(i, value);
        }else hMCnbMatchM1Norm->SetBinContent(i, 0);

        total = hMCMatch->GetBinContent(2);
        if (total != 0) {
            double value = hMCnbMatchM2->GetBinContent(i) / total;
            hMCnbMatchM2Norm->SetBinContent(i, value);
        }else hMCnbMatchM2Norm->SetBinContent(i, 0);

        total = hMCMatch->GetBinContent(3);
        if (total != 0) {
            double value = hMCnbMatchM3->GetBinContent(i) / total;
            hMCnbMatchM3Norm->SetBinContent(i, value);
        }else hMCnbMatchM3Norm->SetBinContent(i, 0);

        total = hMCMatch->GetBinContent(4);
        if (total != 0) {
            double value = hMCnbMatchM4->GetBinContent(i) / total;
            hMCnbMatchM4Norm->SetBinContent(i, value);
        }else hMCnbMatchM4Norm->SetBinContent(i, 0);

        total = hMCMatch->GetBinContent(5);
        if (total != 0) {
            double value = hMCnbMatchM5->GetBinContent(i) / total;
            hMCnbMatchM5Norm->SetBinContent(i, value);
        }else hMCnbMatchM5Norm->SetBinContent(i, 0);

    }
    
    



   
    
    TFile* outputFile = new TFile("output.root", "RECREATE");
    
    // hMatch->Write();
    // hMCMatch->Write();
    // hMC2Match->Write();

    // hMulitplicity->Write();
    // hMCMulitplicity->Write();
    // hMC2Mulitplicity->Write();

    // hEfficiency->Write();
    // hMCEfficiency->Write();
    // hMC2Efficiency->Write();

    hMCnbMatchM1Norm->Write();
    hMCnbMatchM2Norm->Write();
    hMCnbMatchM3Norm->Write();
    hMCnbMatchM4Norm->Write();
    hMCnbMatchM5Norm->Write();


    outputFile->Close();
    
}
