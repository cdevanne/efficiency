#ifndef PROBA_TOOLS_H
#define PROBA_TOOLS_H

#include <vector>
#include <iostream>
#include <TRandom3.h>
#include <cmath>
#include <set>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>

int factorial(int n);
int N_choose_k(int N, int k);

// mean
double efficiencyTrackReco(double efficiencyModule);
double efficiencyTrackRecoUpgrade(double efficiencyModule);
// using each modules efficiency, more accurate with datas
double efficiencyTrackReco(const std::vector<double> &eff, int stationId);
double efficiencyTrackRecoUpgrade(const std::vector<double> &eff, int stationId);
double efficiencyTrackReco1S(const std::vector<double> &eff, int stationId);

bool checkIfAlgoIsWorking(std::vector<int> modulesOFF);

std::vector<double> probabilities_k();
std::vector<double> probabilities_k_given_N(double efficiency, int N);
std::vector<double> probabilities_N(double efficiency, bool renormalisation = false);


#endif