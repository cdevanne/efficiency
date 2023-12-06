#ifndef UTILS_H
#define UTILS_H

#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TList.h>
#include <TObject.h>
#include <iostream>
#include <vector>
#include <TGraphAsymmErrors.h>

void accessAndSaveEfficiencies(const char* fileName, std::vector<TGraphAsymmErrors>& efficiencies);

#endif
