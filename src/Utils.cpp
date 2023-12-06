#include "Utils.h"

using namespace std;

void accessAndSaveEfficiencies(const char* fileName, vector<TGraphAsymmErrors>& efficiencies) 
{
    TFile *file = new TFile(fileName);

    if (!file || file->IsZombie()) {
        cerr << "Error: Could not open file " << fileName << endl;
        return;
    }

    TCanvas *canvasCLlim = dynamic_cast<TCanvas*>(file->Get("c_efficiency_CLlim"));
    if (!canvasCLlim) {
        cerr << "Error: Could not find TCanvas c_efficiency_CLlim in the file." << endl;
        file->Close();
        return;
    }

    TList *listOfPrimitives = canvasCLlim->GetListOfPrimitives();
    TIter next(listOfPrimitives);
    TObject *obj = nullptr;

    while ((obj = next())) {
        if (obj->IsA() == TPad::Class()) {
            TPad *pad = dynamic_cast<TPad*>(obj);
            TList *listOfPrimitives2 = pad->GetListOfPrimitives();
            TIter next2(listOfPrimitives2);
            TObject *obj2 = nullptr;

            while ((obj2 = next2())) {
                if (obj2->IsA() == TGraphAsymmErrors::Class()) {
                    TGraphAsymmErrors *graph = dynamic_cast<TGraphAsymmErrors*>(obj2);
                    TGraphAsymmErrors copy(*graph);
                    efficiencies.push_back(copy);
                }
            }
        }
    }
    

    file->Close();
}