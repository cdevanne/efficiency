#include "probaTools.h"
#include "andreKatieResults.h"

using namespace std;

int factorial(int n) 
{
    if (n <= 1) return 1;
    return n * factorial(n - 1);   
}

int N_choose_k(int N, int k) 
{
    if (N < k)  return 0;
    return factorial(N) / (factorial(k) * factorial(N - k));
}

double efficiencyTrackReco(double efficiencyModule) {
	return pow(efficiencyModule, 6) + 2*(pow(efficiencyModule, 5)*(1-efficiencyModule));
}

double efficiencyTrackRecoUpgrade(double efficiencyModule) {
	return pow(efficiencyModule, 6) + N_choose_k(6, 5)*(pow(efficiencyModule, 5)*(1-efficiencyModule));
}

double efficiencyTrackReco(const vector<double> &eff, int stationId) {
	int id = 0; //id of the module
	if (stationId == 1) id = 6;
	double proba = 1.;
	for (int i = 0; i < 6; ++i) { // 6 stubs
		proba *= eff[i+id]; 
		} 
	for (int idMissing = 2; idMissing < 4; ++idMissing) { // 5stubs
		double helper = 1;
		for (int i = 0; i < 6; ++i) {
			if (i==idMissing) helper *= (1-eff[i+id]);
			else helper *= eff[i+id];
		}
		proba += helper;
	}
	return proba;
}


double efficiencyTrackRecoUpgrade(const vector<double> &eff, int stationId) {
	int id = 0; //id of the module
	if (stationId == 1) id = 6;
	double proba = 1.;
	for (int i = 0; i < 6; ++i) { // 6 stubs
		proba *= eff[i+id]; 
		} 
	for (int idMissing = 0; idMissing < 6; ++idMissing) { // 5stubs
		double helper = 1;
		for (int i = 0; i < 6; ++i) {
			if (i==idMissing) helper *= (1-eff[i+id]);
			else helper *= eff[i+id];
		}
		proba += helper;
	}
	return proba;
}

double efficiencyTrackReco1S(const vector<double> &eff, int stationId) {
	/* 
	I consider Xa Ya Xa Ya U V Xb Yb Xb Yb, and the algo need at least 1Xa + 1Ya + 1Xb + 1Yb +(U or V)
	efficiency of Xa is sqrt(efficiency X 2S)
	*/
	int id = 0; //id of the module
	if (stationId == 1) id = 6;

	// 10 hits
	double proba = 1.;
	for (int i = 0; i < 12; ++i) { 
		if (i==5 || i==7) continue;
		proba *= sqrt(eff[i/2+id]); 
	} 

	// 9hits
	for (int idMissing = 0; idMissing < 12; ++idMissing) { 
		if (idMissing==5 || idMissing==7) continue;
		double helper = 1;
		for (int i = 0; i < 12; ++i) {
			if (i==5 || i==7) continue;
			if (i==idMissing) helper *= (1-sqrt(eff[(i/2)+id]));
			else helper *= sqrt(eff[(i/2)+id]);
		}
		proba += helper;
	}

	// 8 hits
	for (int id1 = 0; id1 < 12; ++id1) { 
		if (id1==5 || id1==7) continue;
		for (int id2 = id1+1; id2 < 12; ++id2) { 
			if (id2==5 || id2==7) continue;
			vector<int> check = {id1, id2};
			if (!checkIfAlgoIsWorking(check)) continue;
			double helper = 1;
			for (int i = 0; i < 12; ++i) {
				if (i==5 || i==7) continue;
				if (i==id1 || i==id2) helper *= (1-sqrt(eff[(i/2)+id]));
				else helper *= sqrt(eff[(i/2)+id]);
			}
		proba += helper;
		}
	}

	// 7 hits
	for (int id1 = 0; id1 < 12; ++id1) { 
	if (id1==5 || id1==7) continue;
	for (int id2 = id1+1; id2 < 12; ++id2) { 
	if (id2==5 || id2==7) continue;
	for (int id3 = id2+1; id3 < 12; ++id3) { 
	if (id3==5 || id3==7) continue;
		vector<int> check = {id1, id2, id3};
		if (!checkIfAlgoIsWorking(check)) continue;
		double helper = 1;
		for (int i = 0; i < 12; ++i) {
			if (i==5 || i==7) continue;
			if (i==id1 || i==id2 || i==id3) helper *= (1-sqrt(eff[(i/2)+id]));
			else helper *= sqrt(eff[(i/2)+id]);
		}
		proba += helper;
	}}}

	// 6 hits
	for (int id1 = 0; id1 < 12; ++id1) {
	if (id1==5 || id1==7) continue; 
	for (int id2 = id1+1; id2 < 12; ++id2) { 
	if (id2==5 || id2==7) continue;
	for (int id3 = id2+1; id3 < 12; ++id3) {
	if (id3==5 || id3==7) continue;
	for (int id4 = id3+1; id4 < 12; ++id4) {
	if (id4==5 || id4==7) continue;
		vector<int> check = {id1, id2, id3, id4};
		if (!checkIfAlgoIsWorking(check)) continue;
		double helper = 1;
		for (int i = 0; i < 12; ++i) {
			if (i==5 || i==7) continue;
			if (i==id1 || i==id2 || i==id3 || i==id4) helper *= (1-sqrt(eff[(i/2)+id]));
			else helper *= sqrt(eff[(i/2)+id]);
		}
		proba += helper;
	}}}}

	// 5 hits
	for (int id1 = 0; id1 < 12; ++id1) { 
	if (id1==5 || id1==7) continue;
	for (int id2 = id1+1; id2 < 12; ++id2) { 
	if (id2==5 || id2==7) continue;
	for (int id3 = id2+1; id3 < 12; ++id3) {
	if (id3==5 || id3==7) continue;
	for (int id4 = id3+1; id4 < 12; ++id4) {
	if (id4==5 || id4==7) continue;
	for (int id5 = id4+1; id5 < 12; ++id5) {
	if (id5==5 || id5==7) continue;
		vector<int> check = {id1, id2, id3, id4, id5};
		if (!checkIfAlgoIsWorking(check)) continue;
		double helper = 1;
		for (int i = 0; i < 12; ++i) {
			if (i==5 || i==7) continue;
			if (i==id1 || i==id2 || i==id3 || i==id4 || i==id5) helper *= (1-sqrt(eff[(i/2)+id]));
			else helper *= sqrt(eff[(i/2)+id]);
		}
		proba += helper;
	}}}}}

	return proba;
}

bool checkIfAlgoIsWorking(vector<int> modulesOFF)
{
	set<int> XaModules = {0, 1};
	set<int> XbModules = {8, 9};
	set<int> YaModules = {2, 3};
	set<int> YbModules = {10, 11};
	set<int> UVModules = {4, 6}; //5 and 7 are dead

	int Xa = 0;
	int Xb = 0;
	int Ya = 0;
	int Yb = 0;
	int UV = 0;

	for (int i = 0; i < (int)modulesOFF.size(); ++i)
	{
		if (XaModules.count(i)) Xa++;
		if (XbModules.count(i)) Xb++;
		if (YaModules.count(i)) Ya++;
		if (YbModules.count(i)) Yb++;
		if (UVModules.count(i)) UV++;
	}
	if (Xa==0 || Xb==0 || Ya==0 || Yb==0 || UV==0) {
		return false;
	}
	return true;
}

vector<double> probabilities_k()
{ //from datas, it's just missing p(k=0)
	vector<double> proba_k(Nmax+1,0); //k can be 0 (so +1)
	int size = incomingTracks.size();
	for (int i = 1; i < size; ++i)
	{
		proba_k[i] = ((double)incomingTracks[i]/i)/(double)total_trakcs_incoming;
	}
	return proba_k;
}

vector<double> probabilities_k_given_N(double efficiency, int N) 
{
    vector<double> proba_k_given_N(Nmax+1,0);
    double sum = 0;
    for (int k = 0; k <= N; ++k) {
        double p = (pow(efficiency, k) * pow(1.0 - efficiency, N - k))* N_choose_k(N, k);
	    proba_k_given_N[k] = p;
	    sum += p;
    }
    for (int i = 0; i <= Nmax; ++i) proba_k_given_N[i];
    return proba_k_given_N;
}

vector<double> probabilities_N(double efficiency , bool renormalisation )
{
	vector<double> proba_N(Nmax+1, 0);
	vector<double> proba_k = probabilities_k();
	double renorm = 0;
	for (int k = Nmax; k > 0; --k) {
	//each k help to determine p(N=k) and so can be use as the N associated parameter is some functions (ex: proba_k_given_N);
		vector<double> proba_k_given_N = probabilities_k_given_N(efficiency, k);
		proba_N[k] = proba_k[k];
		for (int N = k+1; N <= Nmax; ++N) {
			proba_N[k] -= proba_k_given_N[N] * proba_N[N]; 
		}
		proba_N[k] /= proba_k_given_N[k];
		if (renormalisation) renorm += proba_N[k];
	}
	if (renormalisation) {
		for (int i = 0; i <= Nmax; ++i) {
			proba_N[i]/= renorm;
		}
	}
	return proba_N;
}