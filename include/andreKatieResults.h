#ifndef ANDRE_KATIE_RESULTS_H
#define ANDRE_KATIE_RESULTS_H

#include <vector>
#include <numeric> //acumulate vector

using namespace std;

//Andrea & Katie datas


int k0events = 1; 
// this value doesn't matter because we don't considere events when k1 = 0. The p(k) probability change but at the end of the algo, the observable probabilitues still the same (p(k=0) is not observable). if we renormalise P(N) at the end it will always give the same value. 
// with this method I have great results but I am not sure I can do it (and I feel I can't). so I have to investivate a bit about that

const vector<int> incomingTracks = {k0events, 595782, 165354, 30078, 4144, 255}; 
const vector<int> outgoingTracks = {k0events, 599955, 159102, 26970, 3748, 160};
const vector<int> match = {0, 515224, 132045, 22118, 2704, 160};
const vector<double> matchEffNew = {0, 0.825, 0.765, 0.705, 0.605, 0.555};

const int total_trakcs_incoming = accumulate(incomingTracks.begin(), incomingTracks.end(), 0.0);
const int total_trakcs_outgoing = accumulate(outgoingTracks.begin(), outgoingTracks.end(), 0.0);

const double zijMax = 2.7;





//my Approximations:
const int Nmax = 5;

#endif

