#include "tracks.h"

using namespace std;

// Builder
Tracks::Tracks(double x, double y, double errX, double errY, int id, bool k1, bool k2, bool k1_upgrade, bool k2_upgrade)
    : x(x), y(y), errX(errX), errY(errY), id(id), k1(k1), k2(k2), k1_upgrade(k1_upgrade), k2_upgrade(k2_upgrade), isMatched(false),  isMatched_upgrade(false) {}

double Tracks::distance(const Tracks& other) const {
    return sqrt(pow(other.x - x, 2) + pow(other.y - y, 2));
}

// Methods
double Tracks::absoluteError(const Tracks& other) const {
    return distance(other) * sqrt(pow(errX / (other.x - x), 2) + pow(errY / (other.y - y), 2));
}

double Tracks::z(const Tracks& other) const {
    return distance(other)/absoluteError(other);
}