#ifndef TRACKS_H
#define TRACKS_H

#include <vector>
#include <cmath>

class Tracks {
public:
    // Constructeur
    Tracks(double x, double y, double errX, double errY, int id, bool k1, bool k2, bool k1_upgrade, bool k2_upgrade);

    // Accesseurs
    double GetX() const {return x;};
    double GetY() const {return y;};
    double GetErrX() const {return errX;};
    double GetErrY() const {return errY;};
    int GetId() const {return id;};
    bool Getk1() const {return k1;};
    bool Getk2() const {return k2;};
    bool Getk1_upgrade() const {return k1_upgrade;};
    bool Getk2_upgrade() const {return k2_upgrade;};
    bool GetMatched() const {return isMatched;};
    bool GetMatched_upgrade() const {return isMatched_upgrade;};

    void SetMatched() {isMatched = true;};
    void SetMatched_upgrade() {isMatched_upgrade = true;};

    double distance(const Tracks& other) const;
    double absoluteError(const Tracks& other) const;
    double z(const Tracks& other) const;

private:
    double x;
    double y;
    double errX;
    double errY;
    int id;
    bool k1;
    bool k2;
    bool k1_upgrade;
    bool k2_upgrade;
    bool isMatched;
    bool isMatched_upgrade;
};

#endif 
