#ifndef PARTICLE_TYPE_H
#define PARTICLE_TYPE_H
#include <iostream>
class ParticleType {

public:
ParticleType (const char* name, double mass, int charge);
const char* GetName() const; 
double GetMass() const;
int GetCharge() const;
virtual void Print() const;
virtual double GetWidth() const;

private:
const char* fName_;
const double fMass_;
const int fCharge_;
};
#endif