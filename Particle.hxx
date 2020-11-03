#ifndef PARTICLE_H
#define PARTICLE_H
#include "ParticleType.hxx"
#include "ResonanceType.hxx"

class Particle {
public:
Particle(const char* name, double Px, double Py, double Pz);
Particle();
//int Get_maxParticleNumber() const;
static void AddParticleType(const char* name, double mass, int charge, double width = 0);
inline void Set_indexParticle(const char* name);
inline void Set_indexParticle(int index);
inline const int Get_indexParticle() const;
static void PrintArray(); 
void PrintParticle() const;
const double Get_Px() const;
const double Get_Py() const;
const double Get_Pz() const;
const double Get_Mass() const;
const double Energy() const;
inline const double InvMass(const Particle& p) const;
void Set_P(double px,double py,double pz);
int Decay2body(Particle &dau1,Particle &dau2) const;

private:
static const int maxParticleNumber_=10; // number of partycles generated
static ParticleType* ParticleType_[maxParticleNumber_]; //array di puntatori con info dei tipi di particella generati
static int generatedParticle_; //particelle generate: numero elementi non nulli dell'array
int indexParticle_;
double fPx_;
double fPy_;
double fPz_;
inline static int FindParticle(const char* name);
void Boost(double bx, double by, double bz);
};

#endif