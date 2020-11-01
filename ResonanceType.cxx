#include "ResonanceType.hxx"
#include "ParticleType.hxx"
ResonanceType::ResonanceType(const char* name, double mass, int charge, double width) : ParticleType (name, mass, charge), fWidth_{width} {}
double ResonanceType::GetWidth() const {return fWidth_;}
void ResonanceType::Print() const {
    ParticleType::Print(); 
    std::cout<<"Resonance Width = " << GetWidth() <<'\n';
};

