#include "ParticleType.hxx"
ParticleType::ParticleType (const char* name, double mass, int charge) :
    fName_{name}, 
    fMass_{mass},
    fCharge_{charge} {}
const char* ParticleType::GetName() const {return fName_;}
double ParticleType::GetMass() const {return fMass_;}
int ParticleType::GetCharge() const {return fCharge_;}
void ParticleType::Print() const {
    std::cout<< "++++++++++++++++++++++++++++++++++\n"
             << "Name = " << GetName() << ";  "
             << "Mass = " << GetMass() << ";  "
             << "Charge = " << GetCharge() << '\n'; 
}
double ParticleType::GetWidth() const {return 0;};