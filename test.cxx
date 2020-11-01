#include "ParticleType.cxx"
#include "ResonanceType.cxx"
#include "Particle.cxx"

int main() {
/*ParticleType particle("Pione",27.03,0);
ResonanceType resonance("Electron",55.0,2,8);
std::cout << "27.03 =" << particle.GetMass() << '\n'
          << "55.0  =" << resonance.GetMass() << '\n';
std::cout <<"base Print()" << '\n';
particle.Print();
std::cout <<"derived Print()" << '\n';
resonance.Print();
          
ParticleType* a[2];
a[0] = &particle;
a[1] = &resonance; 
for (int i=0; i<2; i++) { a[i]->Print(); }*/


Particle::AddParticleType("K+", 0.05, 1);
Particle::AddParticleType("K-", 0.05, -1);
Particle::AddParticleType("P+", 0.07, 1);
Particle::AddParticleType("P-", 0.07, -1);
Particle::AddParticleType("K+", 0.05, 1, 0.5);
//Particle::PrintArray();

Particle e("K-",3,6,1);
//Particle particle[2];

//e.Set_indexParticle("P-");
//std::cout<<particle[0].Get_indexParticle() ;

e.PrintParticle();

//std::cout<<e.Energy();

//particle.Set_indexParticle("elettrone");
//std::cout<< particle.Get_indexParticle()<< '\n';

}

