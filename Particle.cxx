#include "Particle.hxx"
#include <cmath>
#include <cstdlib> //for RAND_MAX

const int Particle::maxParticleNumber_; 
ParticleType* Particle::ParticleType_[maxParticleNumber_];
int Particle::generatedParticle_ = 0;

int Particle::FindParticle(const char* name) {
    for (int i = 0; i < generatedParticle_; ++i) {
      if ((ParticleType_[i]->GetName()) == name) { return i; } 
    }
        //std::cout << "FindParticle failed, index not found\n";
        return -1; 
  }

Particle::Particle(){ indexParticle_ = -2;}

Particle::Particle(const char* name, double Px = 0, double Py = 0, double Pz = 0)
    : fPx_{Px}, fPy_{Py}, fPz_{Pz} {
      int i = FindParticle(name);
      if(i != -1) {indexParticle_ = i;} //errore: legge l'indice sempre a -1
      else {
        indexParticle_ = -1;
        std::cout<<"the "<<name<<" type is not yet contemplated\n";}
    }

void Particle::AddParticleType(const char *name, double mass, int charge, double width) {
  // controllo se numero massimo di particelle già generato
  if (generatedParticle_ == maxParticleNumber_) {
    std::cout << "Maximum number of particles already reached\n";
  }
  // controllo se già presente
  if(generatedParticle_ > 0){
    if (FindParticle(name) > 0) {
      std::cout << "Particle already exists\n";
    }
  }
  if (width != 0) {
    ParticleType_[generatedParticle_] = new ResonanceType(name, mass, charge, width);
    ++ generatedParticle_;
    std::cout<< "ResonanceType "<<name<<" generated with success.\n";
  } else {
    ParticleType_[generatedParticle_] = new ParticleType(name, mass, charge);
    ++ generatedParticle_;
    std::cout<< "ParticleType "<<name<<" generated with success.\n";
  }
}

void Particle::Set_indexParticle(const char *name) {
  if(FindParticle(name) >= 0) { indexParticle_ = FindParticle(name);}
  else {std::cout<< "Set_indexParticle failed\n";}
} 

void Particle::Set_indexParticle(int index) {
  if (index >= 0 && index <= generatedParticle_) {indexParticle_ = index;}
  else {std::cout<< "Set_indexParticle failed\n";}
}

const int Particle::Get_indexParticle() const { return indexParticle_; }

void Particle::PrintArray() {
    for (int i = 0; i < generatedParticle_; i++)
        ParticleType_[i]->Print();
}

void Particle::PrintParticle() const {
  if (indexParticle_ == -1) { std::cout<<"Printing not allowed, type not in storage\n";}
  else {
    std::cout << "++++++++++++++++++++++++++++++++++\n"
            << "Particle's Index = " << indexParticle_
            << "\nName = " << ParticleType_[indexParticle_]->GetName()
            << "\nImpulse components are :\nPx = " << fPx_ 
            << ", Py = " << fPy_ << ", Pz = " << fPz_ << "\n\n";
  }
}

const double Particle::Get_Px() const { return fPx_;}

const double Particle::Get_Py() const { return fPy_;}

const double Particle::Get_Pz() const { return fPz_;}

const double Particle::Get_Mass() const { return ParticleType_[indexParticle_]->GetMass();}

const double Particle::Energy() const {
return sqrt(pow(Get_Mass(),2) + fPx_*fPx_ + fPy_*fPy_ + fPz_*fPz_);
}

const double Particle::InvMass(const Particle& p) const {
  return sqrt(pow(Energy() + p.Energy(), 2) - 
         pow(sqrt(fPx_*fPx_ + fPy_*fPy_ + fPz_*fPz_) + 
         sqrt(p.fPx_*p.fPx_ + p.fPy_*p.fPy_ + p.fPz_*p.fPz_),2));
}
void Particle::Set_P(double Px,double Py,double Pz) {
  fPx_ = Px;
  fPy_ = Py;
  fPz_ = Pz;
}

//**************************************************************************************

int Particle::Decay2body(Particle &dau1,Particle &dau2) const {
  if(Get_Mass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = Get_Mass();
  double massDau1 = dau1.Get_Mass();
  double massDau2 = dau2.Get_Mass();

  if(indexParticle_ > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += ParticleType_[indexParticle_]->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.Set_P(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.Set_P(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx_*fPx_ + fPy_*fPy_ + fPz_*fPz_ + massMot*massMot);

  double bx = fPx_/energy;
  double by = fPy_/energy;
  double bz = fPz_/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}  
//****************************************************************************
void Particle::Boost(double bx, double by, double bz)
{

  double energy = Energy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx_ + by*fPy_ + bz*fPz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx_ += gamma2*bp*bx + gamma*bx*energy;
  fPy_ += gamma2*bp*by + gamma*by*energy;
  fPz_ += gamma2*bp*bz + gamma*bz*energy;
}