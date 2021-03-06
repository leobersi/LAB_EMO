#include "Particle.hxx"
#include "ParticleType.hxx"
#include "ResonanceType.hxx"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TFile.h"

int main() {
  ///////////////////////////////////////////////
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
  ///////////////////////////////////////////////
  gRandom->SetSeed();
  constexpr int nGen = 1e05; //events
  constexpr int nParticle = 100; //number of particles generated per event
  constexpr int N = 130; //dimension of static array

  Particle::AddParticleType("P+", 0.13957, +1); //0
  Particle::AddParticleType("P-", 0.13957, -1); //1
  Particle::AddParticleType("K+", 0.49367, +1); //2
  Particle::AddParticleType("K-", 0.49367, -1); //3
  Particle::AddParticleType("p+", 0.93827, +1); //4
  Particle::AddParticleType("e-", 0.93827, -1); //5
  Particle::AddParticleType("K*", 0.89166, 0, 0.050); //6
  Particle particle[N];

  TFile *file = new TFile("/home/leonardo/root-build/macros/my_testroot.root","RECREATE"); //CAMBIA DIRECTORY ALLA FINE

  //Histograms definition
         TH1F *types  =  new TH1F("types", "Types of Particles generated", 7, 0, 7);
  TH1F *correlationT  =  new TH1F("correlation_of_Theta", "distribution of Azimutal Angle", 100, 0, TMath::Pi());
  TH1F *correlationP  =  new TH1F("correlation_of_Phi", "distribution of Polar Angle", 100, 0, 2 * TMath::Pi());
       TH1F *impulse  =  new TH1F("impulse", "Impulse distribution", 1000, 0, 6);
     TH1F *trasv_imp  =  new TH1F("trasv_imp", "Traverse impulse distribution", 1000, 0, 2);
        TH1F *energy  =  new TH1F("energy", "Energy of Particles distribution", 1000, 0, 2);
         TH1F *mass1  =  new TH1F("mass1", "Invariant Mass discord charge", 80, 0, 2);
         TH1F *mass2  =  new TH1F("mass2", "Invariant Mass concord charge", 80, 0, 2);
         TH1F *mass3  =  new TH1F("mass3", "Invariant Mass P and K discord charge", 80, 0, 2);
         TH1F *mass4  =  new TH1F("mass4", "Invariant Mass P and K concord charge", 80, 0, 2);
         TH1F *mass5  =  new TH1F("mass5", "Invariant Mass couples decayed from K*", 80, 0, 2);
         TH1F *mass6  =  new TH1F("mass6", "Invariant Mass all particles generated", 80, 0, 2);

  double Theta, Phi, Impulse, Px, Py, Pz, x = 0;
  
  for (int i = 0; i < nGen; ++i) { //centomila eventi
  int decays_counter = 0;
    for (int j = 0; j < nParticle; ++j) { //cento particelle
      Theta = gRandom->Uniform(0, TMath::Pi());
      Phi = gRandom->Uniform(0, 2 * TMath::Pi());
      correlationT->Fill(Theta);
      correlationP->Fill(Phi);

      Impulse = gRandom->Exp(1);
      impulse->Fill(Impulse);

      Px = Impulse * sin(Theta) * cos(Phi);
      Py = Impulse * sin(Theta) * sin(Phi);
      Pz = Impulse * cos(Theta);
      trasv_imp->Fill(sqrt(Px * Px + Py * Py));
      particle[j].Set_P(Px, Py, Pz);
      
      x = gRandom->Rndm();

      if (x < 0.4) {
        particle[j].Set_indexParticle(0);
      } else if (x < 0.8) {
        particle[j].Set_indexParticle(1);
      } else if (x < 0.85) {
        particle[j].Set_indexParticle(2);
      } else if (x < 0.9) {
        particle[j].Set_indexParticle(3);
      } else if (x < 0.945) {
        particle[j].Set_indexParticle(4);
      } else if (x < 0.99) {
        particle[j].Set_indexParticle(5);
      } else {
        particle[j].Set_indexParticle(6);
        x = gRandom->Rndm();
          if(x < 0.5) {
          particle[100 + decays_counter].Set_indexParticle(0);
          particle[100 + decays_counter + 1].Set_indexParticle(3);
          particle[j].Decay2body(particle[100 + decays_counter], particle[100 + decays_counter + 1]);
          } else {
            particle[100 + decays_counter].Set_indexParticle(1);
            particle[100 + decays_counter + 1].Set_indexParticle(2);
            particle[j].Decay2body(particle[100 + decays_counter], particle[100 + decays_counter + 1]);
            } 
        decays_counter++;
        decays_counter++;
      }
      if (particle[j].Get_indexParticle() >= 0) {
        types->Fill(particle[j].Get_indexParticle()); // filling of generated particles
      }
    }

    for (int j = 0; j < 100; j++) {
      energy->Fill(particle[j].Energy()); // filling of generated energies
    }

    for (int j = 0; j < 100 + decays_counter; ++j) {
      for (int h = j + 1; h < 100 + decays_counter; h++) {
        mass6->Fill(particle[j].InvMass(particle[h])); // filling mass6
        if ((particle[j].Get_Charge() == 1 && particle[h].Get_Charge() == -1) ||
           (particle[j].Get_Charge() == -1 && particle[h].Get_Charge() == 1)) {
            mass1->Fill(particle[j].InvMass(particle[h])); // filling mass1
          }
        if (particle[j].Get_Charge() == particle[h].Get_Charge()) {
            mass2->Fill(particle[j].InvMass(particle[h])); // filling mass2
          }
        if (((particle[j].Get_indexParticle() == 0 && particle[h].Get_indexParticle() == 3) || //discordi
            (particle[j].Get_indexParticle() == 1 && particle[h].Get_indexParticle() == 2))) { 
          mass3->Fill(particle[j].InvMass(particle[h])); // filling mass3
        }
        if (((particle[j].Get_indexParticle() == 1 && particle[h].Get_indexParticle() == 3) || //concordi
            (particle[j].Get_indexParticle() == 0 && particle[h].Get_indexParticle() == 2))) {
          mass4->Fill(particle[j].InvMass(particle[h])); // filling mass4
        }
        if((j >= nParticle && h >= nParticle) && (j%2 == 0 && h == j+1)) { // couple of contigous decayed particles 
					mass5-> Fill(particle[j].InvMass(particle[h]));
				}  
      }
    }
  }
  
  file->Write();
  file->Close();
}
