#include "ParticleType.hxx"
#include "ResonanceType.hxx"
#include "Particle.hxx"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TCanvas.h"

int main(){
///////////////////////////////////////////////
 gROOT->SetStyle("Plain");
  gStyle-> SetPalette(57);
  gStyle-> SetOptTitle(0);
  TCanvas* mycanvas = new TCanvas();
///////////////////////////////////////////////
gRandom->SetSeed();
constexpr int N = 130;
constexpr int  n_gen = 1e4;
Particle::AddParticleType("P+",0.13957,+1);
Particle::AddParticleType("P-",0.13957,-1);
Particle::AddParticleType("K+",0.49367,+1);
Particle::AddParticleType("K-",0.49367,-1);
Particle::AddParticleType("p+",0.93827,+1);
Particle::AddParticleType("e-",0.93827,-1);
Particle::AddParticleType("K*",0.89166,0,0.050);
Particle particle[N];

TH1F *types = new TH1F("types","Types of Particles generated",7,0,7);
TH2F *correlation = new TH2F("correlation","distribution of angles",1000,0,TMath::Pi(),1000,0,2*TMath::Pi());
TH1F *impulse = new TH1F("impulse","Impulse distribution",1000,0,10);
TH1F *trasv_imp = new TH1F("trasv_imp","Traverse impulse distribution",1000,0,10);
TH1F *energy = new TH1F("energy","Energy of Particles distribution",1000,0,10);
TH1F *mass1 = new TH1F("mass1","Invariant Mass discord charge",1000,0,10);
TH1F *mass2 = new TH1F("mass2","Invariant Mass concord charge",1000,0,10);
TH1F *mass3 = new TH1F("mass3","Invariant Mass P and K discord charge",1000,0,10);
TH1F *mass4 = new TH1F("mass4","Invariant Mass P and K concord charge",1000,0,10);
TH1F *mass5 = new TH1F("mass7","Invariant Mass couples decayed from K*",1000,0,10);
TH1F *mass6 = new TH1F("mass8","Invariant Mass all particles generated",1000,0,10);

double Theta, Phi, Impulse, Px, Py, Pz, x = 0;

for(int i =0; i<n_gen; ++i){
    int decays_counter = 0;			
 for(int j=0;j<100 + decays_counter;++j){
  Theta = gRandom->Uniform(0,TMath::Pi());
  Phi = gRandom->Uniform(0,2*TMath::Pi());
  correlation->Fill(Theta,Phi);

  Impulse = gRandom->Exp(1);
  impulse->Fill(Impulse);

  Pz = Impulse*cos(Theta);
  Px = Impulse*sin(Theta)*cos(Phi);
  Py = Impulse*sin(Theta)*sin(Phi);
  trasv_imp->Fill(sqrt(Px*Px + Py*Py));			
  
  particle[j].Set_P(Px,Py,Pz);
  x = gRandom->Rndm();
   
  if(x < 0.4){
    particle[j].Set_indexParticle(0);
  } else if (x < 0.8){
    particle[j].Set_indexParticle(1);
  } else if (x < 0.85){
    particle[j].Set_indexParticle(2);
  } else if (x < 0.9){
    particle[j].Set_indexParticle(3);
  } else if (x < 0.945){
    particle[j].Set_indexParticle(4);
  } else if (x < 0.99){
    particle[j].Set_indexParticle(5);
  } else if (x < 0.995){
    particle[j].Set_indexParticle(6);
    particle[100 + decays_counter].Set_indexParticle(0);
    ++ decays_counter;
    particle[100 + decays_counter].Set_indexParticle(3);
    ++ decays_counter;
    particle[j].Decay2body(particle[100 + decays_counter-2], particle[100 + decays_counter-1]);
  } else {
    particle[j].Set_indexParticle(6);
    particle[100 + decays_counter].Set_indexParticle(2);
    ++ decays_counter;
    particle[100 + decays_counter].Set_indexParticle(1);
    ++ decays_counter;
    particle[j].Decay2body(particle[100 + decays_counter-2], particle[100 + decays_counter-1]);
  }
  
   if(particle[j].Get_indexParticle() >= 0){
	 types->Fill(particle[j].Get_indexParticle());		//Filling of generated particles		
	}	
       
 }

 for(int j=0; j<N; j++){
    energy->Fill(particle[j].Energy());
    
    if(particle[j].Get_indexParticle() == -2){ ++j;}//<-- in an unitialize
    
    else if (j!= N-1) {
      for (int h = j+1; h<N; h++) {
      mass6->Fill(particle[j].InvMass(particle[h])); //filling mass6
      if((particle[j].Get_indexParticle()%2) == 0 && (particle[j].Get_indexParticle() != 6) 
                                                 && (particle[h].Get_indexParticle() != 6))
      {
       if(particle[h].Get_indexParticle()%2 != 0){
           mass1->Fill(particle[j].InvMass(particle[h])); //filling mass1

       } else if(particle[h].Get_indexParticle()%2 == 0) {
           mass2->Fill(particle[j].InvMass(particle[h])); //filling mass2
      }
      }
      if((particle[j].Get_indexParticle() == 0 && particle[j].Get_indexParticle() == 3) ||
         (particle[j].Get_indexParticle() == 1 && particle[j].Get_indexParticle() == 2)){
           mass3->Fill(particle[j].InvMass(particle[h])); //filling mass3
      }
      if((particle[j].Get_indexParticle() == 1 && particle[j].Get_indexParticle() == 3) ||
         (particle[j].Get_indexParticle() == 0 && particle[j].Get_indexParticle() == 2)){
           mass4->Fill(particle[j].InvMass(particle[h])); //filling mass4
      }
    } 
    }
 }
}
types->Draw("APE");
/*correlation->Draw("APE");
impulse->Draw("APE");
trasv_imp->Draw("APE");
energy->Draw("APE");
mass1->Draw("APE");
mass2->Draw("APE");
mass3->Draw("APE");
mass4->Draw("APE");
mass5->Draw("APE");
mass6->Draw("APE");*/
}
