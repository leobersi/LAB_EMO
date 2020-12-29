void runGeneration()  {
   gROOT->LoadMacro("Particle.cxx");
   gROOT->LoadMacro("ParticleType.cxx");
   gROOT->LoadMacro("ResonanceType.cxx");
   gROOT->LoadMacro("main.cxx");
  }
void OperationHisto() {
  //open a root file
  TFile* file = new TFile("/home/leonardo/root-build/macros/my_testroot.root"); 
  if ( file->IsOpen() ) cout << "File opened successfully" << '\n';

  //retrieving histograms
  TH1F* hTheta = (TH1F*)file->Get("correlationT");
  TH1F* hPhi = (TH1F*)file->Get("correlationP");
  TH1F* hImpulse = (TH1F*)file->Get("impulse");
  TH1F* invMass_discord = (TH1F*)file->Get("mass1");
  TH1F* invMass_concord = (TH1F*)file->Get("mass2");
  TH1F* invMass_PK_discord = (TH1F*)file->Get("mass3");
  TH1F* invMass_PK_concord = (TH1F*)file->Get("mass4");
  TH1F* invMass_decayment = (TH1F*)file->Get("mass5");

   auto subtraction_invMass_PK = new TH1F(*invMass_PK_discord);
   subtraction_invMass_PK->Sumw2();
   subtraction_invMass_PK->Add(invMass_PK_discord,invMass_PK_concord,1,-1);
   subtraction_invMass_PK->SetTitle("Subtraction of Invariant Mass P&K");
   auto subtraction_invMass = new TH1F(*invMass_discord);
   subtraction_invMass->Sumw2();
   subtraction_invMass->Add(invMass_discord,invMass_concord,1,-1);
   subtraction_invMass->SetTitle("Subtraction of all discord and concord Invariant Mass");
   //fitting
   hTheta->Fit("pol0","Q","0");
   hPhi->Fit("pol0","Q","0");
   hImpulse->Fit("expo","Q","0");
   
   subtraction_invMass_PK->Fit("gaus","R","",0.6,1.2);
   subtraction_invMass->Fit("gaus","R","",0.6,1.2);
   invMass_decayment->Fit("gaus","Q");

   subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral());
   subtraction_invMass->SetEntries(subtraction_invMass->Integral());
   
   //create First canvas1
   TCanvas *canvas1 = new TCanvas("canvas1","Generations",200,10,700,900);
   canvas1->SetFillColor(18);
   canvas1->Divide(3);
   canvas1->SetGridx();
   canvas1->SetGridy();
   
   //canvas1 first part
   canvas1->cd(1);
   canvas1->cd(1)->SetFillColor(42);
   gStyle->SetOptFit(111);
   hTheta->DrawCopy();
   //TCanvas::Update()
   canvas1->cd(1)->Update();
   canvas1->cd(1)->GetFrame()->SetFillColor(21);
   canvas1->cd(1)->Modified();
   //canvas2 second part
   canvas1->cd(2);
   canvas1->cd(2)->SetFillColor(42);
   gStyle->SetOptFit(111);
   hPhi->DrawCopy();
   //TCanvas::Update()
   canvas1->cd(2)->Update();
   canvas1->cd(2)->GetFrame()->SetFillColor(21);
   canvas1->cd(2)->Modified();
   //canvas2 third part
   canvas1->cd(3);
   canvas1->cd(3)->SetFillColor(42);
   gStyle->SetOptFit(111);
   hImpulse->DrawCopy();
   //TCanvas::Update()
   canvas1->cd(3)->Update();
   canvas1->cd(3)->GetFrame()->SetFillColor(21);
   canvas1->cd(3)->Modified();   
///////////////////////////////////////////////////////////////////////////////////////
   //create Second canvas2
   TCanvas *canvas2 = new TCanvas();
   canvas2->Divide(3);
   canvas2->SetGridx();
   canvas2->SetGridy();
   canvas2->SetTitle("Comparison of Invariant Masses");
   canvas2->SetFillColor(48);
   
   //canvas2 first part
   canvas2->cd(1);
   canvas2->cd(1)->SetFillColor(42);
   subtraction_invMass_PK->GetXaxis()->SetTitle("Invariant Masses");
   subtraction_invMass_PK->GetYaxis()->SetTitle("Number of Particles");
   subtraction_invMass_PK->GetYaxis()->SetTitleOffset(1.4); //move title from Yaxes of 1.4 cm  
   //draw
   gStyle->SetOptFit(111);
   subtraction_invMass_PK->DrawCopy();
   //TCanvas::Update()
   canvas2->cd(1)->Update();
   canvas2->cd(1)->GetFrame()->SetFillColor(21);
   canvas2->cd(1)->Modified();

   //canvas2 second part
   canvas2->cd(2);
   canvas2->cd(2)->SetFillColor(42);
   invMass_decayment->GetXaxis()->SetTitle("Invariant Masses");
   invMass_decayment->GetYaxis()->SetTitle("Number of Particles");
   invMass_decayment->GetYaxis()->SetTitleOffset(1.4); //move title from Yaxes of 1.4 cm
   //draw 
   gStyle->SetOptFit(111);
   invMass_decayment->DrawCopy();
   //TCanvas::Update()
   canvas2->cd(2)->Update();
   canvas2->cd(2)->GetFrame()->SetFillColor(21);
   canvas2->cd(2)->Modified();
   
   //canvas2 third part
   canvas2->cd(3);
   canvas2->cd(3)->SetFillColor(42);
   canvas2->cd(3)->GetFrame()->SetFillColor(21);
   subtraction_invMass->GetXaxis()->SetTitle("Invariant Masses");
   subtraction_invMass->GetYaxis()->SetTitle("Number of Particles");
   subtraction_invMass->GetYaxis()->SetTitleOffset(1.4); //move title from Yaxes of 1.4 cm
   //draw
   gStyle->SetOptFit(111);
   subtraction_invMass->DrawCopy();
   //TCanvas::Update()
   canvas2->cd(3)->Update();
   canvas2->cd(3)->GetFrame()->SetFillColor(21);
   canvas2->cd(3)->Modified(); 
  
   file->Close();
}
