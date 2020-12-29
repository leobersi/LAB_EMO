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
  TH1F* typesGenerated = (TH1F*)file->Get("types");
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
   hTheta->Fit("pol0","Q");
   hPhi->Fit("pol0","Q");
   
   subtraction_invMass_PK->Fit("gaus","R","",0.6,1.2);
   subtraction_invMass->Fit("gaus","R","",0.6,1.2);
   invMass_decayment->Fit("gaus","Q");

   subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral());
   subtraction_invMass->SetEntries(subtraction_invMass->Integral());
   
   ////////////////////////////////////////////////////////////////////
   //create FIRST CANVAS with invariant masses results////////////////
   ////////////////////////////////////////////////////////////////////
   TCanvas *c1 = new TCanvas("c1","Generation Results",200,10,1400,900);
   c1->SetFillColor(18);
   auto pad1 = new TPad("pad1","typesGen",0.0125,0.50625,0.49375,0.9875,21);
   auto pad2 = new TPad("pad2", "Impulse",0.0125,0.0125,0.49375,0.49375,21);
   auto pad3 = new TPad("pad3",   "Theta",0.50625,0.50625,0.9875,0.9875,21);
   auto pad4 = new TPad("pad4",     "Phi",0.50625,0.0125,0.9875,0.49375,21);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad4->Draw();

   //working on pad1+++++++++++++++++types particle generated
   pad1->cd();
   pad1->SetGridx(); //setting a grid
   pad1->SetGridy();
   typesGenerated->SetLineColor(4); //histo cosmetics
   typesGenerated->SetLineWidth(2);
   typesGenerated->SetFillColorAlpha(38,50);
   typesGenerated->DrawCopy(); //draw
   //inner pad cosmetics
   pad1->Update();
   pad1->GetFrame()->SetFillColor(42);
   pad1->GetFrame()->SetBorderMode(-1);  
   pad1->GetFrame()->SetBorderSize(5);
   pad1->Modified();

   //working on pad2+++++++++++++++++Impulse distribution
   pad2->cd();
   pad2->SetGridx(); //setting a grid
   pad2->SetGridy();
   hImpulse->SetLineColor(4); //histo cosmetics
   hImpulse->SetFillColorAlpha(38,50);
   hImpulse->DrawCopy(); //draw
   //inner pad cosmetics
   pad2->Update();
   pad2->GetFrame()->SetFillColor(42);
   pad2->GetFrame()->SetBorderMode(-1);  
   pad2->GetFrame()->SetBorderSize(5);
   pad2->Modified();

   //working on pad3+++++++++++++++++Theta distribution
   pad3->cd();
   pad3->SetGridx(); //setting a grid
   pad3->SetGridy();
   hTheta->SetLineColor(4);//histo cosmetics
   hTheta->SetLineWidth(2);
   hTheta->SetFillColorAlpha(38,50);
   hTheta->DrawCopy(); //draw
   //inner pad cosmetics
   pad3->Update();
   pad3->GetFrame()->SetFillColor(42);
   pad3->GetFrame()->SetBorderMode(-1);  
   pad3->GetFrame()->SetBorderSize(5);
   pad3->Modified();

   //working on pad4+++++++++++++++++Phi distribution
   pad4->cd();
   pad4->SetGridx(); //setting a grid
   pad4->SetGridy();
   hPhi->SetLineColor(4);//histo cosmetics
   hPhi->SetLineWidth(2);
   hPhi->SetFillColorAlpha(38,50);
   hPhi->DrawCopy(); //draw
   //inner pad cosmetics
   pad4->Update();
   pad4->GetFrame()->SetFillColor(42);
   pad4->GetFrame()->SetBorderMode(-1);  
   pad4->GetFrame()->SetBorderSize(5);
   pad4->Modified();

   ////////////////////////////////////////////////////////////////////
   //create SECOND CANVAS with invariant masses results////////////////
   ////////////////////////////////////////////////////////////////////
   TCanvas *c2 = new TCanvas("c2","Comparison of Invariant Masses",200,10,1400,900);
   c2->SetFillColor(18);
   auto pad5 = new TPad("pad5","invMassSubtaction3-4",0.0125,0.50625,0.49375,0.9875,21);
   auto pad6 = new TPad("pad6","invMassSubtaction1-2",0.0125,0.0125,0.49375,0.49375,21);
   auto pad7 = new TPad("pad7",       "invMass decay",0.50625,0.50625,0.9875,0.9875,21);
   //auto pad8 = new TPad("pad8","Theta",0.50625,0.0125,0.9875,0.49375,21);
   pad5->Draw();
   pad6->Draw();
   pad7->Draw();
   //pad8->Draw();

   //working on pad5+++++++++++++++++invMassSubtaction3-4
   pad5->cd();
   pad5->SetGridx(); //setting a grid
   pad5->SetGridy();
   subtraction_invMass_PK->GetXaxis()->SetTitle("Invariant Masses");
   subtraction_invMass_PK->GetYaxis()->SetTitle("Number of Particles");
   subtraction_invMass_PK->GetYaxis()->SetTitleOffset(1.3);
   subtraction_invMass_PK->SetLineColor(4); //histo cosmetics
   subtraction_invMass_PK->SetLineWidth(2);
   subtraction_invMass_PK->SetFillColorAlpha(38,50);
   subtraction_invMass_PK->DrawCopy(); //draw
   //inner pad cosmetics
   pad5->Update();
   pad5->GetFrame()->SetFillColor(42);
   pad5->GetFrame()->SetBorderMode(-1);  
   pad5->GetFrame()->SetBorderSize(5);
   pad5->Modified();

   //working on pad6+++++++++++++++++invMassSubtaction1-2
   pad6->cd();
   pad6->SetGridx(); //setting a grid
   pad6->SetGridy();
   subtraction_invMass->GetXaxis()->SetTitle("Invariant Masses");
   subtraction_invMass->GetYaxis()->SetTitle("Number of Particles");
   subtraction_invMass->GetYaxis()->SetTitleOffset(1.3);
   subtraction_invMass->SetLineColor(4); //histo cosmetics
   subtraction_invMass->SetLineWidth(2);
   subtraction_invMass->SetFillColorAlpha(38,50);
   subtraction_invMass->DrawCopy(); //draw
   //inner pad cosmetics
   pad6->Update();
   pad6->GetFrame()->SetFillColor(42);
   pad6->GetFrame()->SetBorderMode(-1);  
   pad6->GetFrame()->SetBorderSize(5);
   pad6->Modified();

   //working on pad7+++++++++++++++++invMass decay
   pad7->cd();
   pad7->SetGridx(); //setting a grid
   pad7->SetGridy();
   invMass_decayment->GetXaxis()->SetTitle("Invariant Masses");
   invMass_decayment->GetYaxis()->SetTitle("Number of Particles");
   invMass_decayment->GetYaxis()->SetTitleOffset(1.3);
   invMass_decayment->SetLineColor(4); //histo cosmetics
   invMass_decayment->SetLineWidth(2);
   invMass_decayment->SetFillColorAlpha(38,50);
   invMass_decayment->DrawCopy(); //draw
   //inner pad cosmetics
   pad7->Update();
   pad7->GetFrame()->SetFillColor(42);
   pad7->GetFrame()->SetBorderMode(-1);  
   pad7->GetFrame()->SetBorderSize(5);
   pad7->Modified();

   file->Close();
}
