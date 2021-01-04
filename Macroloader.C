void runGeneration()  {
   gROOT->LoadMacro("Particle.cxx");
   gROOT->LoadMacro("ParticleType.cxx");
   gROOT->LoadMacro("ResonanceType.cxx");
   gROOT->LoadMacro("main.cxx");
  }
void OperationHisto() {
   gStyle->SetOptStat(2210);

  //open a root file
  TFile* file = new TFile("/home/leonardo/root-build/macros/my_testroot.root"); 
  if ( file->IsOpen() ) cout << "File opened successfully" << '\n';

  //retrieving histograms
  TH1F* typesGenerated = (TH1F*)file->Get("types");
  TH1F* hTheta = (TH1F*)file->Get("correlation_of_Theta");
  TH1F* hPhi = (TH1F*)file->Get("correlation_of_Phi");
  TH1F* hImpulse = (TH1F*)file->Get("impulse");
  TH1F* invMass_discord = (TH1F*)file->Get("mass1");
  TH1F* invMass_concord = (TH1F*)file->Get("mass2");
  TH1F* invMass_PK_discord = (TH1F*)file->Get("mass3");
  TH1F* invMass_PK_concord = (TH1F*)file->Get("mass4");
  TH1F* invMass_decayment = (TH1F*)file->Get("mass5");
   
   //analysis: subtractions
   auto subtraction_invMass_PK = new TH1F("subtraction invMass PK","Subtraction of Invariant Masses P&K",80, 0, 2);
   subtraction_invMass_PK->Sumw2(); 
   subtraction_invMass_PK->Add(invMass_PK_discord,invMass_PK_concord,1,-1);
   auto subtraction_invMass = new TH1F("subtraction invMass","Subtraction of all discord and concord Invariant Masses",80, 0, 2);
   subtraction_invMass->Sumw2();
   subtraction_invMass->Add(invMass_discord,invMass_concord,1,-1);
   
   //fitting
   hTheta->Fit("pol0","Q");
   hPhi->Fit("pol0","Q");
   hImpulse->Fit("expo","Q");
   subtraction_invMass_PK->Fit("gaus","Q","SAME",0.6,1.2);
   subtraction_invMass->Fit("gaus","Q","SAME",0.6,1.2);
   invMass_decayment->Fit("gaus","Q","SAME",0.6,1.2);

   subtraction_invMass_PK->SetEntries(subtraction_invMass_PK->Integral()); //adjusting entries
   subtraction_invMass->SetEntries(subtraction_invMass->Integral());

   
   ////////////////////////////////////////////////////////////////////
   //create FIRST CANVAS with invariant masses results/////////////////
   ////////////////////////////////////////////////////////////////////
   TCanvas *c1 = new TCanvas("c1","Generation Results",200,10,1400,900);
   auto pad1 = new TPad("pad1","typesGen",0.0125,0.50625,0.49375,0.9875); //creating pads were histos will be drawn
   auto pad2 = new TPad("pad2", "Impulse",0.0125,0.0125,0.49375,0.49375);
   auto pad3 = new TPad("pad3",   "Theta",0.50625,0.50625,0.9875,0.9875);
   auto pad4 = new TPad("pad4",     "Phi",0.50625,0.0125,0.9875,0.49375);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad4->Draw();

   //working on pad1+++++++++++++++++types particle generated
   pad1->cd();
   pad1->SetGridx(); //setting a grid
   pad1->SetGridy();
   typesGenerated->GetXaxis()->SetTitle("Angle Values");
   typesGenerated->GetYaxis()->SetTitle("Number of Particles");
   typesGenerated->GetYaxis()->SetTitleOffset(1.3);
   typesGenerated->SetLineColor(65); //histo cosmetics
   typesGenerated->SetLineWidth(2);
   typesGenerated->SetFillColorAlpha(68,0.7);
   typesGenerated->DrawCopy(); //draw
   //inner pad cosmetics
   pad1->Update();
   pad1->GetFrame()->SetFillColor(19);
   pad1->GetFrame()->SetBorderMode(-1);  
   pad1->GetFrame()->SetBorderSize(5);
   pad1->Modified();

   //working on pad2+++++++++++++++++Impulse distribution
   pad2->cd();
   pad2->SetGridx(); //setting a grid
   pad2->SetGridy();
   hImpulse->GetXaxis()->SetTitle("Impulse (GeV)");
   hImpulse->GetYaxis()->SetTitle("Number of Particles");
   hImpulse->GetYaxis()->SetTitleOffset(1.3);
   hImpulse->SetLineColor(65); //histo cosmetics
   hImpulse->SetLineWidth(2);
   hImpulse->SetFillColorAlpha(68,0.7);
   gStyle->SetOptFit(111);
   hImpulse->DrawCopy(); //draw
   //inner pad cosmetics
   pad2->Update();
   pad2->GetFrame()->SetFillColor(19);
   pad2->GetFrame()->SetBorderMode(-1);  
   pad2->GetFrame()->SetBorderSize(5);
   pad2->Modified();
   

   //working on pad3+++++++++++++++++Theta distribution
   pad3->cd();
   pad3->SetGridx(); //setting a grid
   pad3->SetGridy();
   hTheta->GetXaxis()->SetTitle("Angle Values (rad)");
   hTheta->GetYaxis()->SetTitle("Number of Particles");
   hTheta->GetYaxis()->SetTitleOffset(1.3);
   hTheta->SetLineColor(65);//histo cosmetics
   hTheta->SetLineWidth(2);
   gStyle->SetOptFit(111);
   hTheta->DrawCopy(); //draw
   //inner pad cosmetics
   pad3->Update();
   pad3->GetFrame()->SetFillColor(19);
   pad3->GetFrame()->SetBorderMode(-1);  
   pad3->GetFrame()->SetBorderSize(5);
   pad3->Modified();
   

   //working on pad4+++++++++++++++++Phi distribution
   pad4->cd();
   pad4->SetGridx(); //setting a grid
   pad4->SetGridy();
   hPhi->GetXaxis()->SetTitle("Angle Values (rad)");
   hPhi->GetYaxis()->SetTitle("Number of Particles");
   hPhi->GetYaxis()->SetTitleOffset(1.3);
   hPhi->SetLineColor(65);//histo cosmetics
   hPhi->SetLineWidth(2);
   gStyle->SetOptFit(111);
   hPhi->DrawCopy(); //draw
   //inner pad cosmetics
   pad4->Update();
   pad4->GetFrame()->SetFillColor(19);
   pad4->GetFrame()->SetBorderMode(-1);  
   pad4->GetFrame()->SetBorderSize(5);
   pad4->Modified();
   

   ////////////////////////////////////////////////////////////////////
   //create SECOND CANVAS with invariant masses results////////////////
   ////////////////////////////////////////////////////////////////////
   TCanvas *c2 = new TCanvas("c2","Comparison of Invariant Masses",200,10,1400,900);
   auto pad5 = new TPad("pad5","invMassSubtaction3-4",0.0125,0.5,0.5,0.9875); //creating pads were histos will be drawn
   auto pad6 = new TPad("pad6","invMassSubtaction1-2",0.0125,0.0125,0.5,0.5);
   auto pad7 = new TPad("pad7",       "invMass decay",0.5,0.5,0.9875,0.9875);
   pad5->Draw();
   pad6->Draw();
   pad7->Draw();

   //working on pad5+++++++++++++++++invMassSubtaction3-4
   pad5->cd();
   pad5->SetGridx(); //setting a grid
   pad5->SetGridy();
   subtraction_invMass_PK->GetXaxis()->SetTitle("Invariant Masses");
   subtraction_invMass_PK->GetYaxis()->SetTitle("Number of Particles");
   subtraction_invMass_PK->GetYaxis()->SetTitleOffset(1.3);
   subtraction_invMass_PK->SetLineColor(65); //histo cosmetics
   subtraction_invMass_PK->SetLineWidth(2);
   subtraction_invMass_PK->SetFillColorAlpha(68,0.7);
   gStyle->SetOptFit(1111);
   subtraction_invMass_PK->DrawCopy("hist"); //draw
   subtraction_invMass_PK->GetFunction("gaus")->DrawCopy("SAME"); //clean draw without error bars
   //inner pad cosmetics
   pad5->Update();
   pad5->GetFrame()->SetFillColor(19);
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
   subtraction_invMass->SetLineColor(65); //histo cosmetics
   subtraction_invMass->SetLineWidth(2);
   subtraction_invMass->SetFillColorAlpha(68,0.7);
   gStyle->SetOptFit(1111);
   subtraction_invMass->DrawCopy("hist"); //draw
   subtraction_invMass->GetFunction("gaus")->DrawCopy("SAME"); //clean draw without error bars
   //inner pad cosmetics
   pad6->Update();
   pad6->GetFrame()->SetFillColor(19);
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
   invMass_decayment->SetLineColor(65); //histo cosmetics
   invMass_decayment->SetLineWidth(2);
   invMass_decayment->SetFillColorAlpha(68,0.7);
   gStyle->SetOptFit(1111);
   invMass_decayment->DrawCopy("hist"); //draw
   invMass_decayment->GetFunction("gaus")->DrawCopy("SAME"); //clean draw without error bars
   //inner pad cosmetics
   pad7->Update();
   pad7->GetFrame()->SetFillColor(19);
   pad7->GetFrame()->SetBorderMode(-1);  
   pad7->GetFrame()->SetBorderSize(5);
   pad7->Modified();
   
   //Printing
   for(int i=1; i<=typesGenerated->GetNbinsX();i++){
    std::cout<<"Numbers of types "<< i << " generated is: " << typesGenerated->GetBinContent(i) <<'\n';
    std::cout<<"with an error: "<< typesGenerated->GetBinError(i) <<'\n';
  }
   //c1->Print("Generation Results.gif");
   //c2->Print("Comparison of Invariant Masses.gif");
   file->Close();
}
