#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TH2.h"

	Char_t* filename = Form("Tree-Analysis.root");
	
	TFile *file1 = new TFile(filename);
	TTree *tree1 = (TTree*)file1->Get("Pi0ComptonTree");

void Main() {

	ShowAllDistributions();
	MissingMass();

}

void ShowAllDistributions() {

       	Int_t NPhotTemp;
       	Double_t RawPhotonTheta[100], RawPhotonPhi[100], RawPhotonEnergy[100];

       	Int_t NPhoton;
       	Double_t PhotonTheta[100], PhotonPhi[100], PhotonEnergy[100];
	
       	Int_t NPi0;
       	Double_t Pi0Theta[100], Pi0Phi[100], Pi0Energy[100];

		Int_t NPrompt;
		Int_t TaggerChannelPrompt[100];

		Int_t NRandom;
		Int_t TaggerChannelRandom[100];
	
        tree1->SetBranchAddress("NPhotTemp", 	&NPhotTemp);
        tree1->SetBranchAddress("RawPhotonPhi",    &RawPhotonPhi);
        tree1->SetBranchAddress("RawPhotonTheta",  &RawPhotonTheta);
        tree1->SetBranchAddress("RawPhotonEnergy", &RawPhotonEnergy);

        tree1->SetBranchAddress("NPhoton", 	&NPhoton);
        tree1->SetBranchAddress("PhotonPhi",    &PhotonPhi);
        tree1->SetBranchAddress("PhotonTheta",  &PhotonTheta);
        tree1->SetBranchAddress("PhotonEnergy", &PhotonEnergy);

        tree1->SetBranchAddress("NPi0", 	&NPi0);
        tree1->SetBranchAddress("Pi0Phi",    &Pi0Phi);
        tree1->SetBranchAddress("Pi0Theta",  &Pi0Theta);
        tree1->SetBranchAddress("Pi0Energy", &Pi0Energy);

		tree1->SetBranchAddress("NPrompt", &NPrompt);
		tree1->SetBranchAddress("TaggerChannelPrompt", &TaggerChannelPrompt);

		tree1->SetBranchAddress("NRandom", &NRandom);	
		tree1->SetBranchAddress("TaggerChannelRandom", &TaggerChannelRandom);

	

        TH1D *a0 = new TH1D("NPhotTemp",          "NPhotTemp",	       5,         0,         5);
        TH1D *a1 = new TH1D("RawPhotonPhi",       "RawPhotonPhi",    360,    -180,        180);
        TH1D *a2 = new TH1D("RawPhotonTheta",     "RawPhotonTheta",  180,       0,        180);
        TH1D *a3 = new TH1D("RawPhotonEnergy",    "RawPhotonEnergy", 400,       0,        400);

        TH1D *b0 = new TH1D("NPhoton",          "NPhoton",  	5,         0,         5);
        TH1D *b1 = new TH1D("PhotonPhi",        "PhotonPhi",    360,    -180,        180);
        TH1D *b2 = new TH1D("PhotonTheta",      "PhotonTheta",  180,       0,        180);
        TH1D *b3 = new TH1D("PhotonEnergy",     "PhotonEnergy", 400,       0,        400);

        TH1D *c0 = new TH1D("NPi0",          "NPi0",  	5,         0,         5);
        TH1D *c1 = new TH1D("Pi0Phi",        "Pi0Phi",    360,    -180,        180);
        TH1D *c2 = new TH1D("Pi0Theta",      "Pi0Theta",  180,       0,        180);
        TH1D *c3 = new TH1D("Pi0Energy",     "Pi0Energy", 400,       0,        400);

        TH1D *d0 = new TH1D("NPrompt",       "NPrompt",  	5,         0,         5);
		TH1D *d1 = new TH1D("Tagger-channel-prompt","Tagger-channel-prompt", 352,0,352);        
		
        TH1D *e0 = new TH1D("NRandom",       "NRandom",  	5,         0,         5);
		TH1D *e1 = new TH1D("Tagger-channel-random","Tagger-channel-random", 352,0,352);        


        Int_t nentries = (Int_t)tree1->GetEntries();
        for (Int_t i=0;i<nentries;i++) {

                tree1->GetEntry(i);
                a0->Fill(NPhotTemp);
                
                for (Int_t a = 0; a < NPhotTemp; a++){

                        a1->Fill(RawPhotonPhi[a]);
                        a2->Fill(RawPhotonTheta[a]);
                        a3->Fill(RawPhotonEnergy[a]);
                }       

                b0->Fill(NPhoton);
                
                for (Int_t a = 0; a < NPhoton; a++){

                        b1->Fill(PhotonPhi[a]);
                        b2->Fill(PhotonTheta[a]);
                        b3->Fill(PhotonEnergy[a]);
                }    

                c0->Fill(NPi0);
                
                for (Int_t a = 0; a < NPi0; a++){

                        c1->Fill(Pi0Phi[a]);
                        c2->Fill(Pi0Theta[a]);
                        c3->Fill(Pi0Energy[a]);
                }
                
                d0->Fill(NPrompt);

                for (Int_t a = 0; a < NPrompt; a++){

						d1->Fill(TaggerChannelPrompt[a]);
                }
                
                e0->Fill(NRandom);

                for (Int_t a = 0; a < NRandom; a++){

						e1->Fill(TaggerChannelRandom[a]);
                }                
                     
        }

	canvasRawPhot = new TCanvas("Raw Photon Dist","Distributions of all photons (pre pi0 analysis)");
	canvasRawPhot->Divide(2,2);
        canvasRawPhot->cd(1); a0->Draw();
        canvasRawPhot->cd(2); a1->Draw();
        canvasRawPhot->cd(3); a2->Draw();
        canvasRawPhot->cd(4); a3->Draw();

	canvasPhot = new TCanvas("Photon Dist","Distributions of Photons (post pi0 analysis)");
	canvasPhot->Divide(2,2);
        canvasPhot->cd(1); b0->Draw();
        canvasPhot->cd(2); b1->Draw();
        canvasPhot->cd(3); b2->Draw();
        canvasPhot->cd(4); b3->Draw();

	canvasPi0 = new TCanvas("Pi0 Dist","Distributions of Pi0s (post pi0 analysis)");
	canvasPi0->Divide(2,2);
        canvasPi0->cd(1); c0->Draw();
        canvasPi0->cd(2); c1->Draw();
        canvasPi0->cd(3); c2->Draw();
        canvasPi0->cd(4); c3->Draw();
        
	canvasNPrompt = new TCanvas("NPrompt","NPrompt something-something");
	canvasNPrompt->Divide(2);
        canvasNPrompt->cd(1); d0->Draw();
        canvasNPrompt->cd(2); d1->Draw();

	canvasNRandom = new TCanvas("NRandom","NRandom something-something");
	canvasNRandom->Divide(2);
        canvasNRandom->cd(1); e0->Draw();
        canvasNRandom->cd(2); e1->Draw();

}

void RawPhotonDistributions() {

        Int_t NPhotTemp;
        Double_t RawPhotonTheta[100], RawPhotonPhi[100], RawPhotonEnergy[100];

        tree1->SetBranchAddress("NPhotTemp", 	&NPhotTemp);
        tree1->SetBranchAddress("RawPhotonPhi",    &RawPhotonPhi);
        tree1->SetBranchAddress("RawPhotonTheta",  &RawPhotonTheta);
        tree1->SetBranchAddress("RawPhotonEnergy", &RawPhotonEnergy);

        TH1D *a0 = new TH1D("NPhotTemp",          "NPhotTemp",	       5,         0,         5);
        TH1D *a1 = new TH1D("RawPhotonPhi",       "RawPhotonPhi",    360,    -180,        180);
        TH1D *a2 = new TH1D("RawPhotonTheta",     "RawPhotonTheta",  180,       0,        180);
        TH1D *a3 = new TH1D("RawPhotonEnergy",    "RawPhotonEnergy", 400,       0,        400);

        Int_t nentries = (Int_t)tree1->GetEntries();
        for (Int_t i=0;i<nentries;i++) {

                tree1->GetEntry(i);
                a0->Fill(NPhotTemp);
                
                for (Int_t a = 0; a < NPhotTemp; a++){

                        a1->Fill(RawPhotonPhi[a]);
                        a2->Fill(RawPhotonTheta[a]);
                        a3->Fill(RawPhotonEnergy[a]);
                }      
        }

	canvasRawPhot = new TCanvas("Raw Photon Dist","Distributions of all photons (pre pi0 analysis)");
	canvasRawPhot->Divide(2,2);
        canvasRawPhot->cd(1); a0->Draw();
        canvasRawPhot->cd(2); a1->Draw();
        canvasRawPhot->cd(3); a2->Draw();
        canvasRawPhot->cd(4); a3->Draw();
}


void PhotonDistributions() {

        Int_t NPhoton;
        Double_t PhotonTheta[100], PhotonPhi[100], PhotonEnergy[100];

        tree1->SetBranchAddress("NPhoton", 	&NPhoton);
        tree1->SetBranchAddress("PhotonPhi",    &PhotonPhi);
        tree1->SetBranchAddress("PhotonTheta",  &PhotonTheta);
        tree1->SetBranchAddress("PhotonEnergy", &PhotonEnergy);

        TH1D *a0 = new TH1D("NPhoton",          "NPhoton",  	5,         0,         5);
        TH1D *a1 = new TH1D("PhotonPhi",        "PhotonPhi",    360,    -180,        180);
        TH1D *a2 = new TH1D("PhotonTheta",      "PhotonTheta",  180,       0,        180);
        TH1D *a3 = new TH1D("PhotonEnergy",     "PhotonEnergy", 400,       0,        400);

        Int_t nentries = (Int_t)tree1->GetEntries();
        for (Int_t i=0;i<nentries;i++) {

                tree1->GetEntry(i);
                a0->Fill(NPhoton);
                
                for (Int_t a = 0; a < NPhoton; a++){

                        a1->Fill(PhotonPhi[a]);
                        a2->Fill(PhotonTheta[a]);
                        a3->Fill(PhotonEnergy[a]);
                }        
        }

	canvasPhot = new TCanvas("Photon Dist","Distributions of Photons (post pi0 analysis)");
	canvasPhot->Divide(2,2);
        canvasPhot->cd(1); a0->Draw();
        canvasPhot->cd(2); a1->Draw();
        canvasPhot->cd(3); a2->Draw();
        canvasPhot->cd(4); a3->Draw();
}

void Pi0Distributions() {

        Int_t NPi0;
        Double_t Pi0Theta[100], Pi0Phi[100], Pi0Energy[100];

        tree1->SetBranchAddress("NPi0", 	&NPi0);
        tree1->SetBranchAddress("Pi0Phi",    &Pi0Phi);
        tree1->SetBranchAddress("Pi0Theta",  &Pi0Theta);
        tree1->SetBranchAddress("Pi0Energy", &Pi0Energy);
	

        TH1D *a0 = new TH1D("NPi0",          "NPi0",  	5,         0,         5);
        TH1D *a1 = new TH1D("Pi0Phi",        "Pi0Phi",    360,    -180,        180);
        TH1D *a2 = new TH1D("Pi0Theta",      "Pi0Theta",  180,       0,        180);
        TH1D *a3 = new TH1D("Pi0Energy",     "Pi0Energy", 400,       0,        400);

        Int_t nentries = (Int_t)tree1->GetEntries();
        for (Int_t i=0;i<nentries;i++) {

                tree1->GetEntry(i);
                a0->Fill(NPi0);
                
                for (Int_t a = 0; a < NPi0; a++){

                        a1->Fill(Pi0Phi[a]);
                        a2->Fill(Pi0Theta[a]);
                        a3->Fill(Pi0Energy[a]);
                }        
        }

	canvasPi0 = new TCanvas("Pi0 Dist","Distributions of Pi0s (post pi0 analysis)");
	canvasPi0->Divide(2,2);
        canvasPi0->cd(1); a0->Draw();
        canvasPi0->cd(2); a1->Draw();
        canvasPi0->cd(3); a2->Draw();
        canvasPi0->cd(4); a3->Draw();
}

void MissingMass(Int_t TC_low = 0, Int_t TC_high = 351) {


	gStyle->SetOptStat(1111);

	Int_t NPromptPi0, NRandomPi0;
	Int_t TaggerChannelPromptPi0[1000], TaggerChannelRandomPi0[1000];
	Double_t MissingMassPromptPi0[1000], MissingMassRandomPi0[1000];

	tree1->SetBranchAddress("NPromptPi0",		&NPromptPi0);
	tree1->SetBranchAddress("NRandomPi0",		&NRandomPi0);

	tree1->SetBranchAddress("TaggerChannelPromptPi0",&TaggerChannelPromptPi0);
	tree1->SetBranchAddress("TaggerChannelRandomPi0",&TaggerChannelRandomPi0);

	tree1->SetBranchAddress("MissingMassPromptPi0",	&MissingMassPromptPi0);
	tree1->SetBranchAddress("MissingMassRandomPi0",	&MissingMassRandomPi0);

	TH1D *a0 = new TH1D("Missing-Mass-Pi0-Prompt","Missing-Mass-Pi0-Prompt", 400,800,1200);
	TH1D *a1 = new TH1D("Missing-Mass-Pi0-Random","Missing-Mass-Pi0-Random", 400,800,1200);
	TH1D *a2 = new TH1D("Missing-Mass-Pi0",	      "Missing-Mass-Pi0",	 400,800,1200);

	Int_t nentries = (Int_t)tree1->GetEntries();
	for (Int_t i=0;i<nentries;i++) {

		tree1->GetEntry(i);
		
		for (Int_t a = 0; a < NPromptPi0; a++){
			if (TaggerChannelPromptPi0[a] >= TC_low && TaggerChannelPromptPi0[a] <= TC_high){
				a0->Fill(MissingMassPromptPi0[a]);
			}
		}	

		for (Int_t a = 0; a < NRandomPi0; a++){
			if (TaggerChannelRandomPi0[a] >= TC_low && TaggerChannelRandomPi0[a] <= TC_high){
				a1->Fill(MissingMassRandomPi0[a]);
			}	
		}
	}
	
	canvasMM = new TCanvas("canvasMM","Missing-Mass-Pi0");
	a2->Add(a0,1);
	a2->Add(a1,-0.125); 
	a2->Draw();

}

void PromptPhotonAnalysis()
{
	Int_t NPrompt;
	Int_t TaggerChannelPrompt[100];
	
	Int_t NRandom;
	Int_t TaggerChannelRandom[100];
	
	Int_t Peak;
	
	tree1->SetBranchAddress("NPrompt", &NPrompt);
	tree1->SetBranchAddress("TaggerChannelPrompt", &TaggerChannelPrompt);

	tree1->SetBranchAddress("NRandom", &NRandom);	
	tree1->SetBranchAddress("TaggerChannelRandom", &TaggerChannelRandom);

    TH1D *d0 = new TH1D("NPrompt",       "NPrompt",  	5,         0,         5);
	TH1D *d1 = new TH1D("Tagger-channel-prompt","Tagger-channel-prompt", 352,0,352);        
		
    TH1D *e0 = new TH1D("NRandom",       "NRandom",  	5,         0,         5);
	TH1D *e1 = new TH1D("Tagger-channel-random","Tagger-channel-random", 352,0,352);        

	//variables to count pions
	Int_t NPi0;
    Double_t Pi0Theta[100], Pi0Phi[100], Pi0Energy[100];
       
    tree1->SetBranchAddress("NPi0", 	&NPi0);
    tree1->SetBranchAddress("Pi0Phi",    &Pi0Phi);
    tree1->SetBranchAddress("Pi0Theta",  &Pi0Theta);
    tree1->SetBranchAddress("Pi0Energy", &Pi0Energy);


	//Histogram that counts pions
	TH1D *PionCount = new TH1D("PionCount","PionCount", 6,0,6);

    Int_t nentries = (Int_t)tree1->GetEntries();

    for (Int_t i = 0; i < nentries; i++)
    {
		tree1->GetEntry(i);
                
                
        d0->Fill(NPrompt);

        for (Int_t a = 0; a < NPrompt; a++)
        {
			if(NPi0 > 0)
			{
				//Tagger channels corresponding to photon energies 325 - 225 MeV
				if ((TaggerChannelPrompt[0] > 295) && (TaggerChannelPrompt[0] < 318))
				{
					Peak = 0;
				}
				
				//Tagger channels corresponding to photon energies 325 - 425 MeV
				if ((TaggerChannelPrompt[0] > 272) && (TaggerChannelPrompt[0] < 295))
				{
					Peak = 1;
				}

				//Tagger channels corresponding to photon energies 425 - 525 MeV
				if ((TaggerChannelPrompt[0] > 250) && (TaggerChannelPrompt[0] < 272))
				{
					Peak = 2;
				}

				//Tagger channels corresponding to photon energies 525 - 625 MeV
				if ((TaggerChannelPrompt[0] > 227) && (TaggerChannelPrompt[0] < 250))
				{
					Peak = 3;
				}

				//Tagger channels corresponding to photon energies 625 - 725 MeV
				if ((TaggerChannelPrompt[0] > 204) && (TaggerChannelPrompt[0] < 227))
				{
					Peak = 4;
				}

				//Tagger channels corresponding to photon energies 725 - 825 MeV
				if ((TaggerChannelPrompt[0] > 180) && (TaggerChannelPrompt[0] < 204))
				{
					Peak = 5;
				}
				
				PionCount->Fill(Peak);
			}
				
			d1->Fill(TaggerChannelPrompt[a]);
			
        }

               
        e0->Fill(NRandom);

        for (Int_t a = 0; a < NRandom; a++)
        {
			e1->Fill(TaggerChannelRandom[a]);
        }                
     }
     
   	 canvasNPrompt = new TCanvas("NPrompt","NPrompt something-something");
	 canvasNPrompt->Divide(2);
     canvasNPrompt->cd(1); d0->Draw();
     canvasNPrompt->cd(2); d1->Draw();

	 canvasNRandom = new TCanvas("NRandom","NRandom something-something");
	 canvasNRandom->Divide(2);
     canvasNRandom->cd(1); e0->Draw();
	 canvasNRandom->cd(2); e1->Draw();

	 canvasPionCount = new TCanvas("PionCount","PionCount");
     PionCount->Draw();


}

