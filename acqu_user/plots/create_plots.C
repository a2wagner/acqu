#include <iostream>
#include <stdlib.h>
//#include <cmath>
//#include <string>
#include <vector>
#include <map>

#include <TROOT.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
//#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TPaveStats.h>
//#include <TLatex.h>

#define N_WINDOWS 2
#define PROMPT 0
#define RANDOM 1

typedef std::map<int, const char*> IntCharMap;
typedef std::pair<int, const char*> ICPair;
typedef std::map<int, const char*>::iterator ICIter;


void prepare_hist(TH1 *h)
{
	//h->SetContour(50);
	//h->GetXaxis()->SetTitle("inv. Masse [MeV]");
	//h->GetXaxis()->SetRange(45,250);
	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(0.045);
	h->GetXaxis()->SetTitleSize(0.048);
	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitle("#Events");
	//h->GetYaxis()->SetRange(10,320);
	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelSize(0.045);
	h->GetYaxis()->SetTitleSize(0.048);
	h->GetYaxis()->SetTitleFont(42);
	h->GetZaxis()->SetLabelFont(42);
	h->GetZaxis()->SetLabelSize(0.042);
	h->GetZaxis()->SetTitleSize(0.035);
	h->GetZaxis()->SetTitleFont(42);
	h->SetTitle("");
	h->GetYaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetLabelOffset(.008);  // per cent of pad width; standard is .005
	h->GetYaxis()->SetDecimals();  //show e. g. 1.0 instead of just 1 (same decimals for every label)
	h->GetYaxis()->SetNoExponent(false);  //show exponent for very large or small values (seems only to work for values smaller than e-5)
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetLabelOffset(.01);  // per cent of pad width; standard is .005
}

int create_plots(const char* file = "/home/wagners/acqu/acqu_user/PHYS_CB_2014-06-06_completeSet.root")
{
	/* some initial stuff */
	const char* suffix = "pdf";  // file format for saving the histograms

	gStyle->SetCanvasColor(0);
	//gStyle->SetOptStat(0);

	// Color used for histograms
	Int_t color = TColor::GetColor("#014ED0");

	// enum for numbering of different cuts and cut combinations
	enum cut {
		protE,
		copl,
		balance,
		dAlpha,
		missM,
		invM,
		copl_balance,
		balance_missM,
		protE_copl,
		copl_missM,
		missM_invM,
		balance_dAlpha,
		copl_dAlpha,
		protE_copl_balance,
		copl_balance_dAlpha,
		all,
		unknown
	};

	IntCharMap cuts;
	cuts.insert(ICPair(protE, "protE"));
	cuts.insert(ICPair(copl, "copl"));
	cuts.insert(ICPair(balance, "balance"));
	cuts.insert(ICPair(dAlpha, "dAlpha"));
	cuts.insert(ICPair(missM, "missM"));
	cuts.insert(ICPair(invM, "invM"));
	cuts.insert(ICPair(copl_balance, "copl_balance"));
	cuts.insert(ICPair(balance_missM, "balance_missM"));
	cuts.insert(ICPair(protE_copl, "protE_copl"));
	cuts.insert(ICPair(copl_missM, "copl_missM"));
	cuts.insert(ICPair(missM_invM, "missM_invM"));
	cuts.insert(ICPair(balance_dAlpha, "balance_dAlpha"));
	cuts.insert(ICPair(copl_dAlpha, "copl_dAlpha"));
	cuts.insert(ICPair(protE_copl_balance, "protE_copl_balance"));
	cuts.insert(ICPair(copl_balance_dAlpha, "copl_balance_dAlpha"));
	cuts.insert(ICPair(all, "allCuts"));

	std::vector<const char*> props;
	props.push_back("invM");
	props.push_back("missM");

	const char* windows[2] = {"prompt", "random"};

	// connect to the file to read from
	TFile f(file, "READ");
	if (!f.IsOpen()) {
		printf("Error opening file %s!\n", file);
		return 1;
	}
	printf("Read ");
	f.Print();
	// load the desired folder from file (not needed in case the EndFileMacro was used)
	//TFolder *t = (TFolder*)f.Get("ROOT Memory");

	// Get the Prompt-Random-Ratio from the analysed root file
	TH1F *h_ratio = (TH1F*)f.FindObjectAny("PHYS_promptRandomRatio");
	const double ratio = h_ratio->GetMean();
	std::cout << "[INFO] Prompt-Random-Ratio: " << ratio << std::endl;

	// declare some temporary used variables
	TH1F* h[N_WINDOWS];
	char buffer[50];
	// create an empty canvas for histogram drawing
	TCanvas *c = new TCanvas("c", "plot", 20, 10, 700, 500);
	c->SetFillColor(0);
	c->SetBorderMode(0);
	c->SetBorderSize(2);
	//c->SetLogy();
	c->SetFrameBorderMode(0);
	c->SetLeftMargin(.14);
	c->SetRightMargin(.1);
	c->SetBottomMargin(.125);
	c->SetTopMargin(.05);

	for (std::vector<const char*>::iterator it_prop = props.begin(); it_prop != props.end(); ++it_prop) {
		for (ICIter it_cut = cuts.begin(); it_cut != cuts.end(); ++it_cut) {
			// skip the current combination if the cut contains the current variable, e. g. invM with some cut involving invM
			if (strstr(it_cut->second, *it_prop))
				continue;
			// get prompt and random histograms for every property cut combination
			for (unsigned int i = 0; i < N_WINDOWS; i++) {
				sprintf(buffer, "%s_%s_%s", *it_prop, it_cut->second, windows[i]);
				h[i] = (TH1F*)f.FindObjectAny(buffer);
				if (!h[i]) {
					std::cout << "[ERROR] Couldn't find histogram " << buffer << std::endl;
					exit(1);
				}
			}
			c->Clear();
			// subtract to chosen window size normalized random-background from events inside the prompt-window
			h[PROMPT]->Add(h[RANDOM], -1*ratio);
			prepare_hist(h[PROMPT]);
			/* Rebin histogram */
			int n = 50;
			h[PROMPT]->Rebin(n);  // combine n bins => n MeV binning
			/*char* label;
			if (label)
				free(label);
			char* label = malloc(strlen("#Events / ") + 5 + strlen(" MeV^{2}") + 1);
			//sprintf(label, "%s%.2f%s", "#Events / ", n/1000., " GeV^{2}");
			sprintf(label, "%s%d%s", "#Events / ", n, " MeV^{2}");
			h->GetYaxis()->SetTitle(label);*/
			/* set axes according to plotted values and prepare histogram for saving */
			sprintf(buffer, "%s%d%s", "#Events / ", n, " MeV");
			h[PROMPT]->GetYaxis()->SetTitle(buffer);
			h[PROMPT]->SetLineColor(color);
			sprintf(buffer, "%s [MeV]", *it_prop);
			h[PROMPT]->GetXaxis()->SetTitle(buffer);
			h[PROMPT]->Draw();
			sprintf(buffer, "%s_%s.%s", *it_prop, it_cut->second, suffix);
			c->Print(buffer);
		}
		std::cout << "[INFO] Created all random-subtracted plots for " << *it_prop << std::endl;
	}

	f.Close();

	return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	gROOT->Reset();

	if (argc > 2) {
		printf("Error: Too many arguments!\n");
		return argc;
	}

	if (argc == 2) {  // given parameter is the file name
		if (!strstr(argv[1], ".root")) {  // parameter doesn't contain ".root"
			printf("Error: Input %s is not a root file!\n", argv[1]);
			return -1;
		}

		return create_plots(argv[1]);
	}

	return create_plots();
}
