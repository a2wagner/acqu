/**
 * TA2SaschaPhysics
 *
 * Physics class for the analysis of eta' physics, mainly the Dalitz decay mode. 
 * This class is designed to analyze Geant4 generated data as well as data taken in experiments. 
 */

#include "TA2SaschaPhysics.h"

enum { EPromptWindows = 1000, ERandomWindows};

static const Map_t kPhysics[] = {
	{"Prompt:",	EPromptWindows},
	{"Random:",	ERandomWindows},
	{NULL, -1}
};

ClassImp(TA2SaschaPhysics)

//-----------------------------------------------------------------------------

TA2SaschaPhysics::TA2SaschaPhysics(const char* Name, TA2Analysis* Analysis) : TA2BasePhysics(Name, Analysis)
{
	// Initialise prompt and random windows if they are not defined in the config file
	promptLow = -6;
	promptHigh = 4;
	random1Low = -40;
	random1High = -10;
	random2Low = 6;
	random2High = 36;
/* doesn't make any difference when initialised with NULL
	timeLeptons = NULL;
	timeTagger = NULL;
	dTimeFS = NULL;

	for (int i = 0; i < N_WINDOWS; i++) {
		invM_protE[i] = NULL;
		invM_copl[i] = NULL;
		invM_balance[i] = NULL;
		invM_dAlphaProtTAPS[i] = NULL;
		invM_missM[i] = NULL;
		invM_copl_balance[i] = NULL;
		invM_balance_missM[i] = NULL;
		invM_copl_missM[i] = NULL;
		invM_balance_dAlpha[i] = NULL;
		invM_allCuts[i] = NULL;
		balancePx[i] = NULL;
		balancePy[i] = NULL;
		balancePz[i] = NULL;
		balanceE[i] = NULL;
		missM[i] = NULL;
		protonEnergyExpect[i] = NULL;
		protDAlphaTAPSCl[i] = NULL;
	}
*/
	AddCmdList(kPhysics);  // enables keyword recognition in SetConfig
}

//-----------------------------------------------------------------------------

TA2SaschaPhysics::~TA2SaschaPhysics()
{
	delete KinFit;

	// Arrays concerning simulation
	delete[] trueIDPart;
	delete[] trueP4Gamma;
	delete[] trueP4Elec;
	delete[] trueP4Posi;
	delete[] trueP4Prot;
	delete[] trueP4MuPls;
	delete[] trueP4MuMns;
	delete[] trueP4PiPls;
	delete[] trueP4PiMns;
	delete[] trueProtMass;
	delete[] trueProtEnergy;
	delete[] trueProtTheta;
	delete[] trueProtPhi;
	delete[] trueGammaMass;
	delete[] trueGammaEnergy;
	delete[] trueGammaTheta;
	delete[] trueGammaPhi;
	delete[] truePosiMass;
	delete[] truePosiEnergy;
	delete[] truePosiTheta;
	delete[] truePosiPhi;
	delete[] trueElecMass;
	delete[] trueElecEnergy;
	delete[] trueElecTheta;
	delete[] trueElecPhi;
	delete[] trueMuPlsMass;
	delete[] trueMuPlsEnergy;
	delete[] trueMuPlsTheta;
	delete[] trueMuPlsPhi;
	delete[] trueMuMnsMass;
	delete[] trueMuMnsEnergy;
	delete[] trueMuMnsTheta;
	delete[] trueMuMnsPhi;
	delete[] truePiPlsMass;
	delete[] truePiPlsEnergy;
	delete[] truePiPlsTheta;
	delete[] truePiPlsPhi;
	delete[] truePiMnsMass;
	delete[] truePiMnsEnergy;
	delete[] truePiMnsTheta;
	delete[] truePiMnsPhi;

	delete[] timeLeptons;
	delete[] timeTagger;
	delete[] dTimeFS;

    delete[] particles;
}

//---------------------------------------------------------------------------

void TA2SaschaPhysics::SetConfig(Char_t* line, Int_t key)
{
	/* read prompt and random window times */
	switch (key) {
	case EPromptWindows:
		if (sscanf(line, "%lf %lf\n", &promptLow, &promptHigh) != 2) {
			PrintError(line, "<Error: Prompt windows not set correctly!>");
			return;
		} else
			printf("Prompt window set from %f to %f\n", promptLow, promptHigh);
		break;
	case ERandomWindows:
		if (sscanf(line, "%lf %lf %lf %lf\n", &random1Low, &random1High, &random2Low, &random2High) != 4) {
			PrintError(line, "<Error: Random windows not set correctly!>");
			return;
		} else
			printf("Random windows set from %f to %f and from %f to %f\n", random1Low, random1High, random2Low, random2High);
		break;
	default:
		// default SetConfig()
		TA2BasePhysics::SetConfig(line, key);
		break;
	}
}

//---------------------------------------------------------------------------

void TA2SaschaPhysics::LoadVariable()
{
    //Call default LoadVariable()
    TA2BasePhysics::LoadVariable();

	// number of particles
	TA2DataManager::LoadVariable("nParticles", &nParticles, EISingleX);
	TA2DataManager::LoadVariable("nParticlesCB", &nParticlesCB, EISingleX);
	TA2DataManager::LoadVariable("nParticlesTAPS", &nParticlesTAPS, EISingleX);

	// True information from the simulation
	if (MC) {  // ugly, but using this switch all the true histograms aren'T used while analysing data
		TA2DataManager::LoadVariable("trueNPart", &trueNPart, EISingleX);
		TA2DataManager::LoadVariable("trueIDPart", trueIDPart, EIMultiX);
		TA2DataManager::LoadVariable("trueNProt", &trueNProt, EISingleX);
		TA2DataManager::LoadVariable("trueProtMass", trueProtMass, EDMultiX);
		TA2DataManager::LoadVariable("trueProtEnergy", trueProtEnergy, EDMultiX);
		TA2DataManager::LoadVariable("trueProtTheta", trueProtTheta, EDMultiX);
		TA2DataManager::LoadVariable("trueProtPhi", trueProtPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNGamma", &trueNGamma, EISingleX);
		TA2DataManager::LoadVariable("trueGammaMass", trueGammaMass, EDMultiX);
		TA2DataManager::LoadVariable("trueGammaEnergy", trueGammaEnergy, EDMultiX);
		TA2DataManager::LoadVariable("trueGammaTheta", trueGammaTheta, EDMultiX);
		TA2DataManager::LoadVariable("trueGammaPhi", trueGammaPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNPosi", &trueNPosi, EISingleX);
		TA2DataManager::LoadVariable("truePosiMass", truePosiMass, EDMultiX);
		TA2DataManager::LoadVariable("truePosiEnergy", truePosiEnergy, EDMultiX);
		TA2DataManager::LoadVariable("truePosiTheta", truePosiTheta, EDMultiX);
		TA2DataManager::LoadVariable("truePosiPhi", truePosiPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNElec", &trueNElec, EISingleX);
		TA2DataManager::LoadVariable("trueElecMass", trueElecMass, EDMultiX);
		TA2DataManager::LoadVariable("trueElecEnergy", trueElecEnergy, EDMultiX);
		TA2DataManager::LoadVariable("trueElecTheta", trueElecTheta, EDMultiX);
		TA2DataManager::LoadVariable("trueElecPhi", trueElecPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNMuPls", &trueNMuPls, EISingleX);
		TA2DataManager::LoadVariable("trueMuPlsMass", trueMuPlsMass, EDMultiX);
		TA2DataManager::LoadVariable("trueMuPlsEnergy", trueMuPlsEnergy, EDMultiX);
		TA2DataManager::LoadVariable("trueMuPlsTheta", trueMuPlsTheta, EDMultiX);
		TA2DataManager::LoadVariable("trueMuPlsPhi", trueMuPlsPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNMuMns", &trueNMuMns, EISingleX);
		TA2DataManager::LoadVariable("trueMuMnsMass", trueMuMnsMass, EDMultiX);
		TA2DataManager::LoadVariable("trueMuMnsEnergy", trueMuMnsEnergy, EDMultiX);
		TA2DataManager::LoadVariable("trueMuMnsTheta", trueMuMnsTheta, EDMultiX);
		TA2DataManager::LoadVariable("trueMuMnsPhi", trueMuMnsPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNPiPls", &trueNPiPls, EISingleX);
		TA2DataManager::LoadVariable("truePiPlsMass", truePiPlsMass, EDMultiX);
		TA2DataManager::LoadVariable("truePiPlsEnergy", truePiPlsEnergy, EDMultiX);
		TA2DataManager::LoadVariable("truePiPlsTheta", truePiPlsTheta, EDMultiX);
		TA2DataManager::LoadVariable("truePiPlsPhi", truePiPlsPhi, EDMultiX);
		TA2DataManager::LoadVariable("trueNPiMns", &trueNPiMns, EISingleX);
		TA2DataManager::LoadVariable("truePiMnsMass", truePiMnsMass, EDMultiX);
		TA2DataManager::LoadVariable("truePiMnsEnergy", truePiMnsEnergy, EDMultiX);
		TA2DataManager::LoadVariable("truePiMnsTheta", truePiMnsTheta, EDMultiX);
		TA2DataManager::LoadVariable("truePiMnsPhi", truePiMnsPhi, EDMultiX);
	}

	TA2DataManager::LoadVariable("promptRandomRatio", &promptRandomRatio, EDSingleX);

	TA2DataManager::LoadVariable("invM_allPart", &invM_allPart, EDSingleX);
	TA2DataManager::LoadVariable("invM_2charged", &invM_2charged, EDSingleX);
	TA2DataManager::LoadVariable("timeProton", &timeProton, EDSingleX);
	TA2DataManager::LoadVariable("timePhoton", &timePhoton, EDSingleX);
	TA2DataManager::LoadVariable("timeFS", &timeFS, EDSingleX);
	TA2DataManager::LoadVariable("timeLeptons", timeLeptons, EDMultiX);
	TA2DataManager::LoadVariable("timeTagger", timeTagger, EDMultiX);
	TA2DataManager::LoadVariable("dTimeFS", dTimeFS, EDMultiX);
	TA2DataManager::LoadVariable("coplanarity", &coplanarity, EDSingleX);
	TA2DataManager::LoadVariable("protEnergyReconstr", &protEnergyReconstr, EDSingleX);

	TA2DataManager::LoadVariable("invMass2g", &invMass2g, EDSingleX);
	TA2DataManager::LoadVariable("invMass2g1p", &invMass2g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g", &invMass6g, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g1p", &invMass6g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB", &invMass2CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB1TAPS", &invMass2CB1TAPS, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB", &invMass6CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB1TAPS", &invMass6CB1TAPS, EDSingleX);

	// Prompt
	TA2DataManager::LoadVariable("invM_protE_prompt", invM_protE[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_prompt", invM_copl[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_prompt", invM_balance[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_dAlphaProtTAPS_prompt", invM_dAlphaProtTAPS[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_missM_prompt", invM_missM[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_balance_prompt", invM_copl_balance[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_missM_prompt", invM_balance_missM[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_missM_prompt", invM_copl_missM[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_dAlpha_prompt", invM_balance_dAlpha[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("invM_allCuts_prompt", invM_allCuts[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("balancePx_prompt", balancePx[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("balancePy_prompt", balancePy[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("balancePz_prompt", balancePz[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("balanceE_prompt", balanceE[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("missM_prompt", missM[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("protonEnergyExpect_prompt", protonEnergyExpect[PROMPT], EDMultiX);
	TA2DataManager::LoadVariable("protDAlphaTAPSCl_prompt", protDAlphaTAPSCl[PROMPT], EDMultiX);
	// Random
	TA2DataManager::LoadVariable("invM_protE_random", invM_protE[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_random", invM_copl[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_random", invM_balance[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_dAlphaProtTAPS_random", invM_dAlphaProtTAPS[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_missM_random", invM_missM[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_balance_random", invM_copl_balance[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_missM_random", invM_balance_missM[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_copl_missM_random", invM_copl_missM[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_balance_dAlpha_random", invM_balance_dAlpha[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("invM_allCuts_random", invM_allCuts[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("balancePx_random", balancePx[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("balancePy_random", balancePy[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("balancePz_random", balancePz[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("balanceE_random", balanceE[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("missM_random", missM[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("protonEnergyExpect_random", protonEnergyExpect[RANDOM], EDMultiX);
	TA2DataManager::LoadVariable("protDAlphaTAPSCl_random", protDAlphaTAPSCl[RANDOM], EDMultiX);
}

//---------------------------------------------------------------------------

void TA2SaschaPhysics::PostInit()
{
	std::cout << std::endl;  // insert an empty line
	// Check whether simulated or measured data will be processed
	if (gAR->GetProcessType() == EMCProcess) {
		MC = true;
		std::cout << "Analysis of Monte Carlo data" << std::endl;
	} else {
		MC = false;
		std::cout << "Analysis of experimental data" << std::endl;
	}

	/* Get pointers to the different apparati */
	// Tagger
	fTAG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("EPT");
	if (!fTAG) {
		fTAG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
		if (fTAG)
			std::cout << "Found Tagger" << std::endl;
		else
			PrintError("", "No Tagger found!", EErrFatal);
	} else
		std::cout << "Found EPT" << std::endl;

	// Crystal Ball
	fCB = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");
	if (fCB)
		std::cout << "Found Crystal Ball" << std::endl;
	else
		PrintError("", "No Crystal Ball found!", EErrFatal);

	// TAPS
	fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
	if (fTAPS)
		std::cout << "Found TAPS" << std::endl;
	else
		std::cout << "TAPS will not be included in the analysis" << std::endl;

	// Get maximum particles
	maxBeamPhotons = fTAG->GetMaxParticle();
	maxParticlesCB = fCB->GetMaxParticle();
	maxParticlesTAPS = fTAPS ? fTAPS->GetMaxParticle() : 0;
	maxParticles = maxParticlesCB + maxParticlesTAPS;
	// in principle the values nBeam and nPart from TA2BasePhysics could be used...

	if (MC)
		promptRandomRatio = 0.;
	else {
		promptRandomRatio = (promptHigh - promptLow)/(random1High - random1Low + random2High - random2Low);
		std::cout << "Prompt-Random-Ratio: " << promptRandomRatio << std::endl;
	}


	// Information from the simulation
	trueP4Gamma = new TLorentzVector[maxParticles];
	trueP4Elec = new TLorentzVector[maxParticles];
	trueP4Posi = new TLorentzVector[maxParticles];
	trueP4Prot = new TLorentzVector[maxParticles];
	trueP4MuPls = new TLorentzVector[maxParticles];
	trueP4MuMns = new TLorentzVector[maxParticles];
	trueP4PiPls = new TLorentzVector[maxParticles];
	trueP4PiMns = new TLorentzVector[maxParticles];
	trueIDPart = new Int_t[maxParticles];
	trueProtMass = new Double_t[maxParticles];
	trueProtEnergy = new Double_t[maxParticles];
	trueProtTheta = new Double_t[maxParticles];
	trueProtPhi = new Double_t[maxParticles];
	trueGammaMass = new Double_t[maxParticles];
	trueGammaEnergy = new Double_t[maxParticles];
	trueGammaTheta = new Double_t[maxParticles];
	trueGammaPhi = new Double_t[maxParticles];
	truePosiMass = new Double_t[maxParticles];
	truePosiEnergy = new Double_t[maxParticles];
	truePosiTheta = new Double_t[maxParticles];
	truePosiPhi = new Double_t[maxParticles];
	trueElecMass = new Double_t[maxParticles];
	trueElecEnergy = new Double_t[maxParticles];
	trueElecTheta = new Double_t[maxParticles];
	trueElecPhi = new Double_t[maxParticles];
	trueMuPlsMass = new Double_t[maxParticles];
	trueMuPlsEnergy = new Double_t[maxParticles];
	trueMuPlsTheta = new Double_t[maxParticles];
	trueMuPlsPhi = new Double_t[maxParticles];
	trueMuMnsMass = new Double_t[maxParticles];
	trueMuMnsEnergy = new Double_t[maxParticles];
	trueMuMnsTheta = new Double_t[maxParticles];
	trueMuMnsPhi = new Double_t[maxParticles];
	truePiPlsMass = new Double_t[maxParticles];
	truePiPlsEnergy = new Double_t[maxParticles];
	truePiPlsTheta = new Double_t[maxParticles];
	truePiPlsPhi = new Double_t[maxParticles];
	truePiMnsMass = new Double_t[maxParticles];
	truePiMnsEnergy = new Double_t[maxParticles];
	truePiMnsTheta = new Double_t[maxParticles];
	truePiMnsPhi = new Double_t[maxParticles];

    particles = new TA2Particle[maxParticles];

	timeLeptons = new Double_t[2+1];
	timeTagger = new Double_t[maxBeamPhotons];
	dTimeFS = new Double_t[maxBeamPhotons];

	for (int i = 0; i < N_WINDOWS; i++) {
		invM_protE[i] = new Double_t[maxBeamPhotons];
		invM_copl[i] = new Double_t[maxBeamPhotons];
		invM_balance[i] = new Double_t[maxBeamPhotons];
		invM_dAlphaProtTAPS[i] = new Double_t[maxBeamPhotons];
		invM_missM[i] = new Double_t[maxBeamPhotons];
		invM_copl_balance[i] = new Double_t[maxBeamPhotons];
		invM_balance_missM[i] = new Double_t[maxBeamPhotons];
		invM_copl_missM[i] = new Double_t[maxBeamPhotons];
		invM_balance_dAlpha[i] = new Double_t[maxBeamPhotons];
		invM_allCuts[i] = new Double_t[maxBeamPhotons];
		balancePx[i] = new Double_t[maxBeamPhotons];
		balancePy[i] = new Double_t[maxBeamPhotons];
		balancePz[i] = new Double_t[maxBeamPhotons];
		balanceE[i] = new Double_t[maxBeamPhotons];
		missM[i] = new Double_t[maxBeamPhotons];
		protonEnergyExpect[i] = new Double_t[maxBeamPhotons];
		protDAlphaTAPSCl[i] = new Double_t[maxBeamPhotons];
	}

	KinFit = new TA2CBKinematicFitter(3, 1, 0);

	TFile cutFile("/home/wagners/acqu/acqu_user/data/myCuts.root", "READ");
	if (!cutFile.IsOpen()) {
		printf("Error opening the file containing the cuts!\n");
		exit(1);
	}
    //TODO: seems to lead to a crash when not used as TCutG balance = ... ?
    //cutBalance = *(TCutG*)cutFile.Get("balanceCut");
	cutFile.Close();

	// just for testing
	//timeTagger now works, so this is obsolete
	taggerTime = new TH1F("taggerTime", "taggerTime", 2000, -500, 500);

	// Call default PostInit()
	// Important!: Call this AFTER the previous arrays are defined, otherwise LoadVariable gets called and all arrays are still NULL
	TA2BasePhysics::PostInit();
}

//-----------------------------------------------------------------------------

/* Main method for applying the physics related stuff */
void TA2SaschaPhysics::Reconstruct()
{
	//Perform basic physics tasks
	TA2BasePhysics::Reconstruct();

	bool dbg = false;  // print some output for debugging if set to true

	//Initialise array counters
	VarInit();
	  
	Char_t currentFile[1024];
	gUAN->ReadRunName(currentFile);  // save filename without path and ending in the specified char array

	if (dbg)
		std::cout << "Processing file " << currentFile << std::endl;

	// Get true information from simulation
	//if (strstr(currentFile, "g4_sim")) {  // if the filename contains "g4_sim" than it is one of my generated Monte Carlo event files
	if (MC) {
		GetTrueBeam();
		SetTrueParticles();
		PlotTrueParticles();
	}

	/* Get the number of particles from the detectors */
	nBeamPhotons = fTAG->GetNparticle();
	nParticlesCB = fCB->GetNparticle();
	nParticlesTAPS = fTAPS ? fTAPS->GetNparticle() : 0;
	nParticles = nParticlesCB + nParticlesTAPS;

	if (dbg) {
		std::cout << "Found " << nParticlesCB << " particles in CB";
		if (fTAPS)
			std::cout << " and " << nParticlesTAPS << " particles in TAPS";
		std::cout << ", total " << nParticles << " particles; " << nBeamPhotons << " photons in Tagger" << std::endl;
	}

	/* Collect all particles from CB and TAPS */
	/*int n = 0;
	for (UInt_t i = 0; i < nParticlesCB; i++)
		particles[n++] = fCB->GetParticles(i);
	for (UInt_t i = 0; i < nParticlesTAPS; i++)
		particles[n++] = fTAPS->GetParticles(i);*/

	/* Collect all particles from the different detectors in the particles array */
	if (MC)
		GetTrueParticles();
	else
		GetParticles();

	if (dbg) {
		std::cout << "Gathered " << nParticlesCB << " particles in CB";
		if (fTAPS)
			std::cout << " and " << nParticlesTAPS << " particles in TAPS";
		std::cout << ", total " << nParticles << " particles for further analysis" << std::endl;
	}

	bool cutsPassed = ApplyCuts();
	//bool kinFit = KinematicFit();

	if (cutsPassed /*|| kinFit*/)
		PlotData();


	// fast check for multiple gamma/CB events with/out proton/TAPS-Cluster
	int nPhoton = 0;
	bool hasProton = false;
	TA2Particle photons[maxParticles];
	//if (nParticles == 3)
	for (int i = 0; i < nParticles; i++) {
		if (particles[i].GetParticleID() == kGamma)
			photons[nPhoton++] = particles[i];
		else if (particles[i].GetParticleID() == kProton && !hasProton)
			hasProton = true;
	}

	// now sort the events
	if (nPhoton == 2) {
		invMass2g = (photons[0].GetP4()+photons[1].GetP4()).M();
		if (hasProton)
			invMass2g1p = invMass2g;
	}
	if (nPhoton == 6) {
		TLorentzVector tmp(0., 0., 0., 0.);
		for (int i = 0; i < 6; i++)
			tmp += photons[i].GetP4();
		invMass6g = tmp.M();
		if (hasProton)
			invMass6g1p = invMass6g;
	}
	if (nParticlesCB == 2) {
		invMass2CB = (particles[0].GetP4()+particles[1].GetP4()).M();
		if (nParticlesTAPS == 1)
			invMass2CB1TAPS = invMass2CB;
	}
	if (nParticlesCB == 6) {
		TLorentzVector tmp(0., 0., 0., 0.);
		for (int i = 0; i < 6; i++)
			tmp += particles[i].GetP4();
		invMass6CB = tmp.M();
		if (nParticlesTAPS == 1)
			invMass6CB1TAPS = invMass6CB;
	}

	//Clean up and set EBufferEnd end markers for arrays
	TermArrays();
}

//-----------------------------------------------------------------------------

/* Collect all particles from data in one array */
void TA2SaschaPhysics::GetParticles()
{
	nParticles = 0;
	nParticlesCB = 0;
	nParticlesTAPS = 0;

	if (!trigger) return;  // skip getting particles of this event if it's not triggered

	// Collect detected particles from CB
	for (Int_t nCB = 0; nCB < fCB->GetNparticle(); nCB++) {
		//if (fCB->GetParticles(nCB).GetParticleID() != kProton) {  // If it shouldn't be a proton...
		particles[nParticles++] = fCB->GetParticles(nCB);  // copy to TA2Particle array
		nParticlesCB++;
	}

	// Collect detected particles from TAPS, if available
	if (fTAPS)
		for (Int_t nTAPS = 0; nTAPS < fTAPS->GetNparticle(); nTAPS++) {
		//if (fTAPS->GetParticles(nTAPS).GetParticleID() != kProton) {  // If it shouldn't be a proton...
			particles[nParticles++] = fTAPS->GetParticles(nTAPS);  // copy to TA2Particle array
			nParticlesTAPS++;
		}
}

/* Fill true information of all particles in particles array, for some algorithm testing */
void TA2SaschaPhysics::GetTrueParticles()
{
	nParticles = 0;
	nParticlesCB = 0;
	nParticlesTAPS = 0;

	for (Int_t i = 0; i < trueNGamma; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4Gamma[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNPosi; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4Posi[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNElec; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4Elec[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNMuPls; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4MuPls[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNMuMns; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4MuMns[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNPiPls; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4PiPls[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNPiMns; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4PiMns[i]);
		nParticlesCB++;
	}
	for (Int_t i = 0; i < trueNProt; i++) {
		particles[nParticles] = TA2Particle();
		particles[nParticles++].SetP4(trueP4Prot[i]);
		nParticlesTAPS++;
	}
}

/* Apply the kinematic cuts to suppress the background */
bool TA2SaschaPhysics::ApplyCuts()
{
	bool success = false;  // bool to indicate if the event passed all cuts
	// bools concerning the cuts
	passedProtonEnergy = false;
	passedCoplanarity = false;
	passedBalance = false;
	passedDAlphaProtTAPS = false;
	passedMissMass = false;
	passedInvMass = false;
	passedAllCuts = false;

	// reset all counters
	nPrompt = nRandom = 0;
	for (int i = 0; i < N_WINDOWS; i++) {
		n_protE[i] = 0;
		n_copl[i] = 0;
		n_balance[i] = 0;
		n_dAlphaProtTAPS[i] = 0;
		n_missM[i] = 0;
		n_copl_balance[i] = 0;
		n_balance_missM[i] = 0;
		n_copl_missM[i] = 0;
		n_balance_dAlpha[i] = 0;
		n_allCuts[i] = 0;
	}

	if (!trigger && MC)
		return success;  // return success (which is false until here) if MC event was not triggered, also skip everything

	TLorentzVector tmpState(0., 0., 0., 0.);
	TLorentzVector chargedPart[N_FINAL_STATE];
	TLorentzVector balance, missProton;
	// variables to store values only for cuts
	Double_t mMiss, protEexpect, dAlphaProtTAPS;

	if (nParticlesCB == 3 && nParticlesTAPS == 1) {
		// if there are 3 particles in CB (e+, e-, g) and one in TAPS (proton), apply the cuts
		for (int i = 0; i < N_FINAL_STATE; i++)
			tmpState += particles[i].GetP4();
		invM_allPart = tmpState.M();
		timeProton = particles[nParticlesCB].GetTime();

		nCharged = 0;
		for (int i = 0; i < N_FINAL_STATE; i++) {
			if (particles[i].HasDetector(EDetPID)) {
				chargedPart[nCharged++] = particles[i].GetP4();
				if (nCharged <= 2)
					timeLeptons[nCharged-1] = particles[i].GetTime();
			} else
				timePhoton = particles[i].GetTime();
		}
		/* end event selection when we don't have two charged particles (entry in PID) */
		if (nCharged != 2)
			return success;

		timeFS = (timeLeptons[0] + timeLeptons[1] + timePhoton)/3.;
		coplanarity = (tmpState.Phi() - particles[nParticlesCB].GetP4().Phi())*R2D;
		coplanarity = TMath::Abs(coplanarity);
		protEnergyReconstr = particles[nParticlesCB].GetP4().E() - particles[nParticlesCB].GetP4().M();
		invM_2charged = tmpState.M();
		passedInvMass = true;//invM_2charged > 900. && invM_2charged < 1000. ? true : false;

		/* At this point we have two charged and one neutral particle in the CB as well as one particle in TAPS, probably the proton. Start to examine all tagged photons for prompt and random windows. */
		for (int i = 0; i < nBeamPhotons; i++) {
			timeTagger[i] = Tagged[i].GetTime();
			// for testing, as EDMultiX doesn't do anything...
			taggerTime->Fill(Tagged[i].GetTime());
			dTimeFS[i] = timeTagger[i] - timeFS;
			balance = tmpState + particles[nParticlesCB].GetP4() - fP4target[0] - Tagged[i].GetP4();
			missProton = fP4target[0] + Tagged[i].GetP4() - tmpState;
			dAlphaProtTAPS = particles[nParticlesCB].GetP4().Angle(missProton.Vect())*R2D;  // opening angle between TAPS cluster and expected proton
			mMiss = missProton.M();
			protEexpect = missProton.E() - missProton.M();

			/* now "apply" the kinematic cuts */
			passedProtonEnergy = protEexpect < 350. ? true : false;
			passedCoplanarity = coplanarity > 160. && coplanarity < 200 ? true : false;
			passedBalance = true;//cutBalance.IsInside(balance.E(), balance.Pz()) ? true : false;
			passedDAlphaProtTAPS = dAlphaProtTAPS < 6. ? true : false;
			passedMissMass = true;//mMiss > 900. && mMiss < 980. ? true : false;
			if (passedProtonEnergy && passedCoplanarity && passedBalance && passedDAlphaProtTAPS && passedMissMass && passedInvMass)
				passedAllCuts = true;
			else
				passedAllCuts = false;

//		the different cuts that can be used
//			if (passedProtonEnergy)
//			if (passedCoplanarity)
//			if (passedBalance)
//			if (passedDAlphaProtTAPS)
//			if (passedMissMass)
//			if (passedInvMass)
//			if (passedAllCuts)

			/* first fill the prompt histograms */
			if ((timeTagger[i] >= promptLow && timeTagger[i] <= promptHigh) || MC) {  // if MC events are processed then fill them in the prompt array
				missM[PROMPT][nPrompt] = mMiss;
				balancePx[PROMPT][nPrompt] = balance.Px();
				balancePy[PROMPT][nPrompt] = balance.Py();
				balancePz[PROMPT][nPrompt] = balance.Pz();
				balanceE[PROMPT][nPrompt] = balance.E();
				protonEnergyExpect[PROMPT][nPrompt] = protEexpect;
				protDAlphaTAPSCl[PROMPT][nPrompt] = dAlphaProtTAPS;
				nPrompt++;

				if (passedProtonEnergy) {
					invM_protE[PROMPT][n_protE[PROMPT]] = invM_2charged;
					n_protE[PROMPT]++;
				}
				if (passedCoplanarity) {
					invM_copl[PROMPT][n_copl[PROMPT]] = invM_2charged;
					n_copl[PROMPT]++;
				}
				if (passedBalance) {
					invM_balance[PROMPT][n_balance[PROMPT]] = invM_2charged;
					n_balance[PROMPT]++;
				}
				if (passedDAlphaProtTAPS) {
					invM_dAlphaProtTAPS[PROMPT][n_dAlphaProtTAPS[PROMPT]] = invM_2charged;
					n_dAlphaProtTAPS[PROMPT]++;
				}
				if (passedMissMass) {
					invM_missM[PROMPT][n_missM[PROMPT]] = invM_2charged;
					n_missM[PROMPT]++;
				}
				if (passedCoplanarity && passedBalance) {
					invM_copl_balance[PROMPT][n_copl_balance[PROMPT]] = invM_2charged;
					n_copl_balance[PROMPT]++;
				}
				if (passedBalance && passedMissMass) {
					invM_balance_missM[PROMPT][n_balance_missM[PROMPT]] = invM_2charged;
					n_balance_missM[PROMPT]++;
				}
				if (passedCoplanarity && passedMissMass) {
					invM_copl_missM[PROMPT][n_copl_missM[PROMPT]] = invM_2charged;
					n_copl_missM[PROMPT]++;
				}
				if (passedBalance && passedDAlphaProtTAPS) {
					invM_balance_dAlpha[PROMPT][n_balance_dAlpha[PROMPT]] = invM_2charged;
					n_balance_dAlpha[PROMPT]++;
				}
				if (passedAllCuts) {
					invM_allCuts[PROMPT][n_allCuts[PROMPT]] = invM_2charged;
					n_allCuts[PROMPT]++;

					// if all cuts in the prompt window were passed, set success to true which is returned by this method
					if (success)
						std::cout << "[WARNING] More than one tagged photon combination matched all cuts in current event" << std::endl;
					success = true;
				}
			}

			/* now the random ones */
			if ((timeTagger[i] >= random1Low && timeTagger[i] <= random1High) || (timeTagger[i] >= random2Low && timeTagger[i] <= random2High)) {
				missM[RANDOM][nRandom] = mMiss;
				balancePx[RANDOM][nRandom] = balance.Px();
				balancePy[RANDOM][nRandom] = balance.Py();
				balancePz[RANDOM][nRandom] = balance.Pz();
				balanceE[RANDOM][nRandom] = balance.E();
				protonEnergyExpect[RANDOM][nRandom] = protEexpect;
				protDAlphaTAPSCl[RANDOM][nRandom] = dAlphaProtTAPS;
				nRandom++;

				if (passedProtonEnergy) {
					invM_protE[RANDOM][n_protE[RANDOM]] = invM_2charged;
					n_protE[RANDOM]++;
				}
				if (passedCoplanarity) {
					invM_copl[RANDOM][n_copl[RANDOM]] = invM_2charged;
					n_copl[RANDOM]++;
				}
				if (passedBalance) {
					invM_balance[RANDOM][n_balance[RANDOM]] = invM_2charged;
					n_balance[RANDOM]++;
				}
				if (passedDAlphaProtTAPS) {
					invM_dAlphaProtTAPS[RANDOM][n_dAlphaProtTAPS[RANDOM]] = invM_2charged;
					n_dAlphaProtTAPS[RANDOM]++;
				}
				if (passedMissMass) {
					invM_missM[RANDOM][n_missM[RANDOM]] = invM_2charged;
					n_missM[RANDOM]++;
				}
				if (passedCoplanarity && passedBalance) {
					invM_copl_balance[RANDOM][n_copl_balance[RANDOM]] = invM_2charged;
					n_copl_balance[RANDOM]++;
				}
				if (passedBalance && passedMissMass) {
					invM_balance_missM[RANDOM][n_balance_missM[RANDOM]] = invM_2charged;
					n_balance_missM[RANDOM]++;
				}
				if (passedCoplanarity && passedMissMass) {
					invM_copl_missM[RANDOM][n_copl_missM[RANDOM]] = invM_2charged;
					n_copl_missM[RANDOM]++;
				}
				if (passedBalance && passedDAlphaProtTAPS) {
					invM_balance_dAlpha[RANDOM][n_balance_dAlpha[RANDOM]] = invM_2charged;
					n_balance_dAlpha[RANDOM]++;
				}
				if (passedAllCuts) {
					invM_allCuts[RANDOM][n_allCuts[RANDOM]] = invM_2charged;
					n_allCuts[RANDOM]++;
				}
			}

			// if all cuts were passed, set success to true which is returned by this method
/*			if (passedAllCuts) {
				if (success)
					std::cout << "[WARNING] More than one tagged photon combination matched all cuts in current event" << std::endl;
				success = true;
			}*/
		}
	}

	return success;
}

/* Perform a kinematic fit using constraints like missing mass of the proton or energy conservation */
bool TA2SaschaPhysics::KinematicFit()
{
	bool success = false;  // bool to indicate if a good fit value was achieved

	//...

	return success;
}

/* Plot particle properties of the final events which passed the kinematic cuts or the kinematic fit */
void TA2SaschaPhysics::PlotData()
{
	//...
}

/* Do the analysis for all combinations possible with the detected particles */
/*bool TA2SaschaPhysics::ReconstructCombinatorics()
{
	return 0;
}*/

/**
* For the brave souls who get this far: You are the chosen ones,
* the valiant knights of programming who toil away, without rest,
* fixing our most awful code. To you, true saviors, kings of men,
* I say this: never gonna give you up, never gonna let you down,
* never gonna run around and desert you. Never gonna make you cry,
* never gonna say goodbye. Never gonna tell a lie and hurt you.
*/
int TA2SaschaPhysics::GetPermutations(Int_t perm[][N_FINAL_STATE])
{
	//std::cout << "Pointer:   " << &perm[0] << std::endl;
	if (nParticlesCB < 3) return 1;

	int n = 0, i = 0, j = 1, k = 2;

	while (i < j && i <= nParticlesCB-3) {
		while (j < k && j <= nParticlesCB-2) {
			while (k <= nParticlesCB-1) {
				// doesn't work anymore, c++11 not supported?
				//perm[n++] = {i, j, k++};
				perm[n][0] = i;
				perm[n][1] = j;
				perm[n++][2] = k++;
				//printf("%d %d %d (%d)\n", perm[n-1][0], perm[n-1][1], perm[n-1][2], n); //printf("%d %d %d (%d)\n", i, j, k++, ++n);
			}
			j++;
			k = j+1;
		}
		i++;
		j = i+1;
		k = j+1;
	}

	return 0;
}

/* Magic. Do not touch. */
unsigned int TA2SaschaPhysics::binomial(int n, int k)
{
	unsigned int e;
	if (k > n) return 0;
	if (k == 0 || k == n) return 1;
	if (2*k > n) return binomial(n, n-k);
	else {
		e = n-k+1;
		for (int i = 2; i <= k; i++) {
			e *= n-k+i;
			e /= i;
		}
	}

	return e;
}

//-----------------------------------------------------------------------------

void TA2SaschaPhysics::VarInit()
{
	//Initally, set single value histograms to EBufferEnd (this value will not be plotted)
	invM_allPart = EBufferEnd;
	invM_2charged = EBufferEnd;
	timeProton = EBufferEnd;
	timePhoton = EBufferEnd;
	timeFS = EBufferEnd;
	coplanarity = EBufferEnd;
	protEnergyReconstr = EBufferEnd;
	timeLeptons[0] = EBufferEnd;
	timeTagger[0] = EBufferEnd;
	dTimeFS[0] = EBufferEnd;

	invMass2g = EBufferEnd;
	invMass2g1p = EBufferEnd;
	invMass6g = EBufferEnd;
	invMass6g1p = EBufferEnd;
	invMass2CB = EBufferEnd;
	invMass2CB1TAPS = EBufferEnd;
	invMass6CB = EBufferEnd;
	invMass6CB1TAPS = EBufferEnd;

	for (int i = 0; i < N_WINDOWS; i++) {
		invM_protE[i][0] = EBufferEnd;
		invM_copl[i][0] = EBufferEnd;
		invM_balance[i][0] = EBufferEnd;
		invM_dAlphaProtTAPS[i][0] = EBufferEnd;
		invM_missM[i][0] = EBufferEnd;
		invM_copl_balance[i][0] = EBufferEnd;
		invM_balance_missM[i][0] = EBufferEnd;
		invM_copl_missM[i][0] = EBufferEnd;
		invM_balance_dAlpha[i][0] = EBufferEnd;
		invM_allCuts[i][0] = EBufferEnd;
		balancePx[i][0] = EBufferEnd;
		balancePy[i][0] = EBufferEnd;
		balancePz[i][0] = EBufferEnd;
		balanceE[i][0] = EBufferEnd;
		missM[i][0] = EBufferEnd;
		protonEnergyExpect[i][0] = EBufferEnd;
		protDAlphaTAPSCl[i][0] = EBufferEnd;
	}
}

//-----------------------------------------------------------------------------

/* Terminate all arrays with EBufferEnd markers at the end of the Reconstruct() method */
void TA2SaschaPhysics::TermArrays()
{
	// Simulation stuff
	trueIDPart[trueNPart] = EBufferEnd;
	trueProtMass[trueNProt] = EBufferEnd;
	trueProtEnergy[trueNProt] = EBufferEnd;
	trueProtTheta[trueNProt] = EBufferEnd;  
	trueProtPhi[trueNProt] = EBufferEnd;
	trueGammaMass[trueNGamma] = EBufferEnd;
	trueGammaEnergy[trueNGamma] = EBufferEnd;
	trueGammaTheta[trueNGamma] = EBufferEnd;  
	trueGammaPhi[trueNGamma] = EBufferEnd;
	truePosiMass[trueNPosi] = EBufferEnd;
	truePosiEnergy[trueNPosi] = EBufferEnd;
	truePosiTheta[trueNPosi] = EBufferEnd;  
	truePosiPhi[trueNPosi] = EBufferEnd;
	trueElecMass[trueNElec] = EBufferEnd;
	trueElecEnergy[trueNElec] = EBufferEnd;
	trueElecTheta[trueNElec] = EBufferEnd;  
	trueElecPhi[trueNElec] = EBufferEnd;
	trueMuPlsMass[trueNMuPls] = EBufferEnd;
	trueMuPlsEnergy[trueNMuPls] = EBufferEnd;
	trueMuPlsTheta[trueNMuPls] = EBufferEnd;
	trueMuPlsPhi[trueNMuPls] = EBufferEnd;
	trueMuMnsMass[trueNMuMns] = EBufferEnd;
	trueMuMnsEnergy[trueNMuMns] = EBufferEnd;
	trueMuMnsTheta[trueNMuMns] = EBufferEnd;
	trueMuMnsPhi[trueNMuMns] = EBufferEnd;
	truePiPlsMass[trueNPiPls] = EBufferEnd;
	truePiPlsEnergy[trueNPiPls] = EBufferEnd;
	truePiPlsTheta[trueNPiPls] = EBufferEnd;
	truePiPlsPhi[trueNPiPls] = EBufferEnd;
	truePiMnsMass[trueNPiMns] = EBufferEnd;
	truePiMnsEnergy[trueNPiMns] = EBufferEnd;
	truePiMnsTheta[trueNPiMns] = EBufferEnd;
	truePiMnsPhi[trueNPiMns] = EBufferEnd;

	timeLeptons[nCharged] = EBufferEnd;
	timeTagger[nBeamPhotons] = EBufferEnd;
	dTimeFS[nBeamPhotons] = EBufferEnd;

	for (int i = 0; i < N_WINDOWS; i++) {
		invM_protE[i][n_protE[i]] = EBufferEnd;
		invM_copl[i][n_copl[i]] = EBufferEnd;
		invM_balance[i][n_balance[i]] = EBufferEnd;
		invM_dAlphaProtTAPS[i][n_dAlphaProtTAPS[i]] = EBufferEnd;
		invM_missM[i][n_missM[i]] = EBufferEnd;
		invM_copl_balance[i][n_copl_balance[i]] = EBufferEnd;
		invM_balance_missM[i][n_balance_missM[i]] = EBufferEnd;
		invM_copl_missM[i][n_copl_missM[i]] = EBufferEnd;
		invM_balance_dAlpha[i][n_balance_dAlpha[i]] = EBufferEnd;
		invM_allCuts[i][n_allCuts[i]] = EBufferEnd;
		balancePx[i][nPrompt] = EBufferEnd;
		balancePy[i][nPrompt] = EBufferEnd;
		balancePz[i][nPrompt] = EBufferEnd;
		balanceE[i][nPrompt] = EBufferEnd;
		missM[i][nPrompt] = EBufferEnd;
		protonEnergyExpect[i][nPrompt] = EBufferEnd;
		protDAlphaTAPSCl[i][nPrompt] = EBufferEnd;
	}
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void TA2SaschaPhysics::ParseMisc(char* line)
{
	TA2BasePhysics::ParseMisc(line);
}

//-----------------------------------------------------------------------------

/* Reconstruct the beam information which were used in the Pluto event generator */
void TA2SaschaPhysics::GetTrueBeam()
{
	Float_t* beam = (Float_t*)(fEvent[EI_beam]);

	Double_t x, y, z;

	x = beam[0];
	y = beam[1];
	z = beam[2];
	trueP4Beam.SetT(beam[3]*1000.);
	trueP4Beam.SetX(beam[3]*x*1000.);
	trueP4Beam.SetY(beam[3]*y*1000.);
	trueP4Beam.SetZ(beam[3]*z*1000.);

	trueP4Target.SetXYZT(0., 0., 0., MASS_PROTON);
}

//-----------------------------------------------------------------------------

/* Collect the simulated particle properties from the file */
void TA2SaschaPhysics::SetTrueParticles()
{
	Float_t* dircos = (Float_t*)(fEvent[EI_dircos]);
	Float_t* elab = (Float_t*)(fEvent[EI_elab]);
	Float_t* plab = (Float_t*)(fEvent[EI_plab]);
	Int_t* idpart = (Int_t*)(fEvent[EI_idpart]);
	Int_t npart = *(Int_t*)(fEvent[EI_npart]);

	Double_t x, y, z;
	TLorentzVector p4;
	trueNPart = npart;

	for (Int_t i = 0; i < trueNPart; i++) {
		x = *dircos++;
		y = *dircos++;
		z = *dircos++;
		p4.SetT(elab[i]*1000.);
		p4.SetX(plab[i]*x*1000.);
		p4.SetY(plab[i]*y*1000.);
		p4.SetZ(plab[i]*z*1000.);
		trueIDPart[i] = idpart[i];

		switch (idpart[i]) {
		case 1:
			trueP4Gamma[trueNGamma] = p4;
			trueNGamma++;
			break;
		case 2:
			trueP4Posi[trueNPosi] = p4;
			trueNPosi++;
			break;
		case 3:
			trueP4Elec[trueNElec] = p4;
			trueNElec++;
			break;
		case 5:  // mu+
			trueP4MuPls[trueNMuPls] = p4;
			trueNMuPls++;
			break;
		case 6:  // mu-
			trueP4MuMns[trueNMuMns] = p4;
			trueNMuMns++;
			break;
		case 8:  // pi+
			trueP4PiPls[trueNPiPls] = p4;
			trueNPiPls++;
			break;
		case 9:  // pi-
			trueP4PiMns[trueNPiMns] = p4;
			trueNPiMns++;
			break;
		case 14:
			trueP4Prot[trueNProt] = p4;
			trueNProt++;
			break;
		default:
			break;
		}
	}
}

//-----------------------------------------------------------------------------

/* Prepare the gathered simulation information of the particles for being plotted */
void TA2SaschaPhysics::PlotTrueParticles()
{
	for (Int_t i = 0; i < trueNProt; i++) {
		trueProtMass[i] = trueP4Prot[i].M();
		trueProtEnergy[i] = trueP4Prot[i].E() - trueP4Prot[i].M();
		trueProtTheta[i] = trueP4Prot[i].Theta()*R2D;
		trueProtPhi[i] = trueP4Prot[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNGamma; i++) {
		trueGammaMass[i] = trueP4Gamma[i].M();
		trueGammaEnergy[i] = trueP4Gamma[i].E() - trueP4Gamma[i].M();
		trueGammaTheta[i] = trueP4Gamma[i].Theta()*R2D;
		trueGammaPhi[i] = trueP4Gamma[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNPosi; i++) {
		truePosiMass[i] = trueP4Posi[i].M();
		truePosiEnergy[i] = trueP4Posi[i].E() - trueP4Posi[i].M();
		truePosiTheta[i] = trueP4Posi[i].Theta()*R2D;
		truePosiPhi[i] = trueP4Posi[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNElec; i++) {
		trueElecMass[i] = trueP4Elec[i].M();
		trueElecEnergy[i] = trueP4Elec[i].E() - trueP4Elec[i].M();
		trueElecTheta[i] = trueP4Elec[i].Theta()*R2D;
		trueElecPhi[i] = trueP4Elec[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNMuPls; i++) {
		trueMuPlsMass[i] = trueP4MuPls[i].M();
		trueMuPlsEnergy[i] = trueP4MuPls[i].E() - trueP4MuPls[i].M();
		trueMuPlsTheta[i] = trueP4MuPls[i].Theta()*R2D;
		trueMuPlsPhi[i] = trueP4MuPls[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNMuMns; i++) {
		trueMuMnsMass[i] = trueP4MuMns[i].M();
		trueMuMnsEnergy[i] = trueP4MuMns[i].E() - trueP4MuMns[i].M();
		trueMuMnsTheta[i] = trueP4MuMns[i].Theta()*R2D;
		trueMuMnsPhi[i] = trueP4MuMns[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNPiPls; i++) {
		truePiPlsMass[i] = trueP4PiPls[i].M();
		truePiPlsEnergy[i] = trueP4PiPls[i].E() - trueP4PiPls[i].M();
		truePiPlsTheta[i] = trueP4PiPls[i].Theta()*R2D;
		truePiPlsPhi[i] = trueP4PiPls[i].Phi()*R2D;
	}

	for (Int_t i = 0; i < trueNPiMns; i++) {
		truePiMnsMass[i] = trueP4PiMns[i].M();
		truePiMnsEnergy[i] = trueP4PiMns[i].E() - trueP4PiMns[i].M();
		truePiMnsTheta[i] = trueP4PiMns[i].Theta()*R2D;
		truePiMnsPhi[i] = trueP4PiMns[i].Phi()*R2D;
	}
}

//-----------------------------------------------------------------------------

// Tank you for reading!
