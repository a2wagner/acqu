/**
 * TA2KinFitPhysics
 *
 * Physics class using the class TA2KinFit which serves as a wrapper class for KinFitter
 */

#include "TA2KinFitPhysics.h"

ClassImp(TA2KinFitPhysics)

//-----------------------------------------------------------------------------

TA2KinFitPhysics::TA2KinFitPhysics(const char* Name, TA2Analysis* Analysis) : TA2BasePhysics(Name, Analysis)
{
init = false;
}

//-----------------------------------------------------------------------------

TA2KinFitPhysics::~TA2KinFitPhysics()
{
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

	delete[] particles;
}

//---------------------------------------------------------------------------

void TA2KinFitPhysics::SetConfig(Char_t* line, Int_t key)
{
	// default SetConfig()
	TA2BasePhysics::SetConfig(line, key);
}

//---------------------------------------------------------------------------

void TA2KinFitPhysics::LoadVariable()
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

	TA2DataManager::LoadVariable("invM_2neutral", &invM_2neutral, EDSingleX);
	TA2DataManager::LoadVariable("protEnergyReconstr", &protEnergyReconstr, EDSingleX);

	TA2DataManager::LoadVariable("invMass2g", &invMass2g, EDSingleX);
	TA2DataManager::LoadVariable("invMass2g1p", &invMass2g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g", &invMass6g, EDSingleX);
	TA2DataManager::LoadVariable("invMass6g1p", &invMass6g1p, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB", &invMass2CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass2CB1TAPS", &invMass2CB1TAPS, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB", &invMass6CB, EDSingleX);
	TA2DataManager::LoadVariable("invMass6CB1TAPS", &invMass6CB1TAPS, EDSingleX);
}

//---------------------------------------------------------------------------

void TA2KinFitPhysics::PostInit()
{init = true;
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

	KinFit = TA2KinFit();

	// Call default PostInit()
	// Important!: Call this AFTER the previous arrays are defined, otherwise LoadVariable gets called and all arrays are still NULL
	TA2BasePhysics::PostInit();
}

//-----------------------------------------------------------------------------

/* Main method for applying the physics related stuff */
void TA2KinFitPhysics::Reconstruct()
{if (!init) PrintError("No PostInit() executed!", "Reconstruct()", EErrFatal);
	//Perform basic physics tasks
	TA2BasePhysics::Reconstruct();

	//Initialise array counters
	VarInit();

	Char_t currentFile[1024];
	gUAN->ReadRunName(currentFile);  // save filename without path and ending in the specified char array

	if (dbg)
		std::cout << "Processing file " << currentFile 
			<< ", event number: " << gAN->GetNEvent() 
			<< " (" << gAN->GetNEvAnalysed() << " accepted)" << std::endl;

	// Print speed information if set in config
	if (speed_info) {
		static int START = time(NULL);
		static int start = START;
		if (gAN->GetNEvent()%500000 == 0)
			printf("Processed %.1fM events after %d:%.02d min\n", gAN->GetNEvent()/1000000., (int)(time(NULL) - START)/60, (int)(time(NULL) - START)%60);
		if (dbg) {
		static int nEvents = 0;
			if ((time(NULL) - start) >= 60) {
				start = time(NULL);
				printf("Processed %d events after %d minutes (%d events/s)\n", gAN->GetNEvent(), (int)(time(NULL) - START)/60, nEvents/60);
				nEvents = gAN->GetNEvent();
			}
		}
	}

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



std::cout << "KinFitter status: " << KinFit.getStatus() << std::endl;
KinFit.print();

	TLorentzVector tmpState(0., 0., 0., 0.);
	int chargedPart[N_FINAL_STATE];
	int neutralPart[N_FINAL_STATE];
	TLorentzVector missProton;
	// variables to store values only for cuts
	Double_t mMiss, protEexpect;

	if (nParticles == N_FINAL_STATE) {  // proceed if there are 3 particles in total (2g and proton)

		nCharged = nNeutral = 0;
		for (unsigned int i = 0; i < N_FINAL_STATE; i++) {
			if (particles[i].HasDetector(EDetPID) || particles[i].HasDetector(EDetVeto))
				chargedPart[nCharged++] = i;
			else {
				neutralPart[nNeutral++] = i;
				tmpState += particles[i].GetP4();
			}
		}
		/* end event selection when we don't have one charged particle (entry in PID or Veto) and two neutral ones */
		if (nCharged != 1 && nNeutral != 2)
			return;

		protEnergyReconstr = particles[nParticlesCB].GetP4().E() - particles[nParticlesCB].GetP4().M();
		invM_2neutral = tmpState.M();

		/* At this point we have one charged and two neutral particles. Start to examine all tagged photons for prompt and random windows. */
		for (unsigned int i = 0; i < nBeamPhotons; i++) {
			missProton = fP4target[0] + Tagged[i].GetP4() - tmpState;
			mMiss = missProton.M();
			protEexpect = missProton.E() - missProton.M();
		}


		/* prepare kinematic fit */
		int err;
		int ndf;
		double chisq;
		double prob;

		TVector3 photon1 = particles[neutralPart[0]].GetVect();
		TVector3 photon2 = particles[neutralPart[1]].GetVect();
		TVector3 proton = particles[chargedPart[0]].GetVect();

		TMatrixD covPhoton1;
		TMatrixD covPhoton2;
		TMatrixD covProton;

		TLorentzVector fitPhoton1;
		TLorentzVector fitPhoton2;
		TLorentzVector fitProton;

		TFitParticlePThetaPhi ph1("neutral1", "neutral1", &photon1, 0., &covPhoton1);
		TFitParticlePThetaPhi ph2("neutral2", "neutral2", &photon2, 0., &covPhoton2);
		TFitParticlePThetaPhi pr("charged1", "charged1", &proton, 0., &covProton);

		int rows = 3;  // number of rows equal number of cols
		Double_t errors[3];
		int currPart = neutralPart[0];
		errors[0] = particles[currPart].GetSigmaE();
		errors[0] *= errors[0];
		errors[1] = particles[currPart].GetSigmaPhi();
		errors[1] *= errors[1];
		errors[2] = particles[currPart].GetSigmaTheta();
		errors[2] *= errors[2];
		if (!KinFit.fillSquareMatrixDiagonal(&covPhoton1, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		currPart = neutralPart[1];
		errors[0] = particles[currPart].GetSigmaE();
		errors[0] *= errors[0];
		errors[1] = particles[currPart].GetSigmaPhi();
		errors[1] *= errors[1];
		errors[2] = particles[currPart].GetSigmaTheta();
		errors[2] *= errors[2];
		if (!KinFit.fillSquareMatrixDiagonal(&covPhoton2, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		currPart = chargedPart[0];
		errors[0] = particles[currPart].GetSigmaE();
		errors[0] *= errors[0];
		errors[1] = particles[currPart].GetSigmaPhi();
		errors[1] *= errors[1];
		errors[2] = particles[currPart].GetSigmaTheta();
		errors[2] *= errors[2];
		if (!KinFit.fillSquareMatrixDiagonal(&covProton, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");

		TVector3 beam = Tagged[0].GetVect();
		TVector3 target = fP4target[0].Vect();
		TMatrixD covBeam;
		TMatrixD covTarget;
		TFitParticlePThetaPhi bm("beam", "beam", &beam, 0., &covBeam);
		TFitParticlePThetaPhi trgt("target", "target", &target, 0., &covTarget);
		errors[0] = Tagged[0].GetSigmaE();
		errors[0] *= errors[0];
		errors[1] = Tagged[0].GetSigmaPhi();
		errors[1] *= errors[1];
		errors[2] = Tagged[0].GetSigmaTheta();
		errors[2] *= errors[2];
		if (!KinFit.fillSquareMatrixDiagonal(&covBeam, errors, rows))
			fprintf(stderr, "Error filling covariance matrix with uncertainties\n");
		covTarget.Zero();
		covTarget.ResizeTo(3, 3);

		// energy and momentum constraints have to be defined separately for each component. components can be accessed via enum TFitConstraintEp::component
		TFitConstraintEp energyConservation("energyConstr", "Energy conservation constraint", 0, TFitConstraintEp::E, 0.);
		energyConservation.addParticles1(&bm, &trgt);
		energyConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pxConservation("pxConstr", "Px conservation constraint", 0, TFitConstraintEp::pX, 0.);
		pxConservation.addParticles1(&bm, &trgt);
		pxConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pyConservation("pyConstr", "Py conservation constraint", 0, TFitConstraintEp::pY, 0.);
		pyConservation.addParticles1(&bm, &trgt);
		pyConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintEp pzConservation("pzConstr", "Pz conservation constraint", 0, TFitConstraintEp::pZ, 0.);
		pzConservation.addParticles1(&bm, &trgt);
		pzConservation.addParticles2(&ph1, &ph2, &pr);
		TFitConstraintM massConstrProton("massConstr_proton", "mass constraint proton", 0, 0, MASS_PROTON);
		massConstrProton.addParticle1(&pr);

		//TKinFitter fit;
		KinFit.addMeasParticle(&ph1);
		KinFit.addMeasParticle(&ph1);
		KinFit.addMeasParticle(&pr);
		KinFit.setParamUnmeas(&pr, 0);  // proton energy unmeasured
		KinFit.addMeasParticle(&bm);
		KinFit.addMeasParticle(&trgt);
		KinFit.setParamUnmeas(&trgt, 0);
		KinFit.setParamUnmeas(&trgt, 1);
		KinFit.setParamUnmeas(&trgt, 2);

		KinFit.addConstraint(&massConstrProton);
		KinFit.addConstraint(&energyConservation);
		KinFit.addConstraint(&pxConservation);
		KinFit.addConstraint(&pyConservation);
		KinFit.addConstraint(&pzConservation);

		KinFit.setMaxNbIter(50);  // number of maximal iterations
		KinFit.setMaxDeltaS(5e-5);  // max Delta chi2
		KinFit.setMaxF(1e-4);  // max sum of constraints
		KinFit.setVerbosity(1);  // verbosity level
		KinFit.fit();

		fitPhoton1 = (*ph1.getCurr4Vec());
		fitPhoton2 = (*ph2.getCurr4Vec());
		fitProton = (*pr.getCurr4Vec());

		std::cout << "Fit result: " << std::endl;
		KinFit.print();
		ndf = KinFit.getNDF();
		chisq = KinFit.getS();
		prob = TMath::Prob(chisq, ndf);
		std::cout << "\nProbability: " << prob << "\tchi^2: " << chisq << std::endl;

	}





	// fast check for multiple gamma/CB events with/out proton/TAPS-Cluster
	int nPhoton = 0;
	bool hasProton = false;
	TA2Particle photons[maxParticles];
	//if (nParticles == 3)
	for (unsigned int i = 0; i < nParticles; i++) {
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
void TA2KinFitPhysics::GetParticles()
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
void TA2KinFitPhysics::GetTrueParticles()
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

//-----------------------------------------------------------------------------

void TA2KinFitPhysics::VarInit()
{
	//Initally, set single value histograms to EBufferEnd (this value will not be plotted)
	invM_2neutral = EBufferEnd;
	protEnergyReconstr = EBufferEnd;

	invMass2g = EBufferEnd;
	invMass2g1p = EBufferEnd;
	invMass6g = EBufferEnd;
	invMass6g1p = EBufferEnd;
	invMass2CB = EBufferEnd;
	invMass2CB1TAPS = EBufferEnd;
	invMass6CB = EBufferEnd;
	invMass6CB1TAPS = EBufferEnd;
}

//-----------------------------------------------------------------------------

/* Terminate all arrays with EBufferEnd markers at the end of the Reconstruct() method */
void TA2KinFitPhysics::TermArrays()
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
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void TA2KinFitPhysics::ParseMisc(char* line)
{
	TA2BasePhysics::ParseMisc(line);
}

//-----------------------------------------------------------------------------

/* Reconstruct the beam information which were used in the Pluto event generator */
void TA2KinFitPhysics::GetTrueBeam()
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
void TA2KinFitPhysics::SetTrueParticles()
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
void TA2KinFitPhysics::PlotTrueParticles()
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

