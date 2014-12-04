#ifndef __TA2KinFitPhysics_h__
#define __TA2KinFitPhysics_h__

#include "TA2BasePhysics.h"
#include "TA2KinFit.h"

class TA2Apparatus;

class TA2KinFitPhysics : public TA2BasePhysics
{
  protected:
	// some general constants
	static const int N_FINAL_STATE = 3;
	static const double MASS_ETAPRIME = 957.78;
	static const double MASS_PROTON = 938.272;
	static const double MASS_ETA = 547.862;
	static const double MASS_PIZERO = 134.9766;
	static const double R2D = 180./3.14159265359;

	// Apparati for the different used detectors
	TA2Tagger* fTAG;  // Tagger
	TA2CentralApparatus* fCB;  // Crystal Ball
	TA2Taps* fTAPS;  // TAPS
	TA2Ladder* fFPD;  // Ladder
	TA2Apparatus* fPID;
	TA2CalArray* fNaI;
	TA2TAPS_BaF2* fBaF2;
	TA2Apparatus* fVeto;
bool init;
	// Flag for checking what data type will be processed
	Bool_t MC;

	// Print speed information if specified in config file
	bool speed_info;

	// Print some output for debugging if set to true via config
	bool dbg;

	// Stuff used for a kinematic fit
	TA2KinFit KinFit;

	// ...
	TA2Particle* particles;
	UInt_t nParticles, maxParticles;
	UInt_t nParticlesCB, nParticlesTAPS, maxParticlesCB, maxParticlesTAPS;
	UInt_t nBeamPhotons, maxBeamPhotons;

	//just for testing
	Double_t invMass2g, invMass2g1p, invMass6g, invMass6g1p, invMass2CB, invMass2CB1TAPS, invMass6CB, invMass6CB1TAPS;

	// Variables to store the true information from the simulation
	TLorentzVector trueP4Beam;
	TLorentzVector trueP4Target;
	Int_t trueNPart;
	Int_t* trueIDPart;
	TLorentzVector* trueP4Gamma;
	TLorentzVector* trueP4Elec;
	TLorentzVector* trueP4Posi;
	TLorentzVector* trueP4Prot;
	TLorentzVector* trueP4MuPls;
	TLorentzVector* trueP4MuMns;
	TLorentzVector* trueP4PiPls;
	TLorentzVector* trueP4PiMns;
	Int_t trueNGamma;
	Int_t trueNPosi;
	Int_t trueNElec;
	Int_t trueNProt;
	Int_t trueNMuPls;
	Int_t trueNMuMns;
	Int_t trueNPiPls;
	Int_t trueNPiMns;
	Double_t* trueProtMass;
	Double_t* trueProtEnergy;
	Double_t* trueProtTheta;
	Double_t* trueProtPhi;
	Double_t* trueGammaMass;
	Double_t* trueGammaEnergy;
	Double_t* trueGammaTheta;
	Double_t* trueGammaPhi;
	Double_t* truePosiMass;
	Double_t* truePosiEnergy;
	Double_t* truePosiTheta;
	Double_t* truePosiPhi;
	Double_t* trueElecMass;
	Double_t* trueElecEnergy;
	Double_t* trueElecTheta;
	Double_t* trueElecPhi;
	Double_t* trueMuPlsMass;
	Double_t* trueMuPlsEnergy;
	Double_t* trueMuPlsTheta;
	Double_t* trueMuPlsPhi;
	Double_t* trueMuMnsMass;
	Double_t* trueMuMnsEnergy;
	Double_t* trueMuMnsTheta;
	Double_t* trueMuMnsPhi;
	Double_t* truePiPlsMass;
	Double_t* truePiPlsEnergy;
	Double_t* truePiPlsTheta;
	Double_t* truePiPlsPhi;
	Double_t* truePiMnsMass;
	Double_t* truePiMnsEnergy;
	Double_t* truePiMnsTheta;
	Double_t* truePiMnsPhi;

	UInt_t nCharged, nNeutral;
	Double_t invM_2neutral;
	Double_t protEnergyReconstr;


	//General functions
	void VarInit();  //Clear all used variables
	void TermArrays();  //Terminate arrays with EBufferEnd markers
	// Get (and plot) particle properties from simulation
	void GetTrueBeam();
	void SetTrueParticles();
	void PlotTrueParticles();
	// Get particles for analysis independent of previous identification
	void GetParticles();
	void GetTrueParticles();

  public:
	TA2KinFitPhysics(const char*, TA2Analysis*);
	virtual ~TA2KinFitPhysics();
	virtual void LoadVariable();  //Creates histograms
	virtual void SetConfig(Char_t*, Int_t);  //Parses general information from configuration file
	virtual void ParseMisc(char* line);  //Parses additional information from configuration file
	virtual void Reconstruct();  //Event reconstruction
	virtual void PostInit();  //Initialisation etc.

  ClassDef(TA2KinFitPhysics, 1)
};

//-----------------------------------------------------------------------------

#endif

