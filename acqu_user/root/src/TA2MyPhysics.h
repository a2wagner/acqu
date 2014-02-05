// SVN Info: $Id: TA2MyPhysics.h 780 2011-02-03 13:38:26Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2008-2011
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TA2MyPhysics                                                         //
//                                                                      //
// This is an abstract base physics analysis class.                     //
// Classes inheriting from this class must implement the LoadVariable() //
// and the ReconstructPhysics() method.                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef __TA2MyPhysics_h__
#define __TA2MyPhysics_h__

#include "TOSUtils.h"
#include "TOLoader.h"
#include "TOGlobals.h"
#include "TOMCParticle.h"
#include "TOA2RecParticle.h"

#include "TNtupleD.h"
#include "TChain.h"
#include "TF1.h"
#include "THnSparse.h"

#include "TA2Physics.h"
#include "TA2Tagger.h"
#include "TA2MyCrystalBall.h"
#include "TA2MyTAPS.h"
#include "TA2LongScint.h"
#include "TA2GenericApparatus.h"
#include "TA2GenericDetector.h"
#include "CaLibReader_t.h"


enum {
    EMP_USE_CALIB = 20000,
    EMP_CALIB_MISC,
    EMP_CALIB_TAGG,
    EMP_CALIB_CB,
    EMP_CALIB_TAPS,
    EMP_CALIB_PID,
    EMP_CALIB_VETO,
    EMP_BAD_TAGG_CH,
    EMP_P2_TAGGER_RATIO,
    EMP_TRIGGER_TOTAL_MULT,
    EMP_TRIGGER_TAPS_MULT,
    EMP_TRIGGER_CB_THRESHOLDS,
    EMP_TRIGGER_TAPS_LED1_THRESHOLD,
    EMP_TRIGGER_TAPS_LED2_THRESHOLD,
    EMP_CB_TIME_OFFSET,
    EMP_TAPS_TIME_OFFSET,
    EMP_CB_PID_PHI_LIMIT,
    EMP_BEAM_POL_BITS,
    EMP_CALIB_CBENERGY_PER_RUN,
};


static const Int_t gMyPhysics_MaxTaggerCh = 352;
static const Int_t gMyPhysics_MaxScalers  = 200;
static const Int_t gMyPhysics_MaxNaI      = 720;
static const Int_t gMyPhysics_MaxPID      =  24;
static const Int_t gMyPhysics_MaxBaF2PWO  = 438;
static const Int_t gMyPhysics_MaxTAPSRing =  11;
static const Int_t gMyPhysics_MaxMCPart   =  30;


static const Map_t myPhysicsConfigKeys[] = {
    // General keys
    {"Use-CaLib:"                    , EMP_USE_CALIB},                      // key for CaLib activation
    {"Use-CaLib-Misc:"               , EMP_CALIB_MISC},                     // key for miscellaneous CaLib configuration
    {"Use-CaLib-TAGG:"               , EMP_CALIB_TAGG},                     // key for CaLib tagger configuration
    {"Use-CaLib-CB:"                 , EMP_CALIB_CB},                       // key for CaLib CB configuration
    {"Use-CaLib-TAPS:"               , EMP_CALIB_TAPS},                     // key for CaLib TAPS configuration
    {"Use-CaLib-PID:"                , EMP_CALIB_PID},                      // key for CaLib PID configuration
    {"Use-CaLib-Veto:"               , EMP_CALIB_VETO},                     // key for CaLib Veto configuration
    {"Bad-Tagger-Channels:"          , EMP_BAD_TAGG_CH},                    // key for bad tagger channel read-in
    {"P2-Tagger-Ratio:"              , EMP_P2_TAGGER_RATIO},                // key for the p2-tagger-ratio output
    {"Trigger-Total-Mult:"           , EMP_TRIGGER_TOTAL_MULT},             // key for the trigger total multiplicity
    {"Trigger-TAPS-Mult:"            , EMP_TRIGGER_TAPS_MULT},              // key for the trigger TAPS multiplicity
    {"Trigger-CB-Thresholds:"        , EMP_TRIGGER_CB_THRESHOLDS},          // key for the trigger CB sum and segment thresholds
    {"Trigger-TAPS-LED1-Threshold:"  , EMP_TRIGGER_TAPS_LED1_THRESHOLD},    // key for the trigger TAPS LED1 threshold
    {"Trigger-TAPS-LED2-Threshold:"  , EMP_TRIGGER_TAPS_LED2_THRESHOLD},    // key for the trigger TAPS LED2 threshold
    {"CB-Time-Offset:"               , EMP_CB_TIME_OFFSET},                 // key for the CB time offset
    {"TAPS-Time-Offset:"             , EMP_TAPS_TIME_OFFSET},               // key for the TAPS time offset
    {"CB-PID-Phi-Limit:"             , EMP_CB_PID_PHI_LIMIT},               // key for the CB-PID coincidence phi limit
    {"Beam-Pol-Bits:"                , EMP_BEAM_POL_BITS},                  // key for the beam helicity bit definition
    {"Use-CaLib-CBEnergyPerRun:"     , EMP_CALIB_CBENERGY_PER_RUN},         // key for CaLib CBEnergy per Run configuration
    // Termination
    {NULL        , -1           }
};


class TA2MyPhysics : public TA2Physics
{

private:
    // ----------------------------------- General ----------------------------------- 
    TH1* fH_EventInfo;                                      // general event information
    TH1* fH_Corrected_Scalers;                              // overflow corrected current scalers
    TH1* fH_Corrected_SumScalers;                           // overflow corrected accumulated scalers
    TH1* fH_MCVertX;                                        // MC vertex x-coordinate
    TH1* fH_MCVertY;                                        // MC vertex y-coordinate
    TH1* fH_MCVertZ;                                        // MC vertex z-coordinate
    Double_t* fOldScalerSumValue;                           // used to avoid double sum scaler addition
    Double_t fOldP2ScalerSum;                               // used to avoid double sum scaler addition
    Int_t fNBadTaggerChannels;                              // number of bad tagger channels
    Int_t* fBadTaggerChannels;                              // list of bad tagger channels
    TNtupleD** fNt_P2TaggerRatio;                           // ntuple for the ratio p2/tagger
    Double_t fTargetPosition;                               // target position in z-direction    

    // ------------------------------- Time correction ------------------------------- 
    Double_t fCBTimeOffset;                                 // global CB time offset
    Double_t fTAPSTimeOffset;                               // global TAPS time offset
    
    // ------------------------------------- PID ------------------------------------- 
    TH1** fH_CB_PID_Coinc_Hits;                             // CB-PID coincidence hits
    TH2** fH_CB_PID_dE_E;                                   // CB-PID dE vs. E plots
    Double_t fCBPIDPhiLimit;                                // CB-PID coincidence phi limit

    // ------------------------------------ CaLib ------------------------------------ 
    CaLibReader_t* fCaLibReader;                            // CaLib reader
    
    // --------------------------- CBEnergy Correction per Run -----------------------
    Bool_t			CBEnergyPerRunCorrection;
	Double_t		CBEnergyPerRunCorrectionFactor;

    // ------------------------ Utility methods for this class ----------------------- 
    Int_t ExtrtactRunNumber();
    void LoadDetectors(TA2DataManager* parent, Int_t depth); 
    void UpdateCorrectedScaler(Int_t sc);
    Int_t GetPatternIndex(const Char_t* patternName, TA2BitPattern* pattern);
    void ApplyCaLib();

protected:
    // ----------------------------------- General ----------------------------------- 
    Int_t fAnalysisMode;                                    // analysis mode (raw, MC, ...)
    Bool_t fIsMC;                                           // flag for MC analysis mode
    Double_t fMCVertX;                                      // MC vertex x-coordinate
    Double_t fMCVertY;                                      // MC vertex y-coordinate
    Double_t fMCVertZ;                                      // MC vertex z-coordinate
    Int_t fRunNumber;                                       // run number
    Long64_t fEventCounter;                                 // event counter
    Long64_t fEventOffset;                                  // event offset when analyzing multiple files
    Int_t fSaveEvent;                                       // if 1 : save current event in reduced AcquRoot ROOT file

    // ----------------------------------- Tagger ----------------------------------- 
    TA2Tagger* fTagger;                                     // pointer to the Tagger
    TA2Ladder* fLadder;                                     // pointer to the Ladder
    UInt_t fTaggerPhotonNhits;                              // number of photons in the tagger
    Int_t* fTaggerPhotonHits;                               // pointer to the Tagger photon hits
    Double_t* fTaggerPhotonEnergy;                          // pointer to the Tagger photon energy array
    Double_t* fTaggerPhotonTime;                            // pointer to the Tagger photon time array

    // ------------------------------------- CB ------------------------------------- 
    TA2MyCrystalBall* fCB;                                  // pointer to the Crystal Ball
    TA2ClusterDetector* fNaI;                               // pointer to the NaI elements
    UInt_t fNaINhits;                                       // number of NaI hits
    Int_t* fNaIHits;                                        // pointer to the NaI hits
    Double_t* fNaIEnergy;                                   // pointer to the NaI energy array
    Double_t* fNaITime;                                     // pointer to the NaI time array
    Double_t* fCBClTime;                                    // pointer to the CB cluster time
    Int_t* fCBPID_Index;                                    // pointer to the PID index of the CB cluster
    Double_t* fCBPID_dE;                                    // pointer to the PID energy deposition of the CB cluster
    UInt_t* fCBClusterHit;                                  // pointer to the CB cluster hit indices
    UInt_t fCBNCluster;                                     // number of clusters in CB
    HitCluster_t** fCBCluster;                              // pointer to the list of CB clusters
    
    // ------------------------------------- PID ------------------------------------- 
    TA2Detector* fPID;                                      // pointer to the PID
    UInt_t fPIDNhits;                                       // number of hits in the PID
    Int_t* fPIDHits;                                        // pointer to the PID hits
    Double_t* fPIDEnergy;                                   // pointer to the PID energy array
    Double_t* fPIDTime;                                     // pointer to the PID time array
    TVector3** fPIDHitPos;                                  // pointer to the PID hit position list

    // ------------------------------------ TAPS ------------------------------------ 
    TA2MyTAPS* fTAPS;                                       // pointer to TAPS
    TA2MyTAPS_BaF2PWO* fBaF2PWO;                            // pointer to the BaF2 (or the BaF2/PWO) array
    UInt_t fBaF2PWONhits;                                   // number of BaF2/PWO hits
    Int_t* fBaF2PWOHits;                                    // pointer to the BaF2/PWO hits
    Double_t* fBaF2PWOEnergy;                               // pointer the the BaF2/PWO energy array
    Double_t* fBaF2PWOTime;                                 // pointer to the BaF2/PWO time
    
    TA2BitPattern* fBaF2Pattern;                            // pointer to the TAPS bit pattern
    Int_t fBaF2CFDPattern;                                  // index of the CFD pattern in the TAPS bit pattern
    Int_t fBaF2LED1Pattern;                                 // index of the LED1 pattern in the TAPS bit pattern
    Int_t fBaF2LED2Pattern;                                 // index of the LED2 pattern in the TAPS bit pattern
    UInt_t fBaF2CFDNhits;                                   // number of CFD hits in the BaF2
    Int_t* fBaF2CFDHits;                                    // pointer to the BaF2 CFD hits
    UInt_t fBaF2LED1Nhits;                                  // number of LED1 hits in the BaF2
    Int_t* fBaF2LED1Hits;                                   // pointer to the BaF2 LED1 hits
    UInt_t fBaF2LED2Nhits;                                  // number of LED2 hits in the BaF2
    Int_t* fBaF2LED2Hits;                                   // pointer to the BaF2 LED2 hits
    
    UInt_t* fTAPSClusterHit;                                // pointer to the TAPS cluster hit indices
    UInt_t fTAPSNCluster;                                   // number of clusters in TAPS
    HitCluster_t** fTAPSCluster;                            // pointer to the list of TAPS clusters
    
    // ------------------------------------ Veto ------------------------------------ 
    TA2Detector* fVeto;                                     // pointer to the TAPS Vetos
    UInt_t fVetoNhits;                                      // number of Veto hits
    Int_t* fVetoHits;                                       // pointer to the Veto hits
    Double_t* fVetoEnergy;                                  // pointer to the Veto energy array
    Double_t* fVetoTime;                                    // pointer to the Veto time array
    
    // ------------------------------------ PbWO4 ----------------------------------- 
    TA2GenericDetector* fPbWO4;                             // pointer to the PbWO4 ring of TAPS
    UInt_t fPbWO4Nhits;                                     // number of PbWO4 hits
    Int_t* fPbWO4Hits;                                      // pointer to the PbWO4 hits array
    Double_t* fPbWO4Energy;                                 // pointer to the PbWO4 energy array
    Double_t* fPbWO4Time;                                   // pointer to the PbWO4 time array

    // -------------------------------- Polarisation -------------------------------- 
    UShort_t* fCP;                                          // Circular polarisation: pointer to the electron beam helicity state
    Int_t fEBeamBitPos;                                     // positive helicity bit value
    Int_t fEBeamBitNeg;                                     // negative helicity bit value
    Int_t fEBeamHelState;                                   // electron beam helicity state value (+1/-1)
    
    // ------------------------------------ TOF ------------------------------------- 
    TA2LongScint* fTOF;                                     // pointer to the TOF wall
    UInt_t fTOFBarNhits;                                    // number of TOF bar hits
    Int_t* fTOFBarHits;                                     // pointer to the TOF bar hits
    Double_t* fTOFBarMeanEnergy;                            // pointer to the TOF bar mean energy
    Double_t* fTOFBarMeanTime;                              // pointer to the TOF bar mean time
  
    // ----------------------------- Pb glass detector ------------------------------ 
    TA2GenericApparatus* fPbGlassApp;                       // pointer to the Pb glass detector apparatus
    TA2GenericDetector* fPbGlass;                           // pointer to the Pb glass detector
    UInt_t fPbGlassNhits;                                   // number of hits in the Pb glass detector
    Double_t* fPbGlassTime;                                 // pointer to the Pb glass detector time
    
    // ----------------------------- Particle collection ---------------------------- 
    Int_t fNNeutral;                                        // number of detected neutral particles
    Int_t fNCharged;                                        // number of detected charged particles
    TOA2DetParticle** fPartCB;                              // array of particles detected in CB
    TOA2DetParticle** fPartTAPS;                            // array of particles detected in TAPS
    TOA2DetParticle** fPartNeutral;                         // array of detected neutral particles
    TOA2DetParticle** fPartCharged;                         // array of detected charged particles
    
    // ---------------------------- MC Particle collection --------------------------- 
    Int_t fNMC;                                             // number of MC particles
    TOMCParticle* fPartMC[gMyPhysics_MaxMCPart];            // array of MC particles
    
    // ------------------------------- Hardware trigger  ----------------------------- 
    Bool_t fL1CB;                                         // Level 1 CB sum trigger
    Bool_t fL1TAPS_OR;                                    // Level 1 TAPS OR trigger
    Bool_t fL1TAPS_M2;                                    // Level 1 TAPS M2 trigger
    Bool_t fL1Pulser;                                     // Level 1 Pulser trigger
    Bool_t fL1TAPS_Pulser;                                // Level 1 TAPS Pulser trigger
    Bool_t fL2M1;                                         // Level 2 M1+ trigger
    Bool_t fL2M2;                                         // Level 2 M2+ trigger
    Bool_t fL2M3;                                         // Level 2 M3+ trigger
    Bool_t fL2M4;                                         // Level 2 M4+ trigger
    
    // ------------------------------- Software trigger  ----------------------------- 
    Int_t fTrigTotalMult;                                   // total multiplicity
    Int_t fTrigTAPSMult;                                    // TAPS alone multiplicity
    Double_t fTrigCBSum;                                    // CB energy sum
    Double_t fTrigCBSegThr;                                 // CB NaI segment high discriminator threshold
    Double_t fTrigTAPSLED1Thr;                              // TAPS LED1 thresholds
    Double_t fTrigTAPSLED2Thr;                              // TAPS LED2 thresholds

    // -------------------- Utility methods for analysis classes -------------------- 
    Bool_t IsBadTaggerChannel(Int_t ch);
    Bool_t IsTrigger(TOA2RecParticle& mesons);
    void ClearCBTimeWalk();
    Double_t GetBGSubtractionWeight(TOA2RecParticle& meson,
                                    Double_t taggerTime,
                                    Double_t limitPromptLowCB, Double_t limitPromptHighCB,
                                    Double_t limitBG1LowCB, Double_t limitBG1HighCB,
                                    Double_t limitBG2LowCB, Double_t limitBG2HighCB,
                                    Double_t limitPromptLowTAPS, Double_t limitPromptHighTAPS,
                                    Double_t limitBG1LowTAPS, Double_t limitBG1HighTAPS,
                                    Double_t limitBG2LowTAPS, Double_t limitBG2HighTAPS);
    Int_t GetVetoInFrontOfElement(Int_t id) const;
    Int_t GetTAPSRing(Int_t id) const;
    Int_t GetTAPSBlock(Int_t id) const;
    Double_t CheckClusterVeto(HitCluster_t* inCluster, Int_t* outVetoIndex = 0) const;
    Double_t CheckClusterPID(HitCluster_t* inCluster, Int_t* outPIDIndex = 0) const;
    Bool_t CheckPSA(TOA2DetParticle* p, TCutG* c);
    Int_t GetForeignVetoHits(TOA2DetParticle* p);

    void SetP4(TLorentzVector& p4, Double_t Ekin, Double_t mass, TVector3* dir);
    void SetP4(TLorentzVector& p4, Double_t Ekin, Double_t mass, Double_t x, Double_t y, Double_t z);
    
    Bool_t GetPSA(TOA2DetParticle* p, Double_t* psaR, Double_t* psaA);
    
    void FillPSA(TH2* h, TOA2DetParticle* p, Double_t w = 1.);
    void FilldE_E(TH2* hCB, TH2* hTAPS, TOA2DetParticle* p, Double_t w = 1.);
    void FillTOF(TH2* hCB, TH2* hTAPS, TOA2DetParticle* p, Double_t t, Bool_t isMC, Double_t w = 1.);
    void FillTaggerCoincidence(TOA2RecParticle& meson, Double_t taggerTime, TH1* hCB, TH1* hTAPS);
    void FillTaggerCoincidence2D(TOA2RecParticle& meson, Double_t taggerTime, Int_t taggerElement, TH2* hCB, TH2* hTAPS);

public:
    TA2MyPhysics(const char* name, TA2Analysis* analysis);
    virtual ~TA2MyPhysics();

    virtual void SetConfig(Char_t* line, Int_t key);
    virtual void PostInit();
    virtual void LoadVariable();
    virtual void ReconstructPhysics() = 0;
    void Reconstruct();
    virtual TA2DataManager* CreateChild( const char*, Int_t ) { return NULL; } 
    
    Int_t GetRunNumber() const { return fRunNumber; }
    
    ClassDef(TA2MyPhysics, 1)
};

#endif

