// SVN Info: $Id: CaLibReader_t.h 767 2011-01-26 16:54:13Z werthm $

/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller, 2010
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CaLibReader_t                                                        //
//                                                                      //
// CaLib database connection and reader class.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef __CaLibReader_t_h__
#define __CaLibReader_t_h__

#include "TError.h"
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TF1.h"

#include "TA2MyTAPS_BaF2PWO.h"
#include "TA2TAPS_BaF2.h"


class CaLibReader_t 
{

private:
    TSQLServer* fDB;                    // SQL server
    Int_t fRun;                         // run number
    Char_t fCalibration[256];           // calibration identifier
    Bool_t fTargetPosition;             // target position flag
    Bool_t fTAGGtime;                   // tagger time flag
    Bool_t fCBenergy;                   // CB energy flag
    Bool_t fCBquadEnergy;               // CB quadratic energy flag
    Bool_t fCBtime;                     // CB time flag
    Bool_t fCBwalk;                     // CB walk flag
    Bool_t fTAPSenergy;                 // TAPS energy flag
    Bool_t fTAPSquadEnergy;             // TAPS quadratic energy flag
    Bool_t fTAPStime;                   // TAPS time flag
    Bool_t fTAPSled;                    // TAPS LED flag
    Bool_t fPIDphi;                     // PID phi flag
    Bool_t fPIDdroop;                   // PID droop flag
    Bool_t fPIDenergy;                  // PID energy flag
    Bool_t fPIDtime;                    // PID time flag
    Bool_t fVetoenergy;                 // Veto energy flag
    Bool_t fVetotime;                   // Veto time flag
    
    Double_t* fCBQuadEnergyPar0;        // CB quadratic energy correction parameter 0
    Double_t* fCBQuadEnergyPar1;        // CB quadratic energy correction parameter 1
    Double_t* fTAPSQuadEnergyPar0;      // TAPS quadratic energy correction parameter 0
    Double_t* fTAPSQuadEnergyPar1;      // TAPS quadratic energy correction parameter 1
    Double_t* fTAPSLED1Thr;             // TAPS LED1 thresholds
    Double_t* fTAPSLED2Thr;             // TAPS LED2 thresholds
    Int_t fNPIDDroopFunc;               // number of PID droop correction functions
    TF1** fPIDDroopFunc;                // array of PID droop correction functions

    Bool_t ReadParameters(const Char_t* table, Double_t* par, Int_t length);

public:
    CaLibReader_t(const Char_t* dbHost, const Char_t* dbName, 
                  const Char_t* dbUser, const Char_t* dbPass,
                  const Char_t* calibration, Int_t run);
    virtual ~CaLibReader_t();
 
    Bool_t ApplyTargetPositioncalib(Double_t* pos);
    Bool_t ApplyTAGGcalib(TA2Detector* det);
    Bool_t ApplyCBcalib(TA2Detector* det);
    Bool_t ApplyTAPScalib(TA2MyTAPS_BaF2PWO* det);
    Bool_t ApplyTAPScalib(TA2TAPS_BaF2* det);
    Bool_t ApplyPIDcalib(TA2Detector* det);
    Bool_t ApplyVetocalib(TA2Detector* det);

    Bool_t IsConnected() { return fDB ? kTRUE : kFALSE; }
    void Deconnect()
    {
        delete fDB;
        fDB = 0;
    }

    void SetTargetPosition(Bool_t b = kTRUE) { fTargetPosition = b; }
    void SetTAGGtime(Bool_t b = kTRUE) { fTAGGtime = b; }
    void SetCBenergy(Bool_t b = kTRUE) { fCBenergy = b; }
    void SetCBquadEnergy(Bool_t b = kTRUE) { fCBquadEnergy = b; }
    void SetCBtime(Bool_t b = kTRUE) { fCBtime = b; } 
    void SetCBwalk(Bool_t b = kTRUE) { fCBwalk = b; }     
    void SetTAPSenergy(Bool_t b = kTRUE) { fTAPSenergy = b; }
    void SetTAPSquadEnergy(Bool_t b = kTRUE) { fTAPSquadEnergy = b; }
    void SetTAPStime(Bool_t b = kTRUE) { fTAPStime = b; }    
    void SetTAPSled(Bool_t b = kTRUE) { fTAPSled = b; }
    void SetPIDphi(Bool_t b = kTRUE) { fPIDphi = b; }
    void SetPIDdroop(Bool_t b = kTRUE) { fPIDdroop = b; }
    void SetPIDenergy(Bool_t b = kTRUE) { fPIDenergy = b; }
    void SetPIDtime(Bool_t b = kTRUE) { fPIDtime = b; }
    void SetVetoenergy(Bool_t b = kTRUE) { fVetoenergy = b; }
    void SetVetotime(Bool_t b = kTRUE) { fVetotime = b; }
    
    Bool_t GetCBquadEnergy() const { return fCBquadEnergy; }
    Bool_t GetTAPSquadEnergy() const { return fTAPSquadEnergy; }
    Bool_t GetPIDdroop() const { return fPIDdroop; }
    Bool_t GetTAPSled() const { return fTAPSled; }
    Double_t GetCBQuadEnergyPar0(Int_t n) const { return fCBQuadEnergyPar0[n]; }
    Double_t GetCBQuadEnergyPar1(Int_t n) const { return fCBQuadEnergyPar1[n]; }
    Double_t GetTAPSQuadEnergyPar0(Int_t n) const { return fTAPSQuadEnergyPar0[n]; }
    Double_t GetTAPSQuadEnergyPar1(Int_t n) const { return fTAPSQuadEnergyPar1[n]; }
    Double_t GetTAPSLED1Thr(Int_t n) const { return fTAPSLED1Thr[n]; }
    Double_t GetTAPSLED2Thr(Int_t n) const { return fTAPSLED2Thr[n]; }
    TF1* GetPIDDroopFunc(Int_t n) const { return fPIDDroopFunc[n]; }
};

#endif

