// SVN Info: $Id: CaLibReader_t.cc 767 2011-01-26 16:54:13Z werthm $

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


#include "CaLibReader_t.h"


//______________________________________________________________________________
CaLibReader_t::CaLibReader_t(const Char_t* dbHost, const Char_t* dbName, 
                             const Char_t* dbUser, const Char_t* dbPass,
                             const Char_t* calibration, Int_t run)
{
    // Constructor.
    
    // init members
    fDB = 0;
    fRun = run;
    strcpy(fCalibration, calibration);
    fTAGGtime = kFALSE;
    fCBenergy = kFALSE;    
    fCBquadEnergy = kFALSE;
    fCBtime = kFALSE;      
    fCBwalk = kFALSE;     
    fTAPSenergy = kFALSE;  
    fTAPSquadEnergy = kFALSE;
    fTAPStime = kFALSE;    
    fTAPSled = kFALSE;
    fPIDphi = kFALSE;   
    fPIDdroop = kFALSE;   
    fPIDenergy = kFALSE;   
    fPIDtime = kFALSE;     
    fVetoenergy = kFALSE;  
    fVetotime = kFALSE;    
    
    fCBQuadEnergyPar0 = 0;
    fCBQuadEnergyPar1 = 0;
    fTAPSQuadEnergyPar0 = 0;
    fTAPSQuadEnergyPar1 = 0;
    fTAPSLED1Thr = 0;
    fTAPSLED2Thr = 0;
    fNPIDDroopFunc = 0;
    fPIDDroopFunc = 0;

    // connect to the SQL server
    Char_t tmp[256];
    sprintf(tmp, "mysql://%s/%s", dbHost, dbName);
    fDB = TSQLServer::Connect(tmp, dbUser, dbPass);

    // check the connection to the server
    if (IsConnected()) Info("CaLibReader_t", "Successfully connected to CaLib database server");
    else Error("CaLibReader_t", "Could not connect to CaLib database server!");
}

//______________________________________________________________________________
CaLibReader_t::~CaLibReader_t() 
{
    // Destructor.
   
    if (fDB) delete fDB;
    if (fCBQuadEnergyPar0) delete [] fCBQuadEnergyPar0;
    if (fCBQuadEnergyPar1) delete [] fCBQuadEnergyPar1;
    if (fTAPSQuadEnergyPar0) delete [] fTAPSQuadEnergyPar0;
    if (fTAPSQuadEnergyPar1) delete [] fTAPSQuadEnergyPar1;
    if (fTAPSLED1Thr) delete [] fTAPSLED1Thr;
    if (fTAPSLED2Thr) delete [] fTAPSLED2Thr;
    if (fPIDDroopFunc) 
    {
        for (Int_t i = 0; i < fNPIDDroopFunc; i++) delete fPIDDroopFunc[i];
        delete [] fPIDDroopFunc;
    }
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTargetPositioncalib(Double_t* pos)
{
    // Read target positon calibration from database and write it to 'pos'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // target position calibration
    if (fTargetPosition)
    {
        Double_t targetPos;

        // read target position data
        if (!ReadParameters("target_pos", &targetPos, 1)) return kFALSE;
        
        // set target position
        *pos = targetPos;
      
        Info("ApplyTargetPositioncalib", "Target position calibration was applied: "
             "Target is at %f cm", targetPos);
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTAGGcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // get the number of detector elements
    Int_t nElem = det->GetNelement();
    
    // tagger time calibration
    if (fTAGGtime)
    {
        Double_t t0[nElem];
        
        // read TDC offset data
        if (!ReadParameters("tagg_t0", t0, nElem)) return kFALSE;
      
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{   
	    // set TDC offset
	    det->GetElement(i)->SetT0(t0[i]);
	}
      
        Info("ApplyTAGGcalib", "TAGG time calibration was applied");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyCBcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
  
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // CB energy calibration
    if (fCBenergy)
    {
        Double_t e1[nElem];
      
        // read SADC gain data
        if (!ReadParameters("cb_e1", e1, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set SADC gain
	    det->GetElement(i)->SetA1(e1[i]);
        }

        Info("ApplyCBcalib", "CB energy calibration was applied");
    }
  
    // CB time calibration
    if (fCBtime)
    {
        Double_t t0[nElem];
      
        // read TDC offset data
        if (!ReadParameters("cb_t0", t0, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set TDC offset
	    det->GetElement(i)->SetT0(t0[i]);
        }

        Info("ApplyCBcalib", "CB time calibration was applied");
    }

    // CB time walk calibration
    if (fCBwalk)
    {
        Double_t tw0[nElem], tw1[nElem], tw2[nElem], tw3[nElem];
        Double_t tw4[4];
        
        // read time walk data
        if (!ReadParameters("cb_walk0", tw0, nElem)) return kFALSE;
        if (!ReadParameters("cb_walk1", tw1, nElem)) return kFALSE;
        if (!ReadParameters("cb_walk2", tw2, nElem)) return kFALSE;
        if (!ReadParameters("cb_walk3", tw3, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{  
	    tw4[0]=tw0[i];
	    tw4[1]=tw1[i];
	    tw4[2]=tw2[i];
	    tw4[3]=tw3[i];
	
            // set time walk parameters
	    det->GetElement(i)->GetTimeWalk()->SetWalk(tw4);
        }
        
        Info("ApplyCBcalib", "CB time walk calibration was applied");
    }
    
    // CB quadratic energy correction
    if (fCBquadEnergy)
    {
        // create arrays
        if (!fCBQuadEnergyPar0) fCBQuadEnergyPar0 = new Double_t[nElem];
        if (!fCBQuadEnergyPar1) fCBQuadEnergyPar1 = new Double_t[nElem];

        // read quadratic energy correction data
        if (!ReadParameters("cb_equad0", fCBQuadEnergyPar0, nElem) || 
            !ReadParameters("cb_equad1", fCBQuadEnergyPar1, nElem)) 
        {
            fCBquadEnergy = kFALSE;
            return kFALSE;
        }

        Info("ApplyCBcalib", "CB quadratic energy correction was read");
    }
 
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTAPScalib(TA2MyTAPS_BaF2PWO* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // TAPS energy calibration
    if (fTAPSenergy)
    {
        Double_t e0_lg[nElem], e0_sg[nElem];
        Double_t e1_lg[nElem], e1_sg[nElem];
      
        // read ADC LG pedestal data
        if (!ReadParameters("taps_lg_e0", e0_lg, nElem)) return kFALSE;

        // read ADC SG pedestal data
        if (!ReadParameters("taps_sg_e0", e0_sg, nElem)) return kFALSE;
    
        // read ADC LG gain data
        if (!ReadParameters("taps_lg_e1", e1_lg, nElem)) return kFALSE;

        // read ADC SG gain data
        if (!ReadParameters("taps_sg_e1", e1_sg, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set LG ADC pedestal
	    det->GetLGElement(i)->SetA0(e0_lg[i]);
	  
	    // set SG ADC pedestal
	    det->GetSGElement(i)->SetA0(e0_sg[i]);
	  
	    // set LG ADC gain
	    det->GetLGElement(i)->SetA1(e1_lg[i]);
	  
	    // set SG ADC gain
	    det->GetSGElement(i)->SetA1(e1_sg[i]);
        }
        
        Info("ApplyTAPScalib", "TAPS energy calibration was applied");
    } 

    // TAPS time calibration
    if (fTAPStime)
    {  
        Double_t t0[nElem], t1[nElem];
        
        // read TDC offset data
        if (!ReadParameters("taps_t0", t0, nElem)) return kFALSE;
        
        // read TDC gain data
        if (!ReadParameters("taps_t1", t1, nElem)) return kFALSE;

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set LG TDC offset
	    det->GetLGElement(i)->SetT0(t0[i]);
	  
	    // set LG TDC gain
	    det->GetLGElement(i)->SetT1(t1[i]);
        }

        Info("ApplyTAPScalib", "TAPS time calibration was applied");
    }
    
    // TAPS quadratic energy correction
    if (fTAPSquadEnergy)
    {
        // create arrays
        if (!fTAPSQuadEnergyPar0) fTAPSQuadEnergyPar0 = new Double_t[nElem];
        if (!fTAPSQuadEnergyPar1) fTAPSQuadEnergyPar1 = new Double_t[nElem];

        // read quadratic energy correction data
        if (!ReadParameters("taps_equad0", fTAPSQuadEnergyPar0, nElem) ||
            !ReadParameters("taps_equad1", fTAPSQuadEnergyPar1, nElem))
        {
            fTAPSquadEnergy = kFALSE;
            return kFALSE;
        }

        Info("ApplyTAPScalib", "TAPS quadratic energy correction was read");
    }
    
    // TAPS LED thresholds
    if (fTAPSled)
    {
        // create arrays
        if (!fTAPSLED1Thr) fTAPSLED1Thr = new Double_t[nElem];
        if (!fTAPSLED2Thr) fTAPSLED2Thr = new Double_t[nElem];

        // read LED threshold data
        if (!ReadParameters("taps_led1", fTAPSLED1Thr, nElem) ||
            !ReadParameters("taps_led2", fTAPSLED2Thr, nElem))
        {
            fTAPSled = kFALSE;
            return kFALSE;
        }

        Info("ApplyTAPScalib", "TAPS LED thresholds were read");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTAPScalib(TA2TAPS_BaF2* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // TAPS energy calibration
    if (fTAPSenergy)
    {
        Double_t e0_lg[nElem], e0_sg[nElem];
        Double_t e1_lg[nElem], e1_sg[nElem];
      
        // read ADC LG pedestal data
        if (!ReadParameters("taps_lg_e0", e0_lg, nElem)) return kFALSE;

        // read ADC SG pedestal data
        if (!ReadParameters("taps_sg_e0", e0_sg, nElem)) return kFALSE;
    
        // read ADC LG gain data
        if (!ReadParameters("taps_lg_e1", e1_lg, nElem)) return kFALSE;

        // read ADC SG gain data
        if (!ReadParameters("taps_sg_e1", e1_sg, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set LG ADC pedestal
	    det->GetLGElement(i)->SetA0(e0_lg[i]);
	  
	    // set SG ADC pedestal
	    det->GetSGElement(i)->SetA0(e0_sg[i]);
	  
	    // set LG ADC gain
	    det->GetLGElement(i)->SetA1(e1_lg[i]);
	  
	    // set SG ADC gain
	    det->GetSGElement(i)->SetA1(e1_sg[i]);
        }
        
        Info("ApplyTAPScalib", "TAPS energy calibration was applied");
    } 

    // TAPS time calibration
    if (fTAPStime)
    {  
        Double_t t0[nElem], t1[nElem];
        
        // read TDC offset data
        if (!ReadParameters("taps_t0", t0, nElem)) return kFALSE;
        
        // read TDC gain data
        if (!ReadParameters("taps_t1", t1, nElem)) return kFALSE;

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set LG TDC offset
	    det->GetLGElement(i)->SetT0(t0[i]);
	  
	    // set LG TDC gain
	    det->GetLGElement(i)->SetT1(t1[i]);
        }

        Info("ApplyTAPScalib", "TAPS time calibration was applied");
    }
    
    // TAPS quadratic energy correction
    if (fTAPSquadEnergy)
    {
        // create arrays
        if (!fTAPSQuadEnergyPar0) fTAPSQuadEnergyPar0 = new Double_t[nElem];
        if (!fTAPSQuadEnergyPar1) fTAPSQuadEnergyPar1 = new Double_t[nElem];

        // read quadratic energy correction data
        if (!ReadParameters("taps_equad0", fTAPSQuadEnergyPar0, nElem) ||
            !ReadParameters("taps_equad1", fTAPSQuadEnergyPar1, nElem))
        {
            fTAPSquadEnergy = kFALSE;
            return kFALSE;
        }

        Info("ApplyTAPScalib", "TAPS quadratic energy correction was read");
    }
    
    // TAPS LED thresholds
    if (fTAPSled)
    {
        // create arrays
        if (!fTAPSLED1Thr) fTAPSLED1Thr = new Double_t[nElem];
        if (!fTAPSLED2Thr) fTAPSLED2Thr = new Double_t[nElem];

        // read LED threshold data
        if (!ReadParameters("taps_led1", fTAPSLED1Thr, nElem) ||
            !ReadParameters("taps_led2", fTAPSLED2Thr, nElem))
        {
            fTAPSled = kFALSE;
            return kFALSE;
        }

        Info("ApplyTAPScalib", "TAPS LED thresholds were read");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyPIDcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // PID phi calibration
    if (fPIDphi)
    {
        Double_t phi[nElem];
      
        // read phi data
        if (!ReadParameters("pid_phi", phi, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set phi
	    det->GetPosition(i)->SetZ(phi[i]);
        }
        
        Info("ApplyPIDcalib", "PID phi calibration was applied");
    }
    
    // PID droop correction
    if (fPIDdroop)
    {
        Char_t tmp[256];
        Double_t droop0[nElem];
        Double_t droop1[nElem];
        Double_t droop2[nElem];
        Double_t droop3[nElem];

        // read droop data
        if (!ReadParameters("pid_droop0", droop0, nElem) ||  
            !ReadParameters("pid_droop1", droop1, nElem) ||
            !ReadParameters("pid_droop2", droop2, nElem) ||
            !ReadParameters("pid_droop3", droop3, nElem))
        {
            fPIDdroop = kFALSE;
            return kFALSE;
        }
        
        // create droop functions array
        fNPIDDroopFunc = nElem;
        fPIDDroopFunc = new TF1*[fNPIDDroopFunc];

        // create droop functions
        for (Int_t i = 0; i < fNPIDDroopFunc; i++)
        {
            sprintf(tmp,"PIDDroopFunc_%d", i);
            fPIDDroopFunc[i] = new TF1(tmp, "pol3", 0, 180);
            fPIDDroopFunc[i]->SetParameters(droop0[i], droop1[i], droop2[i], droop3[i]);
        }

        Info("ApplyPIDcalib", "PID droop correction functions were created");
    }

    // PID energy calibration
    if (fPIDenergy)
    {
        Double_t e0[nElem];
        Double_t e1[nElem];
      
        // read ADC pedestal data
        if (!ReadParameters("pid_e0", e0, nElem)) return kFALSE;
        
        // read ADC gain data
        if (!ReadParameters("pid_e1", e1, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set ADC pedestal
	    det->GetElement(i)->SetA0(e0[i]);
            
            // set ADC gain
	    det->GetElement(i)->SetA1(e1[i]);
        }
        
        Info("ApplyPIDcalib", "PID energy calibration was applied");
    }

    // PID time calibration
    if (fPIDtime)
    {
        Double_t t0[nElem];
        
        // read TDC offset data
        if (!ReadParameters("pid_t0", t0, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set TDC offset
	    det->GetElement(i)->SetT0(t0[i]);
        }

        Info("ApplyPIDcalib", "PID time calibration was applied");
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyVetocalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // Veto energy calibration
    if (fVetoenergy)
    {
        Double_t e0[nElem];
        Double_t e1[nElem];
      
        // read ADC pedestal data
        if (!ReadParameters("veto_e0", e0, nElem)) return kFALSE;
        
        // read ADC gain data
        if (!ReadParameters("veto_e1", e1, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set ADC pedestal
	    det->GetElement(i)->SetA0(e0[i]);
            
            // set ADC gain
	    det->GetElement(i)->SetA1(e1[i]);
        }
        
        Info("ApplyVetocalib", "Veto energy calibration was applied");
    }

    // Veto time calibration
    if (fVetotime)
    {
        Double_t t0[nElem];
        Double_t t1[nElem];
        
        // read TDC offset data
        if (!ReadParameters("veto_t0", t0, nElem)) return kFALSE;
        
        // read TDC gain data
        if (!ReadParameters("veto_t1", t1, nElem)) return kFALSE;
        
        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{
	    // set TDC offset
	    det->GetElement(i)->SetT0(t0[i]);
            
            // set TDC gain
	    det->GetElement(i)->SetT1(t1[i]);
        }

        Info("ApplyVetocalib", "Veto time calibration was applied");
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ReadParameters(const Char_t* table, Double_t* par, Int_t length)
{
    // Read 'length' parameters of the calibration table 'table'
    // from the database to the value array 'par'.
    // Return kTRUE on success and kFALSE if an error occured.

    Char_t query[256];
    
    // check if database connection is open
    if (!IsConnected())
    {
        Error("ReadParameters", "No connection to database!");
        return kFALSE;
    }
    
    // create the query
    sprintf(query, 
	    "SELECT * FROM %s WHERE "
	    "calibration = '%s' AND "
	    "first_run <= %i AND "
	    "last_run >= %i "
            //"ORDER BY filled DESC LIMIT 1",
            "ORDER BY changed DESC LIMIT 1", // changed by PM 26.09.2011
	    table, fCalibration, fRun, fRun);
    
    // read from database
    TSQLResult* res = fDB->Query(query); 

    // check result
    if (!res)
    {
        Error("ReadParameters", "No calibration found for run %d in table %s!", fRun, table);
        return kFALSE;
    }
    else if (!res->GetRowCount())
    {
        Error("ReadParameters", "No calibration found for run %d in table %s!", fRun, table);
        delete res;
        return kFALSE;
    }

    // get data (parameters start at field 5)
    TSQLRow* row = res->Next();
    for (Int_t i = 0; i < length; i++) par[i] = atof(row->GetField(i+5));
  
    // clean-up
    delete row;
    delete res;
    
    return kTRUE;
}

