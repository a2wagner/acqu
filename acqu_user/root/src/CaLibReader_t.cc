/*************************************************************************
 * Authors: Irakli Keshelashvili, 2010
 *          Dominik Werthmueller, 2010-2012
 *          Thomas Strub, 2013
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
    fWasError = kFALSE;
    fRun = run;
    strcpy(fCalibration, calibration);
    fTAGGtime = kFALSE;
    fTAGGeff = kFALSE;
    fCBenergy = kFALSE;    
    fCBquadEnergy = kFALSE;
    fCBtime = kFALSE;      
    fCBwalk = kFALSE;     
    fCBled = kFALSE;
    fTAPSenergy = kFALSE;  
    fTAPSquadEnergy = kFALSE;
    fTAPStime = kFALSE;    
    fTAPSled = kFALSE;
    fTAPScfd = kFALSE;
    fPIDphi = kFALSE;   
    fPIDenergy = kFALSE;   
    fPIDtime = kFALSE;     
    fVetoenergy = kFALSE;  
    fVetotime = kFALSE;    
    fVetoled = kFALSE;

    fBadScRlist = kFALSE;
    fBadScRfirst = kFALSE;
    fBadScRlast = kFALSE;

    fTAGGEffVal = 0;
    fCBQuadEnergyPar0 = 0;
    fCBQuadEnergyPar1 = 0;
    fCBLEDThr = 0;
    fTAPSQuadEnergyPar0 = 0;
    fTAPSQuadEnergyPar1 = 0;
    fTAPSLED1Thr = 0;
    fTAPSLED2Thr = 0;

    fTotNScR = 0;
    fNBadScR = 0;
    fBadScR = 0;
    
    // connect to the SQL server
    Char_t tmp[256];
    sprintf(tmp, "mysql://%s/%s", dbHost, dbName);
    Info("CaLibReader_t", "Trying to connect to the database server '%s' ...", tmp);
    while (!(fDB = TSQLServer::Connect(tmp, dbUser, dbPass))) gSystem->Sleep(5000); 

    // check the connection to the server
    if (IsConnected()) Info("CaLibReader_t", "Successfully connected to CaLib database server");
    else 
    {
        Error("CaLibReader_t", "Could not connect to CaLib database server!");
        fWasError = kTRUE;
    }
}

//______________________________________________________________________________
CaLibReader_t::~CaLibReader_t() 
{
    // Destructor.
   
    if (fDB) delete fDB;
    if (fTAGGEffVal) delete [] fTAGGEffVal;
    if (fCBQuadEnergyPar0) delete [] fCBQuadEnergyPar0;
    if (fCBQuadEnergyPar1) delete [] fCBQuadEnergyPar1;
    if (fCBLEDThr) delete [] fCBLEDThr;
    if (fTAPSQuadEnergyPar0) delete [] fTAPSQuadEnergyPar0;
    if (fTAPSQuadEnergyPar1) delete [] fTAPSQuadEnergyPar1;
    if (fTAPSLED1Thr) delete [] fTAPSLED1Thr;
    if (fTAPSLED2Thr) delete [] fTAPSLED2Thr;
    if (fBadScR) delete [] fBadScR;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ReadScalerReads()
{
    // Sets the total number of scaler reads 'fTotNScR' from the database value
    // 'scr_n'. Fills the bad scaler read list 'fBadScR' with the database
    // value 'scr_bad'. If the flag 'fBadScRfirst' is set '0' is added to the
    // list, if the flag 'fBadScRlast' is set fTotNScR is added to the list.
    // Finally the number of bad scaler reads 'fNBadScR' will be set. Returns
    // kTRUE on success, kFALSE otherwise.

    // reset members
    fTotNScR = 0;
    fNBadScR = 0;
    if (fBadScR) delete [] fBadScR;
    fBadScR = 0;

    // init helpers
    Char_t tmp[128];
    Int_t nbadscr_tmp = 0;
    Int_t badscr_tmp[1204];

    // check database connection
    if (!IsConnected())
    {
        Error("ReadScalerReads", "BadScR calibration was not applied! No connection to database!");
        fWasError = kTRUE;
        return kFALSE;
    }

    // get total number of scaler reads from database
    if (!SearchRunEntry("scr_n", tmp))
    {
        // total number of scaler reads is neccessary for option fBadScRlast
        if (fBadScRlast)
        {
            Error("ReadScalerReads", "Could not find number of scaler reads in database.");
            fWasError = kTRUE;
            return kFALSE;
        }
        else
        {
            Warning("ReadScalerReads", "Could not find number of scaler reads in database.");
        }
    }
    else
    {
        // set number of scaler reads
        fTotNScR = atoi(tmp);
    }

    // add first to list
    if (fBadScRfirst)
    {
        badscr_tmp[nbadscr_tmp++] = 0;
        Info("ReadScalerReads", "Added first scaler read ('0') to bad scaler reads list.");
    }

    // add list to list
    if (fBadScRlist) 
    {
        // init bad scaler read input string
        Char_t badscr_str[1024];

        // read list
        if (!SearchRunEntry("scr_bad", badscr_str))
        {
            Error("ReadScalerReads", "Could not load bad scaler reads list from database.");
            fWasError = kTRUE;
            return kFALSE;
        }

        // get first token (comma separated)
        Char_t* t = strtok(badscr_str, ",");

        // reset counter if '0' is already in list
        if (fBadScRfirst && t && atoi(t) == 0) nbadscr_tmp = 0;

        // loop over tokens (comma separated)
        while (t)
        {
            // detect series (colon separated)
            Char_t* loc = strchr(t, (Int_t) ':');

            // check for series
            if (!loc)
            {
               // set single bad scr
               badscr_tmp[nbadscr_tmp++] = atoi(t);
            }
            else
            {
                // get start value of series
                Int_t min = atoi(t);

                // get stop value of series
                Int_t max = atoi(&loc[1]);
 
                // fill series to list
                for (Int_t i = min; i <= max; i++)
                    badscr_tmp[nbadscr_tmp++] = i;
            }

            // get next token
            t = strtok(0, " ,");
        }

        Info("ReadScalerReads", "Added bad scaler reads list from database to bad scaler reads list.");
    }

    // add last to list
    if (fBadScRlast)
    { 
        if (!nbadscr_tmp || badscr_tmp[nbadscr_tmp - 1] != fTotNScR) badscr_tmp[nbadscr_tmp++] = fTotNScR;
        Info("ReadScalerReads", "Added last scaler read ('%d') to bad scaler reads list.", fTotNScR);
    }

    // copy result to member
    fNBadScR = nbadscr_tmp;
    fBadScR = new Int_t[fNBadScR];
    for (Int_t i = 0; i < fNBadScR; i++) fBadScR[i] = badscr_tmp[i];

    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTAGGcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
 
    // check database connection
    if (!IsConnected())
    {
        Error("ApplyTAGGcalib", "TAGG calibration was not applied! No connection to database!");
        DisableTAGGcalib();
        fWasError = kTRUE;
        return kFALSE;
    }

    // get the number of detector elements
    Int_t nElem = det->GetNelement();
    
    // tagger time calibration
    if (fTAGGtime)
    {
        Double_t t0[nElem];
        
        // read TDC offset data
        if (!ReadParameters("tagg_t0", t0, nElem)) 
        {
            DisableTAGGcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
	{   
	    // set TDC offset
	    det->GetElement(i)->SetT0(t0[i]);
	}
      
        Info("ApplyTAGGcalib", "TAGG time calibration was applied");
    }
    
    // tagging efficiency
    if (fTAGGeff)
    {
        // create arrays
        if (!fTAGGEffVal) fTAGGEffVal = new Double_t[nElem];

        // read tagging efficiency data
        if (!ReadParameters("tagg_eff", fTAGGEffVal, nElem))
        {
            DisableTAGGcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        Info("ApplyTAGGcalib", "Tagging efficiency values were read");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyCBcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // check database connection
    if (!IsConnected())
    {
        Error("ApplyCBcalib", "CB calibration was not applied! No connection to database!");
        DisableCBcalib();
        fWasError = kTRUE;
        return kFALSE;
    }
 
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // CB energy calibration
    if (fCBenergy)
    {
        Double_t e1[nElem];
      
        // read SADC gain data
        if (!ReadParameters("cb_e1", e1, nElem)) 
        {
            DisableCBcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
        if (!ReadParameters("cb_t0", t0, nElem)) 
        {
            DisableCBcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
        if (!ReadParameters("cb_walk0", tw0, nElem) ||
            !ReadParameters("cb_walk1", tw1, nElem) ||
            !ReadParameters("cb_walk2", tw2, nElem) ||
            !ReadParameters("cb_walk3", tw3, nElem))
        {
            DisableCBcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
            DisableCBcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        Info("ApplyCBcalib", "CB quadratic energy correction was read");
    }
    
    // CB LED thresholds
    if (fCBled)
    {
        // create arrays
        if (!fCBLEDThr) fCBLEDThr = new Double_t[nElem];

        // read LED threshold data
        if (!ReadParameters("cb_led", fCBLEDThr, nElem))
        {
            DisableCBcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        Info("ApplyCBcalib", "CB LED thresholds were read");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyTAPScalib(TA2TAPS_BaF2* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // check database connection
    if (!IsConnected())
    {
        Error("ApplyTAPScalib", "TAPS calibration was not applied! No connection to database!");
        DisableTAPScalib();
        fWasError = kTRUE;
        return kFALSE;
    }
  
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // TAPS energy calibration
    if (fTAPSenergy)
    {
        Double_t e0_lg[nElem], e0_sg[nElem];
        Double_t e1_lg[nElem], e1_sg[nElem];
      
        // read energy data
        if (!ReadParameters("taps_lg_e0", e0_lg, nElem) ||
            !ReadParameters("taps_sg_e0", e0_sg, nElem) ||
            !ReadParameters("taps_lg_e1", e1_lg, nElem) ||
            !ReadParameters("taps_sg_e1", e1_sg, nElem))
        {
            DisableTAPScalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
        
        // read TDC offset and gain data
        if (!ReadParameters("taps_t0", t0, nElem) ||
            !ReadParameters("taps_t1", t1, nElem)) 
        {
            DisableTAPScalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
            DisableTAPScalib();
            fWasError = kTRUE;
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
            DisableTAPScalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        Info("ApplyTAPScalib", "TAPS LED thresholds were read");
    }
    
    // TAPS CFD thresholds
    if (fTAPScfd)
    {
        Double_t cfd[nElem];
      
        // read CFD threshold data
        if (!ReadParameters("taps_cfd", cfd, nElem)) 
        {
            DisableTAPScalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set CFD threshold
	    det->GetElement(i)->SetEnergyLowThr(cfd[i]);
        }

        Info("ApplyTAPScalib", "TAPS CFD thresholds were applied");
    }
  
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::ApplyPIDcalib(TA2Detector* det)
{
    // Read calibration from database using the detector pointer 'det'.
    // Return kTRUE on success and kFALSE if an error occured.
    
    // check database connection
    if (!IsConnected())
    {
        Error("ApplyPIDcalib", "PID calibration was not applied! No connection to database!");
        DisablePIDcalib();
        fWasError = kTRUE;
        return kFALSE;
    }
     
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // PID phi calibration
    if (fPIDphi)
    {
        Double_t phi[nElem];
      
        // read phi data
        if (!ReadParameters("pid_phi", phi, nElem)) 
        {
            DisablePIDcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set phi
	    det->GetPosition(i)->SetZ(phi[i]);
        }
        
        Info("ApplyPIDcalib", "PID phi calibration was applied");
    }
    
    // PID energy calibration
    if (fPIDenergy)
    {
        Double_t e0[nElem];
        Double_t e1[nElem];
      
        // read energy data
        if (!ReadParameters("pid_e0", e0, nElem) ||
            !ReadParameters("pid_e1", e1, nElem)) 
        {
            DisablePIDcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
        if (!ReadParameters("pid_t0", t0, nElem)) 
        {
            DisablePIDcalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
    
    // check database connection
    if (!IsConnected())
    {
        Error("ApplyVetocalib", "Veto calibration was not applied! No connection to database!");
        DisableVetocalib();
        fWasError = kTRUE;
        return kFALSE;
    }
  
    // get the number of detector elements
    Int_t nElem = det->GetNelement();

    // Veto energy calibration
    if (fVetoenergy)
    {
        Double_t e0[nElem];
        Double_t e1[nElem];
      
        // read energy data
        if (!ReadParameters("veto_e0", e0, nElem) ||
            !ReadParameters("veto_e1", e1, nElem)) 
        {
            DisableVetocalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
        
        // read time data
        if (!ReadParameters("veto_t0", t0, nElem) ||
            !ReadParameters("veto_t1", t1, nElem)) 
        {
            DisableVetocalib();
            fWasError = kTRUE;
            return kFALSE;
        }

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
    
    // Veto LED thresholds
    if (fVetoled)
    {
        Double_t led[nElem];
      
        // read LED threshold data
        if (!ReadParameters("veto_led", led, nElem)) 
        {
            DisableVetocalib();
            fWasError = kTRUE;
            return kFALSE;
        }

        // loop over detector elements
        for (Int_t i = 0; i < nElem; i++)
        {
            // set LED threshold
	    det->GetElement(i)->SetEnergyLowThr(led[i]);
        }

        Info("ApplyVetocalib", "Veto LED thresholds were applied");
    }
 
    return kTRUE;
}

//______________________________________________________________________________
Bool_t CaLibReader_t::SearchRunEntry(const Char_t* name, Char_t* out)
{
    // Search the information 'name' for the run 'fRun' and write it to 'out'.
    // Return kTRUE when the information was found, otherwise kFALSE.

    Char_t query[256];

    // create the query
    sprintf(query,
            "SELECT %s FROM run_main "
            "WHERE run = %d",
            name, fRun);

    // read from database
    TSQLResult* res = fDB->Query(query);

    // check result
    if (!res)
    {
        Error("SearchRunEntry", "Could not find the information '%s' for run %d!", name, fRun);
        return kFALSE;
    }
    if (!res->GetRowCount())
    {
        Error("SearchRunEntry", "Could not find the information '%s' for run %d!", name, fRun);
        delete res;
        return kFALSE;
    }

    // get row
    TSQLRow* row = res->Next();

    // write the information
    const Char_t* field = row->GetField(0);
    if (!field) strcpy(out, "");
    else strcpy(out, field);

    // clean-up
    delete row;
    delete res;

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
            "ORDER BY changed DESC LIMIT 1",
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

//______________________________________________________________________________
void CaLibReader_t::DisableTAGGcalib()
{
    // Disable all tagger calibrations.

    fTAGGtime = kFALSE;
    fTAGGeff = kFALSE;
}

//______________________________________________________________________________
void CaLibReader_t::DisableCBcalib()
{
    // Disable all CB calibrations.

    fCBenergy = kFALSE;
    fCBquadEnergy = kFALSE;
    fCBtime = kFALSE;
    fCBwalk = kFALSE;
    fCBled = kFALSE;
}

//______________________________________________________________________________
void CaLibReader_t::DisableTAPScalib()
{
    // Disable all TAPS calibrations.

    fTAPSenergy = kFALSE;
    fTAPSquadEnergy = kFALSE;
    fTAPStime = kFALSE;
    fTAPSled= kFALSE;
}

//______________________________________________________________________________
void CaLibReader_t::DisablePIDcalib()
{
    // Disable all PID calibrations.

    fPIDphi = kFALSE;
    fPIDenergy = kFALSE;
    fPIDtime = kFALSE;
}

//______________________________________________________________________________
void CaLibReader_t::DisableVetocalib()
{
    // Disable all Veto calibrations.

    fVetoenergy = kFALSE;
    fVetotime = kFALSE;
}

