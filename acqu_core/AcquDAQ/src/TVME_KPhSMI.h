//--Author	JRM Annand   2nd Sep 2009  Adapt from acqu 3v7
//--Rev         B. Oussena   1st Dec 2009  Start Adapt from acqu 3v7
//--Rev         B. Oussena   25th Dec 2009 Add cmd MapReg, mod PostInit()
//--Rev         B. Oussena        Jan 2010 Add AddSMI()
//--Rev         B. Oussena        Feb 2010 Add Readout ADCs and Scalers
//--Rev         JRM Annand  27th Jul 2010  Clean up a bit
//--Rev         B. Oussena   6th Aug 2010  add  fNADC and fNScaler in ReadIRQ
//--Rev         JRM Annand   4th Sep 2010  ReadIRQScaler....the name
//--Rev         JRM Annand   9th Sep 2010  Slow read/write via usleep(1)
//--Rev         JRM Annand  14th Sep 2010  Scaler read..pipe suppress disabled
//--Update      JRM Annand  30th Mar 2014  Get rid of usleep in read/write
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_KPhSMI
// VMEbus interface to 8 control registers of LRS1821 SMI
// Interface can handle up to 16 SMIs
//

#ifndef __TVME_KPhSMI_h__
#define __TVME_KPhSMI_h__

#include "TVMEmodule.h"
#include "TFBmodule.h"

enum{ EKPhSMITimeout = 200 };

class TVME_KPhSMI : public TVMEmodule {
 protected:
  UShort_t** fSMIreg;     // Registers for currently selected SMI
  TList *fSmiList;        // List of SMIs controlled by SMI-Interface
  UInt_t* fSMIbuff;       // buffer to hold raw data from SMI
  Int_t fNsmi;            // No of SMIs
 public:
  TVME_KPhSMI( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_KPhSMI();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual void S_KPH(Int_t);
  virtual void R_KPH(Int_t, UShort_t*);
  virtual void W_KPH(Int_t, UShort_t);
  virtual void AddSMI(TDAQmodule*);
  virtual void ReadIRQ(void **);
  virtual void ReadIRQScaler(void  **);
  virtual void Pause(){ for(Int_t i=0; i<EKPhSMITimeout; i++){i++;i--;}}
  ClassDef(TVME_KPhSMI,1)   

    };

//-----------------------------------------------------------------------
inline void TVME_KPhSMI::S_KPH( Int_t smi )
{
  // Select SMI by offset in KPh register list
  fSMIreg = (UShort_t**)fReg + smi*8;
}

//---------------------------------------------------------------------
inline void TVME_KPhSMI::R_KPH(Int_t a, UShort_t *pdata)
{
  // Read SMI register a
  // a "pause" is necessary afterwards as the SMI-register bus is slow
  *pdata = *fSMIreg[a];
  Pause();
  //usleep(1);
}

//-----------------------------------------------------------------------
inline void TVME_KPhSMI::W_KPH( Int_t a, UShort_t data )
{  
  // Write SMI register a
  // a "pause" is necessary afterwards as the SMI-register bus is slow
  *fSMIreg[a] = data;
  Pause();
  //usleep(1);
}

#endif
