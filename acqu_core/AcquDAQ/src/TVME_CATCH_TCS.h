//--Author	JRM Annand	 7th Jun 2007
//--Rev 	JRM Annand
//--Rev 	JRM Annand..28th Apr 2009..remove TDAQmemmap.h
//--Rev 	JRM Annand..25th Oct 2009..further development
//--Rev         B.Oussena   22nd Jul 2010 Added new cmd EcatchMapSize
//--Rev         B.Oussena                 fix an error kVME_CATCH_TCSKeys[]
//--Rev         B.Oussena                 set access to Virtual adr
//--Rev         B.Oussena                 modified PostInit()
//--Rev         JRM Annand  12th Sep 2010 Add mickey-mouse method ReConfig
//--Update      JRM Annand   6th Jul 2011 gcc4.6-x86_64 warning fix
//
//--Description
//                *** AcquDAQ++ <-> Root ***
// DAQ for Sub-Atomic Physics Experiments.
//
// TVME_CATCH_TCS
// COMPASS Trigger Controller for CATCH and GeSiCA
//

#ifndef __TVME_CATCH_TCS_h__
#define __TVME_CATCH_TCS_h__

#include "TVME_CATCH.h"

enum {
  ETCS_BinLine = 32,
  ETCS_MaxInitTry = 10000,
};
enum {
  ETCS_ConReg0 = 0xa,
  ETCS_ConReg1 = 0x40,
  ETCS_ConReg2 = 0x2,
  ETCS_Program = 0x1,
  ETCS_CClk = 0x2,
  ETCS_Data = 0x4,
  ETCS_Init = 0x1,
  ETCS_Done = 0x2,
  ETCS_Low = 0x0,
  ETCS_High = 0x1
};

class TVME_CATCH_TCS : public TVME_CATCH {
 protected:
  Char_t* fTCSfile;
  Char_t* fTCSfile1;
  Char_t* fTCSfile2;
 public:
  TVME_CATCH_TCS( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_CATCH_TCS();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  virtual Int_t ProgTCS( Char_t*, UInt_t );
  virtual Bool_t* FillBuffer(Char_t*, Int_t*);
  virtual void SetRunMode(UInt_t runmode);  
  virtual void ReConfig();
  ClassDef(TVME_CATCH_TCS,1)   
    };


#endif
