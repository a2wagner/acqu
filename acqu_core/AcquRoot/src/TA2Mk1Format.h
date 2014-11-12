//--Author	JRM Annand	12th Jan 2004   Convert Acqu++ -> AcquRoot
//--Rev 	JRM Annand
//--Rev 	JRM Annand       5th Nov 2006   fStartOffset, nsrc
//--Rev 	JRM Annand      22nd Jan 2007   4v0 update
//--Rev 	JRM Annand      27th Apr 2007   Standard header interface
//--Rev 	JRM Annand      13th Jul 2007   Bug Mk1 merge header
//--Rev 	JRM Annand      29th Nov 2007   GetModNChan() returns fAmax + 1
//--Rev 	JRM Annand      12th Feb 2008   Check if module already ADC
//--Rev 	JRM Annand      21st Apr 2008   UpdateInfo(),ConstructHeader()
//--Update	JRM Annand       6th Sep 2008   UpdateInfo() m[j].fAmax
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
//	TA2Mk1Format
//	Specify standard Mk1 (old) ACQU formatted data stream
//
//---------------------------------------------------------------------------

#ifndef _TA2Mk1Format_h_
#define _TA2Mk1Format_h_

#include "TA2DataFormat.h"

class TA2Mk1Format : public TA2DataFormat {
 protected:
  AcquExptInfo_t* fExpInfo;   // experimental info from header
  ModuleInfo_t* fModInfo;     // experimental hardware modules
 public:
  TA2Mk1Format( Char_t*, Int_t, Int_t=0 );
  virtual ~TA2Mk1Format(){ }
  virtual void ConstructHeader(void* = NULL);
  virtual void MergeHeader(void*, void*, Int_t, Int_t);
  virtual void CheckHeader();
  virtual void SwapHeader(void*);	                // byte-swap header
  virtual void UpdateHeader(void*, void*);	        // start of new file
  virtual void UpdateInfo();	                        // new fModInfo
  virtual void ResetOut(){ }			        // reset output buffer

  virtual Bool_t IsHeadBuff( void* b ){         // Header buffer ?
    if( *(UInt_t*)b == EHeadBuff ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsDataBuff( void* b ){         // Data buffer ?
    if( *(UInt_t*)b == EDataBuff ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsEndBuff( void* b ){          // Trailer buffer ?
    if( *(UInt_t*)b == EEndBuff ) return ETrue;
    else return EFalse;
  }
  virtual Bool_t IsTermBuff( void* b ){         // Terminator buffer ?
    if( *(UInt_t*)b == EKillBuff ) return ETrue;
    else return EFalse;
  }
  virtual UInt_t GetBufferType( void* b ){ return *(UInt_t*)b; }

  // Retrieve header info
  virtual const Char_t* GetTime() {return fExpInfo->fTime;}
  virtual const Char_t* GetDescription() {return fExpInfo->fDescription;}
  virtual const Char_t* GetRunNote() {return fExpInfo->fRunNote;}
  virtual const Char_t* GetOutFile() {return fExpInfo->fOutFile;}
  virtual const Int_t GetRunNumber() {return fExpInfo->fRun;}
  virtual Int_t GetNScaler(){ return fExpInfo->fNscaler; }
  virtual Int_t GetNADC(){ return fExpInfo->fNspect; }
  virtual Int_t GetNModule(){ return fExpInfo->fNmodule; }
  //  virtual Int_t GetModID(Int_t i){ return (fModInfo[i]).fModType; }
  //  virtual Int_t GetModIndex(Int_t i){ return (fModInfo[i]).fStation; }
  virtual Int_t GetModType(Int_t i){ return (fModInfo[i]).fVicCrate; }
  virtual Int_t GetModAmin(Int_t i){
    if( fModInfo[i].fVicCrate == EForeignADC ) return fModInfo[i].fCrate;
    else return (fModInfo[i]).fAmin;
  }
  virtual Int_t GetModNChan(Int_t i){ 
    return fModInfo[i].fAmax + 1;
  }
  virtual Int_t GetModBits(Int_t i){ return (fModInfo[i]).fBits; }
  virtual Int_t GetErrorSize( ){ return sizeof(ReadError_t); }
  //  

  ClassDef(TA2Mk1Format,1)       		// for ROOT
};


#endif
