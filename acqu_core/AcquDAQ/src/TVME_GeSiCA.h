//--Author      JRM Annand    19th Jun 2003
//--Rev         JRM Annand...
//--Rev         ss,dk         When ??
//--Rev         JRM Annand     7th Dec 2003
//--Rev         D. Krambrich  17th Jan 2004
//--Rev         JRM Annand    30th Apr 2011
//--Rev         B. Oussena    10th Jun 2011  removed htons I2cWriteChk()
//                                           usleep() after a GeSiCA write
//--Rev         B. Oussena    13th Jun 2011  split (if Read(GeSiCA-Reg) & Value)
//--Rev         B. Oussena    27th Jun 2011  sleep(2) in TVME_GeSiCA::ProgSADC
//--Rev         B. Oussena    4th  Jul 2011  usleep(2) in TVME_GeSiCA::ProgSADC
//--Rev         B. Oussena    5th  Jul 2011  split (if Read(GeSiCA-Reg) & Value)
//                                           TVME_GeSiCA::SpyRead()
//--Rev         JRM Annand    9th  Jul 2011  Tidy up and further delay fiddle
//--Rev         JRM Annand   10th  Jul 2011  Use EIBase register Spy reset
//--Rev         JRM Annand   11th  Jul 2011  Try to improve speed
//--Rev         JRM Annand   24th  Jan 2012  Save TCS event ID
//--Rev         JRM Annand   25th  Jan 2012  Constructor set bits = 13
//--Rev         B. Oussena   22nd  Nov 2012  Add Send Event ID in SpyRead()
//--Rev         JRM Annand   17th  Sep 2013  Spy buff timeout 200 us
//--Update      JRM Annand   27th  Sep 2013  Try pause() spy buff wait
//
//--Description
//
//                *** Acqu++  ***
// Data Acquisition for Sub-Atomic Physics Experiments.
//
// GeSiCA
// Procedures for GeSiCA driver boards for SG-ADC
// Converted from ACQU 3v7 to AcquDAQ format JRMA 30th April 2011
// Initialisation Debugging by B. Oussena June-July 2011
//

#ifndef TVME_GeSiCA_h_
#define TVME_GeSiCA_h_

#include "TVME_CATCH.h"

//
// Options for SG-ADC setup
//
enum{ ESG_Single=1, ESG_Triple=3, ESG_Sample=4, ESG_All=5, ESG_Sum1=6, 
      ESG_Sum2=7 };

enum{ ENADCBoards=300, EGeSiCAFile, EGeSiCAReg, EGeSiCAMode, 
      ESADCReg, ESADCFile, ESADCThresh, ESADCSum, EGeSiCALog };
// 32 channel SADC has 32/2 = 16 channels per block
enum{ EChanPerBlock = 16 };

// Size of internal data buffer
enum{ ESizeData = 16384, EGeSiCATimeout = 200, EI2CTimeout = 2000000 }; 

enum{ EErrBufferEmpty = 1, EErrNoSLinkStart = 2, EErrSLinkErr = 3,
      EErrBuffersizeMismatch = 4, EErrNoCatchTrailer = 5};
// Indices of register list
// see VMEreg_t GeSiCAReg[]
enum { EIBase, EIFPGAReg, EIMStatus, EIDStatus, EIDatum, EIHlPort,
       EI2CAddrCtl, EI2CStatus, EI2CdwLow, EI2CdwHigh, EI2CdrLow,
       EI2CdrHigh, EI2CPortSelect, EIGenLD8 };



//-----------------------------------------------------------------------
// Register Structure of GeSiCA board (SG-ADC readout)
//----------------------------------------------------------------------

class TVME_GeSiCA : public TVME_CATCH {
 protected:
  Int_t fNSample;
  Int_t fNSampleWd;
  Int_t fOpMode;
  Int_t fADCmode;
  Int_t fLatency;
  Double_t fPedFactor;
  Int_t fNSADC;
  Int_t fNChip;
  Int_t fNthresh;
  Char_t* fSADCfile;
  Int_t* fSADCport;
  Int_t* fSADCchan;
  Int_t* fSADCthresh;
  UInt_t fSamStart[3];
  UInt_t fSamWidth[3];
 public:
  TVME_GeSiCA( Char_t*, Char_t*, FILE*, Char_t* );
  virtual ~TVME_GeSiCA();
  virtual void SetConfig( Char_t*, Int_t );   // configure VME module
  virtual void PostInit( );
  void ReadIRQ( void** );
  void ReadIRQorig( void** );
  virtual void SpyReset( );
  virtual Int_t SpyRead( void** );
  Char_t* LoadRBT(char*, Int_t*);
  void ProgFPGA();
  void ProgSADC();
  void ProgOpMode();
  void ProgSampleSum();
  void ProgThresh();
  void ProgramSgFPGA();
  void i2cWait();
  void i2cCheckACK();
  void i2cReset();
  void i2cSwitchPort(UInt_t);
  void i2cEnablePort(UInt_t);     
  void i2cEnableAll();
  void i2cDisablePort(UInt_t);
  void i2cDisableAll();
  void i2cWrite(UShort_t, UShort_t);
  void i2cWriteB(UShort_t, UShort_t);
  UShort_t i2cRead(UShort_t);
  UShort_t i2cReadB(UShort_t);
  void i2cWriteChk(Int_t, Int_t, UInt_t );
  void pause(){ Int_t j=0; for(Int_t i=0; i<10000; i++)j++; }

  ClassDef(TVME_GeSiCA,1)   
    };

//---------------------------------------------------------------------------
inline Int_t TVME_GeSiCA::SpyRead( void** outBuffer )
{
  // Read the GeSiCA spy (VMEbus) buffer
  // Check that the module is active..ie has received a trigger
  // If any error detected return 0 words read, reset spy buffer
  // and write error block into data stream
  
  UInt_t datum = 0;
  UInt_t* pStatus = (UInt_t*)fReg[EIDStatus];
  for( Int_t k=0; k<=EGeSiCATimeout; k++ ){  // less equal is important here!
    datum = *pStatus;
    if( datum & 0x1 )break;
    if( k == EGeSiCATimeout ){
      ErrorStore(outBuffer, 7);                 // error code 7
      SpyReset();
      return 0;
    }
  }
  UInt_t* pDatum = (UInt_t*)fReg[EIDatum];
  
  // Examine header, the first data word...should be 0x0
  UInt_t header = *pDatum;
  if( header != 0 ){
    ErrorStore( outBuffer, 1 );                // error code 1
    SpyReset();
    return 0;
  }
  
  // Check for error flag in 1st non-zero header datum
  header = *pDatum;
  if( header & ECATCH_ErrFlag ){
    ErrorStore( outBuffer,0x1 );
    SpyReset();
    return 0;
  }
  
  // Check consistency of data-buffer header and data-status register
  UInt_t nWordHeader = header & 0xffff;
  UInt_t nWordStatus = 0;
  UInt_t nWordTries = 0;
  do {
    nWordStatus = (*pStatus >> 16) & 0xffff;    
    if(nWordTries>200) {
      // reached maximum number of tries, this is 
      // the NEW error 4 (previously it only checked one time...)
      ErrorStore( outBuffer, 4 );                  // error code 4
      SpyReset();
      return 0;
    }
    nWordTries++;
  }
  while(nWordHeader != nWordStatus);
  
  // Check too many data words
  if( nWordHeader > fMaxSpy ){                            // overflow ?
    ErrorStore( outBuffer, 2 );                     // error code 2
    SpyReset();
    return 0;
  }
  
  // Make reads from buffer until we 
  // find the magic ECATCH_Trailer word, 
  // or fail after fMaxSpy reads
  UInt_t nWordRead = 0;
  while(true) {
    datum = *pDatum;
    fSpyData[nWordRead] = datum;
    nWordRead++;
    // be careful not to use *pDatum but datum, 
    // since this reads another word!
    if(datum == ECATCH_Trailer) {
      break;
    }
    
    
    if(nWordRead == fMaxSpy) {
      // this is a pretty bad error, 
      // one should probably better stop the DAQ somehow if this occurs...
      ErrorStore( outBuffer, 5 );                     // error code 5
      SpyReset();
      return 0;
    }
    
  }
  
  // Check if expected number of words match
  // number of words in spybuffer
  if(nWordRead != nWordHeader) {
    ErrorStore( outBuffer, 6 );                     // error code 6
    SpyReset();
    return 0;
  }
  
  
  // Check and buffer status reg. is 0
  datum = *pStatus;
  if( datum != 0 ) { 
    // we should make sure that the buffer is empty
    // so read until the status reg becomes zero
    UInt_t n = 0;
    do {
      datum = *pDatum;
      n++;
      if(n == fMaxSpy) {
        // mark this as an error
        ErrorStore( outBuffer, 3 );                     // error code 8
        SpyReset();
        return 0;
      }
    }
    while( *pStatus != 0);
    // mark this as an error
    ErrorStore( outBuffer, 3 );                     // error code 3
    SpyReset();
    return 0;
  }
  
  // Don't forget this: without saving the event ID there is no synch
  fTCSEventID = fSpyData[0];                      // 1st word = TCS event ID

  // if the system specifies that this GeSiCA sends the event ID to a remote 
  // system do it here
  if( fEventSendMod ) 
    fEventSendMod->SendEventID( fTCSEventID );  

  // no errors at all, return number of read words from SpyBuffer
  return (Int_t)nWordRead;
}

//---------------------------------------------------------------------------
inline void TVME_GeSiCA::SpyReset( )
{
  // Reset the GeSiCA module
  // Igor tells us, not to do this because he does
  // not know, what happens to the onboard TCS-Receiver. But anyway,
  // to have an empty buffer after BOS/EOS-procedures, I clear it, to
  // start catching "good data" in later events. We will need a "read the
  // buffer empty"-procedure on BOS/EOS.
  Write(EIBase,1);
  return;
}

#endif
