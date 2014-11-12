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
       EI2CdrHigh, EI2CPortSelect };



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
  UInt_t fNSADC;
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
  virtual void SpyReset( );
  virtual Int_t SpyRead( void** );
  void ProgFPGA();
  void ProgSADC();
  void ProgOpMode();
  void ProgSampleSum();
  void ProgThresh();
  bool ProgramSgFPGA(std::vector<UInt_t> rbt_data);
  
  bool i2c_wait(bool check_ack);  
  bool i2c_reset();  
  void i2c_set_port(UInt_t port_id, bool broadcast = false);  
  bool i2c_read(UInt_t bytes, UInt_t addr, UInt_t& data);  
  bool i2c_write(UInt_t bytes, UInt_t addr, UInt_t data);  
  bool i2c_read_reg(UInt_t adc_side, UInt_t reg, UInt_t& data);  
  bool i2c_write_reg(UInt_t adc_side, UInt_t reg, UInt_t data, bool check = false);  
  bool load_rbt(const char* rbt_filename, std::vector<UInt_t>& data);  
  bool init_gesica(bool skip_i2c = false);
  
  ClassDef(TVME_GeSiCA,1)   
};



#endif
