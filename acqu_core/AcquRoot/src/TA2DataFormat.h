//--Author	JRM Annand   12th Jan 2004   Convert Acqu++ -> AcquRoot
//--Rev 	JRM Annand
//--Rev 	JRM Annand   22nd Apr 2004   Bug-fix header merge
//--Rev	        JRM Annand   12th Apr 2006   Add SetEventIndex(Int_t)
//--Rev 	JRM Annand    5th Nov 2006   fStartOffset, fNsrc
//--Rev 	JRM Annand...22nd Jan 2007...4v0 update
//--Rev 	JRM Annand...27th Apr 2007...Flexi header info getters
//--Rev 	JRM Annand...11th Jun 2007...Error sort flag
//--Rev 	JRM Annand...21st Apr 2008...UpdateInfo(),ConstructHeader()
//--Rev 	JRM Annand...21st Nov 2012...Stuff for Mk2 merging
//--Update	JRM Annand... 8th Feb 2013...Add GetNEPICS() prototype
//--Update  A Neiser...   6th June 2013  Make char_t* const
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
//	TA2DataFormat
//	Basis of all event-by-event data stream formats
//	recognised by AcquRoot for Data Serving purposes
//
//---------------------------------------------------------------------------

#ifndef _TA2DataFormat_h_
#define _TA2DataFormat_h_

#include "TA2System.h"
#include "DataFormats.h"		// specifies data structure
#include "TA2RingBuffer.h"
#include "SwapBytes.h"
#include "ModuleIndex.h"

// Foreign data setup file keys
enum { EForeignTot, EForeignModule };

class TA2DataFormat : public TA2System {
protected:
  TA2RingBuffer* fIn;		// -> buffer where data is read in
  TA2RingBuffer* fOut;		// -> buffer where data is transferred
  void* fHeader;	 	// -> ACQU-style header buffer
  void* fPartEvent;             // partial event holder 
  Char_t* fDataFileName;	// name of data file
  Int_t fNevent;	       	// event counter
  Int_t fRecLen;		// Input record length
  Int_t fOnline;	      	// online/offline data taking
  Int_t fNadc;			// no. of adc's (foreign formats)
  Int_t fNscaler;	       	// no, of scalers (foreign formats)
  Int_t fADCoffset;             // ADC offset (foreign formats)
  Int_t fReadyServer;           // DataServer ready semaphore
  Int_t fReadySort;             // Sort ready semaphore
  Int_t fPartEventSize;         // length (bytes) of any partial event
  UInt_t fStartOffset;          // allows for any header info at start of event
  Int_t fNsrc;                  // data source # (multi source running)
  Bool_t fIsHeader;             // is there a header record?
  Bool_t fIsSwap;		// Flag to determine byte order
  Bool_t fIsBoundary;           // formats that don't respect buff boundaries 
  Bool_t fIsPartEvent;          // event buffer overflow flag
  Bool_t fIsErrorSort;          // support for analysis of error blocks
  UShort_t fEventID;            // Event ID stamp
  UShort_t fEventIndex;         // "adc" index of marker channel
  UInt_t fEventNumber;          // Event Number  
public:
  TA2DataFormat(const Char_t*, Int_t, Int_t=0);
  virtual ~TA2DataFormat(){}

  virtual void Initialise( void*, TA2RingBuffer*, TA2RingBuffer*, UInt_t );
  virtual void CheckHeader(){}                  // Is the header buffer OK?
  virtual void InitForeign();			// for foreign formats
  virtual void ConstructHeader(void* = NULL){}	// ACQU header construction
  virtual void MergeHeader(void*, void*, Int_t, Int_t) = 0;
                                                // ACQU header merging
  virtual void UpdateHeader(void*, void*){};    // start of new file
  virtual void UpdateInfo(){}                   // header info update prototype
  virtual void SwapHeader(void*){}	        // byte-swap header buffer
  virtual void SendData();		        // process data buffer
  virtual void SendEOF();			// signal end-of-file
  virtual void SendTerminate();			// signal end-of-analysis
  virtual void SetFileName( AcquExptInfo_t*, const char* = "/", int = 0 );
  virtual char* GetFileName(){ return fDataFileName; }
  virtual Int_t HeaderSize( void*);		// no. bytes in header
  virtual Int_t GetEventLength( void* );

  virtual Bool_t IsHeadBuff( void* ) = 0;       // is buffer a header
  virtual Bool_t IsDataBuff( void* ) = 0;       // is buffer a data buffer
  virtual Bool_t IsEndBuff( void* ) = 0;        // is buffer a trailer
  virtual Bool_t IsTermBuff( void* ) = 0;       // is buffer a terminator
  virtual void SetHeader( void* b ){ *(UInt_t*)b = EHeadBuff; } // Set header
  virtual void SetData( void* b ){ *(UInt_t*)b = EDataBuff; }
  virtual void SetTrailer( void* b ){ *(UInt_t*)b = EEndBuff; }
  virtual void SetTerminate( void* b ){ *(UInt_t*)b = EKillBuff; }
  virtual UInt_t GetBufferType( void* b ){ return (UInt_t)EDataBuff; }

  Int_t GetReadyServer(){ return fReadyServer; }
  Int_t GetReadySort(){ return fReadySort; }
  Int_t GetADCoffset(){ return fADCoffset; }
  void* GetHeader(){ return fHeader; }
  UShort_t GetEventID(){ return fEventID; }
  UShort_t GetEventIndex(){ return fEventIndex; }
  void SetEventIndex( Int_t index ){ fEventIndex = index; }
  UInt_t GetEventNumber(){ return fEventNumber; }
  virtual UInt_t GetStartOffset(){ return fStartOffset; }
  Int_t GetNsrc(){ return fNsrc; }
  // Retrieve header info
  Bool_t IsMaybeMk2(); // looks for the Mk2 marker...
  virtual const Char_t* GetTime() {return "";}
  virtual const Char_t* GetDescription() {return "";}
  virtual const Char_t* GetRunNote() {return "";}
  virtual const Char_t* GetOutFile() {return "";}
  virtual const Int_t GetRunNumber() {return 0;}
  virtual Int_t GetNScaler(){ return 0; }
  virtual Int_t GetNADC(){ return 0; }
  virtual Int_t GetNEPICS(){ return 0; }
  virtual Int_t GetNModule(){ return 0; }
  virtual Int_t GetModID(Int_t i){ return 0; }
  virtual Int_t GetModIndex(Int_t i){ return 0; }
  virtual Int_t GetModType(Int_t i){ return EDAQ_Undef; }
  virtual Int_t GetModAmin(Int_t i){ return 0; }
  virtual Int_t GetModNChan(Int_t i){ return 0; }
  virtual Int_t GetModBits(Int_t i){ return 0; }
  virtual Int_t GetErrorSize( ){ return 0; }
  virtual UInt_t GetDataBuffHead(){ return EDataBuff; }
  //
  Bool_t IsHeader(){ return fIsHeader; }
  Bool_t IsSwap(){ return fIsSwap; }
  Bool_t IsBoundary(){ return fIsBoundary; }
  Bool_t IsPartEvent(){ return fIsPartEvent; }
  Bool_t IsErrorSort(){ return fIsErrorSort; }

  virtual void SetConfig(Char_t*, Int_t){ }
  ClassDef(TA2DataFormat,1)              	// for ROOT
};

//-----------------------------------------------------------------------------
inline int TA2DataFormat::HeaderSize(void* h)
{
  // Calculate the no. of bytes in an ACQU-standard header
  //
  AcquExptHeader_t* A = (AcquExptHeader_t*)((char*)h + sizeof(UInt_t));
  
  return( sizeof(AcquExptHeader_t)	+
	  sizeof(Pspect_t) * (A->n_spect + A->n_scaler) +
	  sizeof(ModHeader_t) * A->n_module +
	  sizeof(UInt_t) );
}

//-----------------------------------------------------------------------------
inline void TA2DataFormat::SendEOF()
{
  //  Write End-of-File marker to output buffer
  //
  fOut->WaitEmpty();
  *((UInt_t*)fOut->GetStore()) = EEndBuff;	// End-of-File flag
  fOut->FNext();
}

//-----------------------------------------------------------------------------
inline void TA2DataFormat::SendTerminate()
{
  // Write End-of-File marker to output buffer
  //
  fOut->WaitEmpty();
  *((UInt_t*)fOut->GetStore()) = EKillBuff;   // End-of-File flag
  fOut->FNext();
}

//-----------------------------------------------------------------------------
inline void TA2DataFormat::SendData( )
{
  // Copy input data buffer to output (sort) data buffer
  //
  fOut->WaitEmpty();
  memcpy( fOut->GetStore(), fIn->GetStore(), fIn->GetLenBuff() );
  fOut->FNext();
}

//-----------------------------------------------------------------------------
inline Int_t TA2DataFormat::GetEventLength( void* evbuff )
{
  // Identify event's - worth of data and return event length (# of bytes )
  // Mk1 format....if 0th data source account for event number in data
  //
  Int_t length = 0;
  UInt_t* event = (UInt_t*)evbuff;

  if( !fNsrc ){
    fEventNumber = *event;
    event++; length++;     // step past event #
  }
  else fEventNumber = *(event-1);
  if( fEventNumber == EEndEvent ){
    return 0;
  }

  // If final event read of scalers only
  // For now frig synch by discarding this event
  if( *event == EScalerBuffer ){
    return 0;
  }

  if( fEventIndex == *(UShort_t*)event )
    fEventID = *(UShort_t*)( (UShort_t*)event + 1 );
  else fEventID = 0xffff;

  while( *event != EEndEvent ){
    length++;    event++;
  }
  return (length * sizeof(UInt_t));
}

#endif
