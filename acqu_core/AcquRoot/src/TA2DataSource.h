//--Author	JRM Annand	12th Jan 2004   Convert Acqu++ -> AcquRoot
//--Rev 	JRM Annand
//--Update	JRM Annand      20th Nov 2007   Add is process flag
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2DataSource
//  Basis of all online/offline sources for AcquRoot data streams
//
//---------------------------------------------------------------------------

#ifndef _TA2DataSource_h_
#define _TA2DataSource_h_

#include "TA2System.h"
//#include "TA2DataFormat.h"
#include "TA2RingBuffer.h"
#include "TThread.h"                   // for threaded event loop
//#include "TCondition.h"                // thread synch
#include "SwapBytes.h"

//class Semaphore_t;

class TA2DataSource : public TA2System {
protected:
  TThread* fSourceThread;       // To run data source as separate thread
  TA2RingBuffer* fBuffer;       // To run data source as separate thread
  //  Semaphore_t* fSem;
  //  TCondition* fSourceCond;
  //  TCondition* fServerCond;
  char** fFileName;             // File names
  UInt_t* fStartList;		// 1st record to analyse
  UInt_t* fStopList;		// last record to analyse
  UInt_t fStart;		// 1st record to analyse
  UInt_t fStop;		        // last record to analyse
  Int_t fInRecLen;		// no. of bytes in record
  Int_t fInPath;
  Int_t fOutPath;
  Int_t fNfile;			// counter number files processed
  UInt_t fNrecord;		// current record number
  Int_t fNbyte;			// no. of bytes in current record
  Int_t fRingSize;              // # elements in ring buffer
  //  Int_t fSemSource;             // semaphore id's
  //  Int_t fSemServer;
  Bool_t fIsStore;		// store data ?
  Bool_t fIsSwap;               // byte swapping?
  Bool_t fIsProcess;            // run a process?

public:
  TA2DataSource( const char*, Bool_t=EFalse, Int_t = 1 );
  // sets up inter-process communications
  virtual ~TA2DataSource();
  //  void UpdateHeader( char* = "/" );// new file...update info
  //  virtual void Initialise( Semaphore_t*, Int_t );
  virtual void Initialise( );
  virtual void Process()=0;	   // Process and send data
  virtual void Shutdown() { }	   // Shut down the data source
                                   // data-file info from file
  virtual void InputList(char*, UInt_t, UInt_t){ }
  virtual void SetConfig(char*, Int_t){ }
  virtual void Start();
  virtual void Run();		   // contains all the work

  Int_t GetInRecLen(){ return fInRecLen; }
  Int_t GetNfile(){ return fNfile; }
  Int_t GetNrecord(){ return fNrecord; }
  Int_t GetNbyte(){ return fNbyte; }
  Int_t GetRingSize(){ return fRingSize; }
  TA2RingBuffer* GetBuffer(){ return fBuffer; }    
  //  Int_t GetSemSource(){ return fSemSource; }             // semaphore id's
  //  Int_t GetSemServer(){ return fSemServer; }             // semaphore id's
  Bool_t IsStore(){ return fIsStore; }
  Bool_t IsSwap(){ return fIsSwap; }
  Bool_t IsProcess(){ return fIsProcess; }

  void SetIsStore(Bool_t state){ fIsStore = state; } // turn on/off storage
  void ClearList(){ fNfile = 0; }

  ClassDef(TA2DataSource,1)       	// for ROOT
};

#endif
