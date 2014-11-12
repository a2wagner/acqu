// SVN Info: $Id: LinkDef.h 1633 2013-01-24 14:00:17Z werthm $

/*************************************************************************
 * Author: Dominik Werthmueller, 2007-2010
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// LinkDef.h                                                            //
//                                                                      //
// OSCAR dictionary header file.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifdef __CINT__ 

// turn everything off
#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 
#pragma link off all typedef; 

#pragma link C++ nestedclasses; 
#pragma link C++ nestedtypedef; 

// Module A2
#pragma link C++ enum EA2Detector; 
#pragma link C++ typedef A2Detector_t; 
#pragma link C++ enum EA2TAPSType; 
#pragma link C++ typedef A2TAPSType_t; 
#pragma link C++ enum EA2TrigBitL1; 
#pragma link C++ typedef A2TrigBitL1_t; 
#pragma link C++ enum EA2TrigBitL2; 
#pragma link C++ typedef A2TrigBitL2_t; 
#pragma link C++ enum EA2BeamHelBit; 
#pragma link C++ typedef A2BeamHelBit_t; 
#pragma link C++ class TOA2Detector+; 
#pragma link C++ class TOA2ScalerEvent+; 
#pragma link C++ class TOA2Event+; 
#pragma link C++ class TOA2MCEvent+; 
#pragma link C++ class TOA2Trigger+; 
#pragma link C++ class TOA2TaggerHit+; 
#pragma link C++ class TOA2MCParticle+; 
#pragma link C++ class TOA2BaseDetParticle+; 
#pragma link C++ class TOA2DetParticle+; 
#pragma link C++ class TOA2RecParticle+; 
#pragma link C++ class TOA2RecMultiParticle+; 
#pragma link C++ class TOA2PSA+; 

// Module Analysis
#pragma link C++ namespace TOPhysics; 

// Module Main
#pragma link C++ namespace TOGlobals; 

// Module MC
#pragma link C++ class TOMCParticle+; 

// Module Utils
#pragma link C++ enum EONucleus; 
#pragma link C++ enum EOKinCut; 
#pragma link C++ typedef ONucleus_t; 
#pragma link C++ typedef OKinCut_t; 
#pragma link C++ namespace TOSUtils; 
#pragma link C++ namespace TOLoader; 
#pragma link C++ namespace TOFermiDistribution; 
#pragma link C++ class TOASCIIReader+; 
#pragma link C++ class TONumberReader+; 
#pragma link C++ class TOArray2D-;
#pragma link C++ class TOKinCut+;

// Module Math
#pragma link C++ class TOF1-;
#pragma link C++ class TOF1H+;
#pragma link C++ class TOdFCM+;
#pragma link C++ class TOFuncManager+;
#pragma link C++ class TOA2RecPi0Pi0+;
#pragma link C++ class TOF1Expo+;
#pragma link C++ class TOF1G+;
#pragma link C++ class TOF1GaussTail+;
#pragma link C++ class TOA2RecPi02g+;
#pragma link C++ class TOF1Empty+;
#pragma link C++ class TOA2RecEtaPi0+;
#pragma link C++ class TOF1BW_S11_1535+;
#pragma link C++ class TOA2RecEta2g+;
#pragma link C++ class TOA2RecEta6g+;
#pragma link C++ class TOF1BW+;
#pragma link C++ class TOF1Pol+;
#pragma link C++ class TOA2RecEtaPrime+;
#pragma link C++ class TOF1LegPol+;
#pragma link C++ class TOFCM+;
#pragma link C++ class TOF1Gauss+;
#pragma link C++ class TOPairCombinations+;
#pragma link C++ class TOA2RecK0Lambda+;


#endif

