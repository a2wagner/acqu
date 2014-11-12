//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2KFit
// Kinematical fit routines
//

#include "TA2KFit.h"
//#include <iostream.h>

extern "C"
{
  extern Float_t Vmodf( Float_t*, Int_t );
  extern Double_t Vmod( Double_t*, Int_t );
  extern void Vsubf( Float_t*, Float_t*, Float_t*, Int_t );
  extern void Vsub( Double_t*, Double_t*, Double_t*, Int_t );
  extern void Ptpxyz( Double_t*, Double_t* ); 
  extern void Ptpxyzf( Float_t*, Float_t* ); 
  extern void Pxyztp( Double_t*, Double_t* );
  extern void Pxyztpf( Float_t*, Float_t* );
}			    

//---------------------------------------------------------------------------
TA2KFit::TA2KFit()
{
 //  fPhotons = new TLorentzVector[8];
  //fPulls = new Double_t[50];
  fTypDat=1;
}

//---------------------------------------------------------------------------
TA2KFit::~TA2KFit()
{
  // Free memory allocated to store parameters
  //  if( fPhotons ) delete fPhotons;
  //if( fPulls ) delete fPulls;
}

/*---------------------------------------------------------------------------
Float_t TA2KFit::dEovEProt(Float_t Ecl)
{
  Float_t Epr = EnProt(Ecl);
  Float_t dEcl = 0.011472 + 0.051288*Ecl;
  return dEcl/Epr;
}
---------------------------------------------------------------------------*/
Int_t TA2KFit::Kfilbm(Double_t Amastag, Int_t Ikind,
      Double_t Beampr[4], Double_t Sbeam[3], Double_t Target[2])
{
/*********************************************************************
     it serves to fill the beam information and have to be called first
       of the other subroutines.

       Input parameters:

        Npmea  - number of measured particles;
                 Npmea + number of decaying particles .LE.20 

        Amastag - mass of the target particle
                  (proton for the LH2 target)

        Ikind = GEANT index of the particle

        Beampr[4]  = mass, energy, x, y of the incident particle at z=0.

        Sbeam[3] = dE/E(beam energy), sigma(beam x), sigma(beam y).

        Target[2] = z coordinate of the target center, length of the target

        Idecfl - flag of secondary vertex: = 0 - no sec.vertex != 0 - there is one. 

   Function returns ierr parameters:

   ierr = 0 when all is correct; >0 in case of the error of the input parameters
************************************************************************/
 
  Int_t il, ierr, Npart;
  Float_t einv;
 
  ierr = 0;
  Npart = 0;
  fMastg = Amastag;
  fKind[Npart] = Ikind;
  fMass[Npart] = Beampr[0];
  einv = 1./(Float_t)Beampr[1];

  fPin[Npart][0] = einv;
  fPin[Npart][1] = (Float_t)Beampr[2];  //Xentr
  fPin[Npart][2] = (Float_t)Beampr[3];  //Yentr
  fPin[Npart][3] = (Float_t)Target[0];  //Zentr

  fErrSq[Npart][0] = (Float_t)(Sbeam[0]*Sbeam[0]*einv*einv);
  fErrSq[Npart][1] = (Float_t)(Sbeam[1]*Sbeam[1]);
  fErrSq[Npart][2] = (Float_t)(Sbeam[2]*Sbeam[2]);
  fErrSq[Npart][3] = (Float_t)(0.333*0.333*Target[1]*Target[1]);
  //fErrSq[Npart][3] = (Float_t)(0.4*0.4*Target[1]*Target[1]);

  il = Npart*4;
  fPlim[il][0] = 0.;
  fPlim[il][1] = 1000.;
  if ( fMass[Npart] > 0. ) {
    fPlim[il][0] = 1./(fMass[Npart]*100.);
    fPlim[il][1] = 1./fMass[Npart];
  }      
  fPlim[il+1][0] = 0.;
  fPlim[il+1][1] = 0.;
  fPlim[il+2][0] = 0.;
  fPlim[il+2][1] = 0.;
  fPlim[il+3][0] = -12.5-Target[1]*0.5;
  fPlim[il+3][1] =  12.5+Target[1]*0.5;

  return ierr;
}

//---------------------------------------------------------------------------
Int_t TA2KFit::Kfilcst(Int_t Npart, Int_t Ikind, Double_t Pacst[6] )
                    
{
/*************************************************************************
      serves to fill information about one cluster (usually from photon)
      in CB/TAPS and moreover to fill information about cluster
      from neutron/proton !!! in case it is measured !!!

       Input parameters:

    Npart - sequence number of the measured particle submitted to the KINFIT
            input;  must be > 0, as #0 is reserved  for the beam particle

    Ikind = GEANT index of the particle

    Pacst[5] = mass, kin. energy, theta, phi, cluster depth in the crystals,
               original cluster energy

    Scst[4] = de/e, sigma(theta), sigma(phi), sigma(depth)

   Function returns ierr parameters:

   ierr = 0 when all is correct; >0 in case of the error of the input parameters
************************************************************************/        
  Int_t il, ierr;
  Float_t cstmas, eng, einv, Ecl, Scst[4];

  ierr = 0;
  fMass[Npart] = Pacst[0];
  cstmas = (Float_t)Pacst[0];
  eng = (Float_t)Pacst[1];
  Ecl = (Float_t)Pacst[5];
  if ( Npart==1 ) {
    fLcst = 0;
    fLcstCB = 0;
    fLcstTAPS = 0;
  }
  fLcst++;
  if ( Npart < 1 )  return ierr=1;
  if ( Ikind > 0 ) {
    fCalor[Npart] = 1;  // CB cluster
    fLcstCB++;
    Scst[0] = dEovEclCB( Ecl, abs(Ikind) );         // photon dE/E
    Scst[1] = dThetaCB( Ecl , abs(Ikind) );         // dTheta
    Scst[2] = Scst[1]/sinf(Pacst[2]);         // dPhi
    Scst[3] = dDepthShowCB( Ecl , abs(Ikind) );         // dDepth
  }
  else {
    fCalor[Npart] = 2;  // TAPS cluster
    fLcstTAPS++;
    Scst[0] = dEovEclTAPS( Ecl, abs(Ikind) );         // photon dE/E
    //if ( abs(Ikind)==14 ) printf("%f %f %f\n",Ecl,eng,Scst[0]);
    Scst[1] = dTanThTAPS( Ecl, abs(Ikind) );         // dTheta
    Scst[2] = Scst[1]/Pacst[2];         // dPhi
    Scst[3] = dDepthShowTAPS( Ecl , abs(Ikind) );         // dDepth
  }
  if ( eng > 0. ) {
    fKind[Npart] = abs(Ikind);  //for a totally measured particle.
    einv = 1./eng;
  }
  else {
    fKind[Npart] =-abs(Ikind);  //for a case when the angles are measured only.
    einv = 0.;
  }
  fPin[Npart][0] = einv;
  fPin[Npart][1] = (Float_t)Pacst[2];  //theta for CB, tan(theta)*z_cluster for TAPS (i.e., radius)
  fPin[Npart][2] = (Float_t)Pacst[3];  //phi
  fPin[Npart][3] = (Float_t)Pacst[4];  //cluster depth in the crystals

  fErrSq[Npart][0] = (Float_t)(Scst[0]*Scst[0])*einv*einv;
  fErrSq[Npart][1] = (Float_t)(Scst[1]*Scst[1]);
  fErrSq[Npart][2] = (Float_t)(Scst[2]*Scst[2]);
  fErrSq[Npart][3] = (Float_t)(Scst[3]*Scst[3]);

  il = Npart*4;
  fPlim[il][0] = 0.;
  fPlim[il][1] = 200.;
  fPlim[il+1][0] = 0.;
  Float_t thetup = 3.141;
  if ( fCalor[Npart] == 2 )  thetup = 75.;
  fPlim[il+1][1] = thetup;
  fPlim[il+2][0] = 0.;
  fPlim[il+2][1] = 0.;
  if ( fCalor[Npart] == 1 ) {
     fPlim[il+3][0] = 20.;
     fPlim[il+3][1] = 80.;
  }
  if ( fCalor[Npart] == 2 ) {
     fPlim[il+3][0] = (Float_t)Pacst[4]-26.;
     fPlim[il+3][1] = (Float_t)Pacst[4]+26.;
  }
  return ierr;
}
  
//---------------------------------------------------------------------------
Int_t TA2KFit::Kfilunm(Int_t Npart, Int_t Ikind, Double_t Amas)
{
/************************************************************************
      serves to inform KINFIT about unmeasured particle

     Must be called after all KFILCST calls.

    Input parameters:

    Npart - sequence number of the particle submitted to KINFIT input;
            must be equal   2 + Number of clusters

    Ikind = GEANT index of the particle

    Amas - mass of the unmeasured hadron

   Function returns ierr parameters:

   ierr = 0 when all is correct; >0 in case of the error of the input parameters
************************************************************************/
  Int_t ierr;

  ierr = 0;
  fKind[Npart] = Ikind + 1000;      // FOR UNMEASURED PARTICLE
  fMass[Npart] = Amas;
  return ierr;
}

//---------------------------------------------------------------------------
Int_t TA2KFit::Kinfit(Int_t Nptall, Int_t Ndecay, Int_t Idecfl,
		       Int_t Ikind[NDMAX], Double_t Amsdec[NDMAX],
		       Int_t Jdecay[NDMAX], Int_t Ldecay[NDMAX],
		       Int_t Idecay[NDMAX][10], Int_t Ifzfree )
{
  Int_t i, j, k, ierr;
  Int_t Nunme, NY, NF, NFF, IV, ISV, I, J, IY, K, I2;
  Int_t Isecvt[NPMAX], lclsec, jdecfl, idechk, ideca, Lpart, Iunmea;
  Int_t IFIRSL, NZDIV, IZPOS, IfZfree;
  Int_t ICDEC, IFSCND, LLDEC, IDE, INDE, IEF; 

  Float_t Y[NXMAX], F[NFMAX], WBEST[2], ZBEST[2], ZBESTM;
  Float_t ZV1, ZV2, DZMV, VRT[3], VRTN[3];  
  Float_t RLPED[3], DLPED[3], RTPGAM[3], RDNPED, VRTSEC[3];
  Float_t FF, RDNPBS, RDINI, RDNSEC, YSTP;
  const Int_t NVY=(NXMAX+NXMAX*NXMAX)/2;
  Float_t VY[NVY];
  Double_t EPRIMV, PRIMV[3], PMISS[3], PPMISS, EMISS, MMISS2;
  Double_t EBM, PB, EINIT, PBM[3];
  Double_t EPAR, PPAR, ETOT, PDEC[3], PPDEC, EDEC, PJPAR;
  Double_t PDTP[3], ESECV, PSECV[3], PRES[3], ERES;

  ierr = 0;
  fNptall = Nptall;
  Lpart = Nptall - Ndecay;
  Nunme = 0;
  Iunmea = 0;
  fNpmeas = 1;
  for ( i=1; i<Lpart; i++ ) {
    if ( fKind[i]<0 || fKind[i]>1000 ) {
      Nunme++;
      Iunmea = i;
    }
    if ( fKind[i]<1000 ) fNpmeas++;
  }

  IfZfree = Ifzfree;
  if ( Ifzfree != 0 ) {
    if ( fLcstCB < 2 ) IfZfree = 0; 
    if ( fLcst < 3 ) IfZfree = 0; 
    if ( Ndecay==0 && fNpmeas<Lpart ) IfZfree = 0;
  }

  Float_t fMMiss = fMass[Iunmea];

  if ( Nunme>1 ) {
    printf("Number of unmeasured particles %d > 1\n", Nunme );
    return ierr = 1;
  }
  NY = fNpmeas*4;
  if ( Idecfl>0 ) NY += 4;
 
  NF = 1;
  if ( fNpmeas==Lpart )  NF += 3;
  if ( Idecfl>0 ) NF += 3;   
  NF += Ndecay;
  

//   filling the hypothesis on the decay chain
  for ( i=0; i<NPMAX; i++ )  Isecvt[i]=0;
  if ( Ndecay>0 ) {
    for ( j=0; j<Ndecay; j++ ) { 
      fMass[Lpart+j] = Amsdec[j];
      fKind[Lpart+j] = -Ikind[j] - 1000; // for the decay hypotheses
    }
  }

  lclsec = 0;  
  jdecfl = 0;
  if ( Idecfl>0 ) {
    jdecfl = Jdecay[Idecfl-1];
    Isecvt[ Lpart+Idecfl-1 ] = 1;
    idechk = Idecfl;      
 L111:
    for ( k=0; k<Ldecay[idechk-1]; k++ ) {
      ideca = Idecay[idechk-1][k];
      Isecvt[ideca] = 1;
    }
    while ( idechk>1 ) {
      idechk--;
      if ( Isecvt[ Lpart+idechk-1 ] > 0 )  goto L111;
    } 
 
    for ( i=1; i<=fLcst; i++ ) if ( Isecvt[i] !=.0 )  lclsec++;
  }
//  beam information and cluster information
  for ( i=0; i<=fLcst; i++ ) {
    j = i*4;
    for ( k=0; k<4; k++ )  Y[j+k]   = fPin[i][k];
  }
// index of primary and secondary vertex
  IV = 3;
  ISV = fNpmeas*4;

// search for the initial vertex Z coordinate when it is a free parameter
  WBEST[0]=WBEST[1]=0.;
  ZBEST[0]=ZBEST[1]=0.;

  ZBESTM = 0.;
  IFIRSL = 1;
  ZV1 = fPlim[IV][0];
  ZV2 = fPlim[IV][1];
  DZMV = 0.25;
  if ( IfZfree==0 ) NZDIV = 1;
  else NZDIV = (Int_t)((ZV2-ZV1)/DZMV) + 1;

//   calculating initial values for free parameters of the fit

  if ( Y[0] < 0.000001 )  Y[0] = 0.000001;
  EBM = 1./Y[0];
  PB = sqrt( EBM*EBM - fMass[0]*fMass[0] );
  EINIT = fMastg + EBM;   

  PBM[0] = 0.;
  PBM[1] = 0.;
  PBM[2] = PB;

  VRT[0] = Y[1];
  VRT[1] = Y[2];

  for ( IZPOS=0; IZPOS<NZDIV; IZPOS++ ) {
    if( NZDIV>1 )  Y[IV] = ZV1 + DZMV * (Float_t)IZPOS;
    VRT[2] = Y[IV];

    EPRIMV = 0.;
    for ( i=0; i<3; i++ )  PRIMV[i] = PDEC[i] = 0.;    
//   loop on clusters from a primary vertex
    for ( I=1; I<=fLcst; I++ ) {
      if ( ( Iunmea>0 && Isecvt[Iunmea]==0 ) || Isecvt[I]==0 ) {
        if ( I==Iunmea )  continue;
	IY = I*4;
	//RLPED[0] = fPin[I][3];
	RLPED[0] = Y[IY+3];
	RLPED[1] = Y[IY+1];
	RLPED[2] = Y[IY+2]; 
	Dvpxyz( fCalor[I], RLPED, VRT, DLPED );

	Pxyztpf( DLPED, RTPGAM ); 

	RDNPED = Vmodf(DLPED,3);
	if ( RDNPED < 0.00001 )  RDNPED = 0.00001;

	if ( Y[IY] < 0.000001 ) Y[IY] = 0.000001; 
	EPAR = 1./Y[IY] + fMass[I];
	PPAR = sqrt( EPAR*EPAR - fMass[I]*fMass[I] );
	EPRIMV += EPAR;
	for ( J=0; J<3; J++ ) {
	  PJPAR =  PPAR * DLPED[J]/RDNPED;
	  PRIMV[J] += PJPAR;
	  if ( Isecvt[I] != 0 )  PDEC[J] += PJPAR;
	}
      }
    }

//  the energy constraint calculation

    EMISS = 0.;
    if ( Idecfl>0 || Iunmea > 0 ) {
      Vsub( PBM, PRIMV, PMISS, 3 );
      PPMISS = Vmod( PMISS, 3 );
      if ( PPMISS < 0.000001 ) PPMISS = 0.000001;              
      if ( Isecvt[Iunmea] == 0 )
	EMISS = sqrt( PPMISS*PPMISS + fMMiss*fMMiss );
      else
	EMISS = sqrt( PPMISS*PPMISS + fMass[jdecfl]*fMass[jdecfl] );
    }
    F[0] = EINIT - EPRIMV - EMISS;

    FF = fabsf(F[0])+0.0000001;
    if ( WBEST[0] < 1./FF ) {
      WBEST[0] = 1./FF;
      ZBEST[0] = Y[IV];
      if ( Iunmea > 0 && fKind[Iunmea]<0 && Isecvt[Iunmea] == 0 )
                          Y[Iunmea*4] = 1./(EMISS-fMMiss);
      //if ( Iunmea > 0 && fKind[Iunmea]<0 && Isecvt[Iunmea] == 0 )
      //printf("1/Y[Iunmea*4] %lf %lf %lf\n",1./Y[Iunmea*4],
      //                            Y[Iunmea*4+1],Y[Iunmea*4+2] );
      if ( Idecfl > 0 ) {
	if ( Iunmea > 0 && Isecvt[Iunmea] == 0 )  Pxyztp( PDEC, PDTP );
	else Pxyztp( PMISS, PDTP );
        if ( PDTP[0] < 0.000001 ) PDTP[0] = 0.000001;
	Y[ISV] =
	  1./( sqrt( PDTP[0]*PDTP[0] + fMass[jdecfl]*fMass[jdecfl] ) - fMass[jdecfl]);
	Y[ISV+1] = PDTP[1];
	Y[ISV+2] = PDTP[2];
      }
    }
  }     //    end of the primary vertex loop

  Y[IV] = ZBEST[0];

  if ( Idecfl > 0 ) {              //  secondary vertex 

    for ( I=0; I<3; I++ )  PMISS[I] = Y[ISV+I];
    PMISS[0] += fMass[jdecfl];
    Ptpxyz( PMISS, PDEC );
    EDEC = sqrt( PMISS[0]*PMISS[0] + fMass[jdecfl]*fMass[jdecfl] );

    VRT[2] = Y[IV];
    for ( i=0; i<3; i++ ) VRTN[i] = -VRT[i];
    RDNPBS = 0.;
    IZPOS = 0;
    RDINI = -7.;

 L122:
    Y[ISV+3] = RDINI + DZMV * (Float_t)IZPOS;
    if ( fabsf( Y[ISV+3] ) < 0.00001 ) for ( i=0; i<3; i++ ) VRTSEC[i] = VRT[i];
    else {
      RLPED[0] = Y[ISV+3];
      RLPED[1] = Y[ISV+1];
      RLPED[2] = Y[ISV+2];
      Dvpxyz( 1, RLPED, VRTN, VRTSEC );
    } 
    RDNSEC = Vmodf( VRTSEC, 3 );
    if ( RDNSEC < 20.5 ) {
      ESECV = 0.;
      for ( i=0; i<3; i++ )  PSECV[i] = 0.;    
          
//   loop on clusters from the secondary vertex
      for ( I=1; I<=fLcst; I++ ) {
	if ( Isecvt[I] != 0 ) {
	  IY = I*4;
	  RLPED[0] = Y[IY+3];
	  RLPED[1] = Y[IY+1];
	  RLPED[2] = Y[IY+2]; 
	  Dvpxyz( fCalor[I], RLPED, VRTSEC, DLPED );

	  Pxyztpf( DLPED, RTPGAM ); 

	  RDNPED = Vmodf( DLPED, 3 );
	  if ( RDNPED < 0.00001 )  RDNPED = 0.00001;

	  if ( I==Iunmea && fKind[Iunmea]<0 ) {
	    Vsub( PDEC, PSECV, PMISS, 3 );
	    PPMISS = Vmod( PMISS, 3 );
	    if ( PPMISS < 0.000001 ) PPMISS = 0.000001;              
	    PPAR = PPMISS;
	    EPAR = sqrt( PPAR*PPAR + fMMiss*fMMiss );
	    EMISS = EPAR;
	  }
	  else {
	    if (Y[IY] < 0.000001) Y[IY] = 0.000001; 
	    EPAR = 1./Y[IY] + fMass[I];
	    PPAR = sqrt( EPAR*EPAR - fMass[I]*fMass[I] );
	  }
	  ESECV += EPAR;
	  for ( J=0; J<3; J++ ) {
	    PJPAR =  PPAR * DLPED[J]/RDNPED;
	    PSECV[J] += PJPAR;
	  }
	}
      }

//  the energy constraint calculation for particle decaying inflight 

      if ( Isecvt[Iunmea] != 0 && Iunmea==fNpmeas ) {
	Vsub( PDEC, PSECV, PMISS, 3 );
	PPMISS = Vmod( PMISS, 3 );
	if ( PPMISS < 0.000001 ) PPMISS = 0.000001;
	EMISS = sqrt( PPMISS*PPMISS + fMMiss*fMMiss );
	ETOT = ESECV + EMISS;            
	PPDEC = Vmod( PDEC, 3);
	MMISS2 = ETOT*ETOT - PPDEC*PPDEC;
      }
      else {
	PPDEC = Vmod( PSECV, 3 );
	MMISS2 = ESECV*ESECV - PPDEC*PPDEC;
      }

      F[1] = MMISS2 - fMass[jdecfl]*fMass[jdecfl];

      FF = fabsf( F[1] ) + 0.0000001;
      if ( WBEST[1] < 1./FF ) {
	WBEST[1] = 1./FF;
	ZBEST[1] = Y[ISV+3];
	RDNPBS = RDNSEC; 
	if ( Iunmea > 0 && fKind[Iunmea]<0 && Isecvt[Iunmea] != 0 )
	                   Y[Iunmea*4] = 1./(EMISS-fMMiss);
      }

      IZPOS++;
      goto L122;
    }
    else {
      if (Y[ISV+3] < 0.) {
	IZPOS++;
	goto L122;
      }          
    }           //  the secondary vertex radius is more than 23 cm

    Y[ISV+3] = ZBEST[1];

  }     //    end of the secondary vertex search

  //SIMNIT( NY, NF, 0.001 );
  SIMNIT( NY, NF, 0.01 );

//     The covariance matrix is set to zero

  for ( I=0; I<NVY; I++) VY[I] = 0.;
  for ( I=0; I<fNpmeas; I++) {
    K=I*4;
    fCov[0] = fErrSq[I][0];
    fCov[1] = 0.;
    fCov[2] = fErrSq[I][1];
    fCov[3] = 0.;
    fCov[4] = 0.;
    fCov[5] = fErrSq[I][2];
    fCov[6] = 0.;
    fCov[7] = 0.;
    fCov[8] = 0.;
    fCov[9] = fErrSq[I][3];
    if ( I==0 && IfZfree!=0 ) {
      fCov[9] = 0.;
      SIMSTP( IV, 1.8 );
    } 
    SMTOS( fCov, 1, VY, K+1, 4 );
    for ( j=K; j<K+4; j++ )  SIMLIM( j, fPlim[j][0], fPlim[j][1] );
  }

  if ( Iunmea > 0 && fKind[Iunmea]<0 ) {
    K = Iunmea*4;
    YSTP = 0.06*Y[K];   
    SIMSTP( K, YSTP ); 
  }

  if ( Idecfl > 0 ) {
    SIMSTP( ISV, 0.06*Y[ISV] ); 
    SIMLIM( ISV, 0., 200. );
    SIMSTP( ISV+1, 0.04 ); 
    SIMLIM( ISV+1, 0., 0. );
    SIMSTP( ISV+2, 0.05 ); 
    SIMLIM( ISV+2, 0., 0. );
    SIMSTP( ISV+3, 2.5 );
    SIMLIM( ISV+3, -15., 40. );
  }
   
//     MAIN FITTING LOOP
 L11:
//     The values of the constraint equations F(J) are computed
//     using the corrected values Y(I)

  if ( Y[0] < 0.000001 )  Y[0] = 0.000001;
  EBM = 1./Y[0];
  PB = sqrt( EBM*EBM - fMass[0]*fMass[0] );
  EINIT = fMastg + EBM;   
  
  PBM[0] = 0.;
  PBM[1] = 0.;
  PBM[2] = PB;

  VRT[0] = Y[1];
  VRT[1] = Y[2];
  VRT[2] = Y[IV];

  if ( Idecfl > 0 ) {
    for ( i=0; i<3; i++ ) VRTN[i] = -VRT[i];
    RLPED[0] = Y[ISV+3];
    RLPED[1] = Y[ISV+1];
    RLPED[2] = Y[ISV+2];
    Dvpxyz( 1, RLPED, VRTN, VRTSEC ); 
    RDNSEC = Vmodf( VRTSEC, 3 );
  }

  EPRIMV = 0.;
  for ( i=0; i<3; i++ )  PRIMV[i] = 0.;    
//   loop on clusters from a primary vertex
  if ( lclsec < fLcst ) {  
    for ( I=1; I<=fLcst; I++ ) {
      if ( Isecvt[I]==0 ) {
	IY = I*4;
	RLPED[0] = Y[IY+3];
	RLPED[1] = Y[IY+1];
	RLPED[2] = Y[IY+2]; 
	Dvpxyz( fCalor[I], RLPED, VRT, DLPED );

	Pxyztpf( DLPED, RTPGAM ); 

	fProut[I][1] = RTPGAM[1];
	fProut[I][2] = RTPGAM[2];

	RDNPED = Vmodf(DLPED,3);
	if ( RDNPED < 0.00001 )  RDNPED = 0.00001;

	if ( Y[IY] < 0.000001 ) Y[IY] = 0.000001; 
	EPAR = 1./Y[IY] + fMass[I];
	PPAR = sqrt( EPAR*EPAR - fMass[I]*fMass[I] );
	EPRIMV += EPAR;
	for ( J=0; J<3; J++ ) {
	  PJPAR =  PPAR * DLPED[J]/RDNPED;
	  PRIMV[J] += PJPAR;
	}
      }
    }
  }

  if ( Idecfl > 0 ) {
    if ( Y[ISV] < 0.000001 )  Y[ISV] = 0.000001;
    EDEC = 1./Y[ISV] + fMass[jdecfl];
    EPRIMV += EDEC;
    PDTP[0] = sqrt( EDEC*EDEC - fMass[jdecfl]*fMass[jdecfl] );
    PDTP[1] = Y[ISV+1];
    PDTP[2] = Y[ISV+2];
    for ( j=0; j<3; j++ ) fProut[jdecfl][j] = Y[ISV+j];
    Ptpxyz( PDTP, PDEC );
    for ( J=0; J<3; J++ ) PRIMV[J] += PDEC[J];
  }

  Vsub( PBM, PRIMV, PMISS, 3 );

  if ( Isecvt[Iunmea] == 0 && Iunmea==fNpmeas ) {
//  the energy constraint calculation in case of fully unmeasured particle 
    Pxyztp( PMISS, PDTP );
    if ( PDTP[0] < 0.000001 ) PDTP[0] = 0.000001;
    EMISS = sqrt( PDTP[0]*PDTP[0] + fMMiss*fMMiss );
    fProut[Iunmea][0] = EMISS - fMMiss;
    fProut[Iunmea][1] = PDTP[1];
    fProut[Iunmea][2] = PDTP[2];
    F[0] = EINIT - EPRIMV - EMISS;
    NFF = 1;
  }
  else {
// the energy and 3-momentum constraints calculation in case of all measured particles
    F[0] = EINIT - EPRIMV;
    for ( J=0; J<3; J++ ) F[J+1] = PMISS[J];
    NFF = 4;
    //printf("F %lf %lf %lf %lf\n",F[0],F[1],F[2],F[3]);
  }

  if ( Idecfl > 0 ) {

    ESECV = 0.;
    for ( i=0; i<3; i++ )  PSECV[i] = 0.;    
          
//   loop on clusters from the secondary vertex
    for ( I=1; I<=fLcst; I++ ) {
      if ( Isecvt[I] != 0 ) {
	IY = I*4;
	RLPED[0] = Y[IY+3];
	RLPED[1] = Y[IY+1];
	RLPED[2] = Y[IY+2]; 
	Dvpxyz( fCalor[I], RLPED, VRTSEC, DLPED );

	Pxyztpf( DLPED, RTPGAM ); 

	RDNPED = Vmodf(DLPED,3);
	if ( RDNPED < 0.00001 )  RDNPED = 0.00001;

	if (Y[IY] < 0.000001) Y[IY] = 0.000001; 
	EPAR = 1./Y[IY] + fMass[I];
	ESECV += EPAR;
	PPAR = sqrt( EPAR*EPAR - fMass[I]*fMass[I] );
	for ( J=0; J<3; J++ ) {
	  PJPAR =  PPAR * DLPED[J]/RDNPED;
	  PSECV[J] += PJPAR;
	}
      }
    }

    Vsub( PDEC, PSECV, PMISS, 3 );
    if ( Isecvt[Iunmea] != 0 && Iunmea==fNpmeas ) {
//  the energy constraint calculation in case of unmeasured particle 
      Pxyztp( PMISS, PDTP );
      if ( PDTP[0] < 0.000001 ) PDTP[0] = 0.000001;
      EMISS = sqrt( PDTP[0]*PDTP[0] + fMMiss*fMMiss );
      fProut[Iunmea][0] = EMISS - fMMiss;
      fProut[Iunmea][1] = PDTP[1];
      fProut[Iunmea][2] = PDTP[2];
      F[NFF] = EDEC - ESECV - EMISS;
      NFF++;
    }
    else {
// the energy and 3-momentum constraints calculation in case of measured hadron
      F[NFF] = EDEC - ESECV;
      NFF++;
      for ( J=0; J<3; J++ ) F[NFF+J] = PMISS[J];
      NFF += 3;          
    }
  }

  ICDEC = 0;
// calculation of the constraints due to decaying particles
//  (1 more constraint for every decaying particle)
  IFSCND = 0;
  if ( Idecfl > 0 ) IFSCND = 1;
  if ( Ndecay-IFSCND > 0 ) {
    for ( I=Lpart; I<Nptall; I++ ) {
      if ( I != jdecfl ) {
	K = I - Lpart;
	LLDEC = Ldecay[K];
	ICDEC++;
	ERES = 0.;
	for ( J=0; J<3; J++ ) PRES[J]=0.;
	for ( IDE=0; IDE<LLDEC; IDE++ ) { 
	  INDE = Idecay[K][IDE];
	  I2 = INDE * 4;
	  if ( INDE < fNpmeas ) {
	    if ( Y[I2] < 0.000001 )  Y[I2] = 0.000001;
	    EDEC = 1./Y[I2] + fMass[INDE];
	  }
          else {
	    if ( fProut[INDE][0] < 0.000001 ) fProut[INDE][0] = 0.000001;
	    EDEC = fProut[INDE][0] + fMass[INDE];
	  }
	  PDTP[0] = sqrt( EDEC*EDEC - fMass[INDE]*fMass[INDE] );
	  PDTP[1] = fProut[INDE][1];
	  PDTP[2] = fProut[INDE][2];
	  Ptpxyz( PDTP, PDEC );
	  for ( J=0; J<3; J++ ) PRES[J] += PDEC[J];
	  ERES += EDEC;
	}

	Pxyztp( PRES, PDTP ); 
	if ( PDTP[0] < 0.000001 ) PDTP[0] = 0.000001;
	fProut[I][0] = sqrt( PDTP[0]*PDTP[0] + fMass[I]*fMass[I] ) - fMass[I];
	fProut[I][1] = PDTP[1];
	fProut[I][2] = PDTP[2];
	MMISS2 = ERES*ERES - PDTP[0]*PDTP[0];

        if (MMISS2<0.000001) F[NFF]  = fMass[I]*fMass[I] - MMISS2;
	else F[NFF]  = fMass[I] - sqrt(MMISS2);
        //F[NFF]  = fMass[I]*fMass[I] - MMISS2;
	NFF++;
      }
    }
  }

  IEF = APLCON( Y, VY, F);

  if ( IEF < 0 ) goto L11;
     
  if ( IEF > 0 ) fChisq = fabsf(CHSQ) + 10000.*IEF + 100000.;  //     Fit is unsuccessful
  else {
    // fit is successful
    fProut[0][0] = 1./Y[0]; // beam energy
    fProut[0][1] = Y[1];  //X vertex
    fProut[0][2] = Y[2];  //Y
    fProut[0][3] = Y[3]; //Z

    for ( I=1; I<=fLcst; I++ ) {
      I2 = I*4;
      fProut[I][0] = 1./Y[I2];
      fProut[I][3] = Y[I2+1];
      fProut[I][4] = Y[I2+2];
    }

    for (i=0; i<NY; i++)  fPulls[i] = A[i];

    fChisq = (Double_t)CHSQ;

    if ( Idecfl>0 ) fDeclen = (Double_t)Y[ISV+3]; 
  }
  fNDF = ND; 
  if ( fChisq>=100000. )  return ierr;

  Double_t Eb = (Double_t)fProut[0][0];
  Double_t Pb = sqrt( Eb*Eb - fMass[0]*fMass[0] );
  fBeam.SetPxPyPzE( 0., 0., Pb, Eb );
  Trevcm();
  fBeamcm = Lorencm( fBeam );
  fVertex.SetXYZ( fProut[0][1], fProut[0][2], fProut[0][3] );
  Double_t ph3[3];
  Double_t vm[3];
  for (j=1; j<Nptall; j++) {
    for (i=0; i<3; i++) vm[i] = fProut[j][i];
    vm[0] += fMass[j]; 
    vm[0] = sqrt( vm[0]*vm[0] - fMass[j]*fMass[j] ); 
    Ptpxyz( vm, ph3 );
    fParticles[j-1].SetXYZM( ph3[0], ph3[1], ph3[2], fMass[j] );
    fParticlescm[j-1] = Lorencm( fParticles[j-1] );
  }
  return ierr;
}
   

//---------------------------------------------------------------------------
void TA2KFit::Dvpxyz( Int_t Ical, Float_t RLPED[3], Float_t VERT[3], Float_t DVPED[3] )
{
  Int_t i;
  if ( Ical == 1 ) {
    DVPED[0] = RLPED[0] * sinf(RLPED[1]) * cosf(RLPED[2]);
    DVPED[1] = RLPED[0] * sinf(RLPED[1]) * sinf(RLPED[2]);
    DVPED[2] = RLPED[0] * cosf(RLPED[1]);
  }
  else if ( Ical == 2 ) {
    DVPED[0] = RLPED[1] * cosf(RLPED[2]);
    DVPED[1] = RLPED[1] * sinf(RLPED[2]);
    DVPED[2] = sqrtf( RLPED[0]*RLPED[0] - RLPED[1]*RLPED[1] ); 
  }
  for ( i=0; i<3; i++ )  DVPED[i] -= VERT[i];
}

//---------------------------------------------------------------------------
void TA2KFit::Trevcm()
{
  Double_t Pbm = fBeam.P();
  Double_t Etot = fBeam.E() + fMastg;
  Ecm = sqrt( Etot*Etot - Pbm*Pbm);
  Betacm = Pbm/Etot;
  Velcm[0] = fBeam.Px()/Etot;
  Velcm[1] = fBeam.Py()/Etot;
  Velcm[2] = fBeam.Pz()/Etot;
  Velcm[3] = 1./sqrt( 1. - Betacm*Betacm );  
}

//---------------------------------------------------------------------------
void TA2KFit::Trevcmf( TLorentzVector Beam, Double_t Mastg, Double_t Vel[4] )
{
  Double_t Pbm = Beam.P();
  Double_t Etot = Beam.E() + Mastg;
  Double_t Ecmf = sqrt( Etot*Etot - Pbm*Pbm);
  Double_t Beta = Pbm/Etot;
  Vel[0] = Beam.Px()/Etot;
  Vel[1] = Beam.Py()/Etot;
  Vel[2] = Beam.Pz()/Etot;
  Vel[3] = 1./sqrt( 1. - Beta*Beta );  
}

//---------------------------------------------------------------------------
void TA2KFit::Trevrest( TLorentzVector P4part, Double_t Vel[4] )
{
  Double_t Ppart = P4part.P();
  Double_t Epart = P4part.E();
  Double_t Beta = Ppart/Epart;
  Vel[0] = P4part.Px()/Epart;
  Vel[1] = P4part.Py()/Epart;
  Vel[2] = P4part.Pz()/Epart;
  Vel[3] = 1./sqrt( 1. - Beta*Beta );  
}

//---------------------------------------------------------------------------
Double_t TA2KFit::Zeta()
{
  Double_t zeta = 0., edif;
  Double_t veleta[4];
  const Double_t etamass = 0.5475, pi0mass=0.1349764;
  Double_t zfact = 6./(etamass-3.*pi0mass)/(etamass-3.*pi0mass);
  TLorentzVector p4eta, p4pi0, p4pi0eta;
  p4eta = Particle( 17, 1);
  Trevrest ( p4eta, veleta );
  for ( Int_t ip = 1; ip<=3; ip++ ) {
    p4pi0 = Particle( 7, ip);
    p4pi0eta = Loren( veleta, p4pi0 );
    edif = p4pi0eta.E() - etamass/3.;
    zeta += zfact * edif*edif;
  }
  return zeta;
}

//---------------------------------------------------------------------------
TLorentzVector TA2KFit::Loren( Double_t Vel[4], TLorentzVector LVlab )
{
  TLorentzVector LVvel;
  Double_t vcm[4];
  Double_t betpa = Vel[0] * LVlab.Px()
                 + Vel[1] * LVlab.Py() 
                 + Vel[2] * LVlab.Pz();
  Double_t bpgam = ( betpa * Vel[3]/(Vel[3]+1.) - LVlab.E() ) * Vel[3];
  vcm[0] = LVlab.Px() + bpgam * Vel[0];
  vcm[1] = LVlab.Py() + bpgam * Vel[1];
  vcm[2] = LVlab.Pz() + bpgam * Vel[2];
  vcm[3] = ( LVlab.E() - betpa ) * Vel[3];
  LVvel.SetPxPyPzE( vcm[0], vcm[1], vcm[2], vcm[3] );
  return LVvel;
}

//---------------------------------------------------------------------------
TLorentzVector TA2KFit::Lorencm( TLorentzVector LVlab )
{
  TLorentzVector LVcm;
  Double_t vcm[4];
  Double_t betpa = Velcm[0] * LVlab.Px()
                 + Velcm[1] * LVlab.Py() 
                 + Velcm[2] * LVlab.Pz();
  Double_t bpgam = ( betpa * Velcm[3]/(Velcm[3]+1.) - LVlab.E() ) * Velcm[3];
  vcm[0] = LVlab.Px() + bpgam * Velcm[0];
  vcm[1] = LVlab.Py() + bpgam * Velcm[1];
  vcm[2] = LVlab.Pz() + bpgam * Velcm[2];
  vcm[3] = ( LVlab.E() - betpa ) * Velcm[3];
  LVcm.SetPxPyPzE( vcm[0], vcm[1], vcm[2], vcm[3] );
  return LVcm;
}

//---------------------------------------------------------------------------
TLorentzVector TA2KFit::NewLVec( TVector3 Vcl, Double_t Ekin,
                                                        Double_t Pmass )
{
  TLorentzVector LVnew;
  Double_t vm[3], ph3[3];
  Double_t Etot = Ekin + Pmass;
  vm[0] = sqrt(Etot*Etot - Pmass*Pmass);
  vm[1] = Vcl.Theta(); 
  vm[2] = Vcl.Phi(); 
  Ptpxyz( vm, ph3 );
  LVnew.SetXYZM( ph3[0], ph3[1], ph3[2], Pmass );
  return LVnew;
}

//---------------------------------------------------------------------------
void TA2KFit::SIMNIT(Int_t NI, Int_t NJ, Float_t EPSIL)
{
/************************************************************************
     initialize and define dimension  NI of X, NJ of Y/F, and
                    define debug flag JDEBUG, and epsilon (F)
************************************************************************/

  Int_t I;

  NX    =NI;
  MYF   =NJ;
  NUM   =0;
  IFLG  =0;
  INIT  =0;
  ND    =0;
  CHSQ  =0.;
  EPSF  =EPSIL;
  ISTAT =0;
  ITER = 0;
  IDR  = 0;
  XD[0] = XD[1] = 0.;
  //     ... limits XL and steps ST
  for ( I=0; I<NXMAX; I++ ) XL[I][0]=XL[I][1]=FC[I]=ST[I]=0.;  
  for ( I=0; I<NXMAX+NFMAX; I++ ) H[I]=0.;  

  //     clear derivative matrix A and ...
  for ( I=0; I<1000; I++ ) A[I]=0.;
  for ( I=0; I<100; I++ ) DR[I][0]=DR[I][1]=0.;
}

//---------------------------------------------------------------------------
void TA2KFit::SMTOS(Float_t* V, Int_t I, Float_t* W, Int_t J, Int_t N)
{
/************************************************************************
     Copy symmetric N-by-N matrix or a N-by-N submatrix of a  symmetric
     matrix to another symmetric matrix
     N rows and columns of the matrix V, starting from diagonal element
     (i,i), are copied to the matrix W, starting  in  diagonal  element
     (j,j). thus if a complete symmetric matrix has to be copied, i=j=1
     has to be used.
************************************************************************/

  Int_t IM, JM, K, L;

  IM=(I*I+I)/2-1;
  JM=(J*J+J)/2-1;
  for ( K=0; K<N; K++ ) {
    for ( L=0; L<=K; L++ ) {
      W[JM]=V[IM];
      IM++;
      JM++;
    }
    IM+=I-1;
    JM+=J-1;
  }
}

//---------------------------------------------------------------------------
void TA2KFit::SIMLIM(Int_t IA, Float_t XLOW, Float_t XHIG)
{
/************************************************************************
    define limits for variable X(IA)
************************************************************************/

  if ( IA<0 || IA>=NX ) return;
  //     lower or upper limit of X(IA)
  XL[IA][0] = ( XLOW<XHIG )  ? XLOW : XHIG;
  XL[IA][1] = ( XLOW>=XHIG ) ? XLOW : XHIG;
}

//---------------------------------------------------------------------------
void TA2KFit::SIMSTP(Int_t IA, Float_t STEP)
{
/************************************************************************
     define step size for num.dif. of variable X(IA)
************************************************************************/

  if ( IA<0 || IA>=NX ) return;
  //     step for numerical differentiation of X(IA)
  ST[IA]=fabsf(STEP);
}

//---------------------------------------------------------------------------
Int_t TA2KFit::APLCON(Float_t* X, Float_t* VX, Float_t* F)
{
/************************************************************************
    Apply constraints F(J) = function of X   J=1,NF
    to data X(1)...X(nx) with covariance matrix VX.
    VX(.) in symmetric storage mode
       1,1  1,2  2,2  1,3  2,3  3,3 ...

    Restriction because of limited dimensions of
    internal arrays:
       NX = 50  is maximum
       NF = 20  is maximum

    Usage
    =====

            SIMNIT(NX,NF,JDEBUG,EPSILON)

                   JDEBUG = 0   no printout
                          > 0   more and more printout
                   EPSILON=     precision required for constraints
 
       Now the variables X(1) ... X(NX) have to be defined
       and their covariance matrix V(1) ...
       variables may be measured ones or unmeasured ones.
       unmeasured variables are characterized by zero elements
       of the covariance matrix (at least the corresponding
       diagonal element has to be zero). For unmeasured variables
       the value of x has to be some reasonable initial value.
       in addition it is necessary to define some step size
       for numerical differentiation for the unmeasured
       variables (for measured values the necessary step size
       is taken from the standard deviation in the covar.
       matrix). The call to define step ST for variable X(i) is
               SIMSTP(I,ST)
       The user may optionally define limits for physical regions
       of variable x(i) by
               SIMLIM(I,XLOW,XHIG)
       The following coding example represents the loop, which
       performs the constrained fit. during the loop the
       variables X(1)...X(NX) are modified, to perform the
       numerical differentiation and to apply corrections
       during the fit. The user has to supply the code to
       calculate the constraint equations F(1)...F(NF).
       Finally (at convergence) the covariance matrix is modified
       to the matrix for the fitted values, which (for sufficient
       constraints) has nonzero elements for the previously
       unmeasured variables.

    10 F(1)=function of X(1) ... X(NX)
       ...
       F(NF) = function of X(1) ... X(NX)
            APLCON(X,VX,F,IRET)
       IF(IRET.LT.0) GOTO 10
            SIMNCH(NDEG,CHISQ)      optional to get chisquare

       IRET = 0   convergence reached
       IRET = 1   bad convergence         no convergence
       IRET = 2   too many iterations         - " -
       IRET = 3   unphysical region           - " -
       IRET = 4   ND less or equal zero       - " -

       Non-convergence is assumed for
       chisquare above 4 standard deviations for first 10 iterations
       chisquare above 3 standard deviations for next 10 iterations
       more than 20 iterations

       The user may scale up his covariance matrix to force more
       cases with convergence.

       Remark to precision of the constraints:
       a necessary condition for convergence is the reduction of
       the constraints to about EPSILON. 
       Convergence may be
       difficult due to roundoff-errors, if the required
       accuracy is too high.

       Remark to differentiation:
       By default numerical differentiation is done, to calculate
       the nx*nf elements of the derivative matrix for the
       constraits w.r.t the variables. Usually this works well.
       The user may calculate the elements in his program, for
       reasons of speed or in cases, where the num. diff. fails.
       The derivative matrix is array A
       with the definition
       df(j)/dx(i) = A(i+nx*(j-1))
       and has to be defined by the user before the call of
       APLCON. In the first iteration of the first case the
       program will automatically compare the elements with
       values calculated numerically and will printout elements
       with disagreement.

************************************************************************/

  Int_t I, II, J, IJ, IA, K, NM, M, JK, IK, NRANK, JRET;
  static Int_t MXF, IRET, IUNPH;
  static Int_t ICNT=0, NXF=0, NCST=0;
  static const Int_t NDIMR = NXMAX + NFMAX;
  static const Int_t NDIMW = (NDIMR+NDIMR*NDIMR)/2;
  static Float_t XS[NXMAX], DX[NXMAX], DXP[NXMAX], R[NDIMR], W[NDIMW];
  static Int_t NRD[NXMAX];
  Float_t SWII, SUM, WJK, FTEST;
  static Float_t FTESTP=0.;
  //static Double_t CHSQP=0.;
  static Float_t CHSQP=0.;

  if(ISTAT!=0) goto L40;
  IFLG=ICNT;
  ICNT++;
//     initialization----------------------------------------------------
//     define loop parameters
  NXF=NX+MYF;
  MXF=(NXF*NXF+NXF)/2;

  SIMMAT(VX);
//     count nr of degrees of freedom
  ND=MYF;
  II=0;
  for ( I=1; I<=NX; I++ ) {
    II+=I;
//     save initial X values and reset correction DX
    XS[I-1]=X[I-1];
    DX[I-1]=0.;
//     check unphysical region
    if ( XL[I-1][0] != XL[I-1][1] ) {
      if ( X[I-1]<XL[I-1][0] || X[I-1]>XL[I-1][1] ) {
	IRET=3;
	//printf("unph %d %lf %lf %lf\n",I,X[I-1],XL[I-1][0],XL[I-1][1]);
	goto L99;
      }
    }
    if ( VX[II-1] <= 0.) {
//        unmeasured variable
      ND--;
//        clear remaining elements
      IJ=II-I;
      for ( J=1; J<=NX; J++ ) {
        if( J<=I ) IJ++;
	VX[IJ-1]=0.;
        if( J>=I ) IJ+=J;
      }
    }
  }

  if ( ND <= 0 ) {
//        insufficient information
    IRET=4;
    goto L99;
  }
//     initial value of ftest
  FTESTP=0.;
  for ( J=0; J<MYF; J++ )  FTESTP += fabsf(F[J]);
  FTESTP /= (Float_t)MYF;
  ITER=0;
  NCST=0;
  CHSQ=0.;

//     prepare next iteration--------------------------------------------

 L30:
  ISTAT=1;
  IRET=-1;
//     define right hand side of equation R
  for ( I=0; I<NX; I++ )  R[I]=0.;

  for ( J=0; J<MYF; J++ ) {
//     fc is used in SIMDER
    FC[J]=F[J];
    R[NX+J]=-F[J];
  }

  if ( ITER != 0 ) {
//        define steps = 0.5 sigma from W
    II=0;
    for ( I=1; I<=NX; I++ ) {
      II+=I;
      if ( ST[I-1]==0. || W[II-1]==0.)  continue;
      SWII = 0.5*sqrtf(fabsf(W[II-1]));
      if ( SWII < ST[I-1] )  ST[I-1] = SWII;
    }
  }

//     loop for numerical calculation of derivatives---------------------

 L40:
  if (ISTAT != 1) goto L90;
  JRET = SIMDER(X,F);
  if (JRET < 0) return IRET;

//     construct matrix w and correct vector r---------------------------

//     insert -v and a into w and update r
  IJ=(NX*NX+NX)/2;
  for ( I=0; I<IJ; I++ ) W[I]=-VX[I]; 

  IA=0;
  for ( J=1; J<=MYF; J++ ) {
    for ( I=1; I<=NX; I++ ) {
      R[NX+J-1]+=A[IA+I-1]*DX[I-1];
      W[IJ+I-1]=A[IA+I-1];
    }
    IJ+=NX;
    for ( K=1; K<=J; K++ ) {
      IJ++;
      W[IJ-1]=0.;
    }
    IA+=NX;
  }

//     calculate step delx-----------------------------------------------

  ITER++;

//     First part of matrix inversion, making use of
//     the fact, that all elements corresponding to
//     measured variables are already the inverse elements

  NM=0;
  II=0;
  for ( I=1; I<=NX; I++ ) {
    II+=I;
    if ( W[II-1] < 0.) {
      DR[I-1][0]=0.;
      NM++;
      NRD[NM-1]=I;
    }
    else {
      W[II-1]=0.;
      DR[I-1][0]=1.;
    }
  }

  for ( I=NX+1; I<=NXF; I++ ) {
    DR[I-1][0]=1.;
    II=(I*I-I)/2;
    for ( M=1; M<=NM; M++) {
      J=NRD[M-1];
      SUM=0.;
      JK=(J*J-J)/2;
      for ( K=1; K<=NX; K++ ) {
	if (K<=J) JK++;
	if (DR[K-1][0] == 0.) SUM+=W[JK-1]*W[II+K-1];
	if (K>=J) JK+=K;
      }
      H[J-1]=SUM;
    }
    for ( K=I; K<=NXF; K++ ) {
      IK=(K*K-K)/2;
      WJK=0.;
      for ( M=1; M<=NM; M++) {
	J=NRD[M-1];
	WJK+=W[IK+J-1]*H[J-1];
      }
      W[IK+I-1]+=WJK;
    }
    for ( M=1; M<=NM; M++) {
      J=NRD[M-1];
      W[II+J-1]=-H[J-1];
    }
  }
//     save right hand side for Chi**2 calculation
  for ( J=0; J<MYF; J++) H[J]=R[NX+J];
    
//     complete matrix inversion and calculate chisquare

  NRANK = SMINVv(W,R,-NXF,1);
//     Chi**2 calculation
  CHSQP=CHSQ;
  CHSQ=0.;
  for ( J=0; J<MYF; J++) CHSQ-=H[J]*R[NX+J];

  for ( I=0; I<NX; I++) {
    DXP[I]=DX[I];
    DX[I]=R[I];
  }
  goto L84;
//     make cutstep
 L80:
  for ( I=0; I<NX; I++) DX[I]=0.5*(DX[I]+DXP[I]);
//     correct x and return to test constraints
 L84:
  for ( I=0; I<NX; I++) X[I]=XS[I]+DX[I];
  ISTAT=2;
//     check unphysical region
  IUNPH=0;
  for ( I=0; I<NX; I++) {
    if (XL[I][0] != XL[I][1]) {
      if (X[I] < XL[I][0] || X[I] > XL[I][1]) IUNPH=1;
    }
  }
  if (IUNPH != 0) goto L95;
  IRET=-2;
  return IRET;

//     convergence check-------------------------------------------------

//     calculate value of constraints and compare
 L90:
  FTEST=0.;
  for ( I=0; I<MYF; I++ )  FTEST+=fabsf(F[I]);
  FTEST /= (Float_t)MYF;
  if ( FTEST < FTESTP ) FTESTP=FTEST; 

  IUNPH=0;
  for ( I=0; I<NX; I++ ) {
    if ( XL[I][0] != XL[I][1] ) {
      if ( X[I]<XL[I][0] || X[I]>XL[I][1] )  IUNPH=1;
    }
  }

//     divergence/convergence tests
 L95:
  if ( IUNPH != 0 || FTEST > 1.1*FTESTP+EPSF ) {
//        divergence, make cut steps
    NCST++;
    if ( NCST<5 ) goto L80;
  }
  else if ( NCST==0 ) {
    if ( (ITER>=2 || CHSQ<CHLIM(1,ND)) &&
         (FTEST<EPSF && CHSQ-CHSQP<0.1) ) {
//           convergence

//           pulls
      II=0;
      for ( I=1; I<=NX; I++ ) {
	II+=I;
	A[I-1]=0.;
	if (VX[II-1] > 0.) {
	  if (VX[II-1]-W[II-1] > 0.) A[I-1]=DX[I-1]/sqrtf(VX[II-1]-W[II-1]);
        }
      }
      for ( I=0; I<(NX*NX+NX)/2; I++ )  VX[I]=W[I];

      IRET =0;
      ISTAT=0;
      return IRET;
    }
  }
//     continue with iteration or stop
  NCST=0;
  //if (ITER >= 2 && CHSQ > CHLIM(4,ND) )  IRET=1;
  if (ITER >= 3 && CHSQ > CHLIM(4,ND) )  IRET=1;
  if (ITER > 10 && CHSQ > CHLIM(3,ND) )  IRET=1;
  if (ITER > 20                       )  IRET=2;
  if ( IRET<0 ) goto L30;
 L99:
  ISTAT=0;
  return IRET;
    
}

//---------------------------------------------------------------------------
void TA2KFit::SIMMAT(Float_t* VX)
{
/************************************************************************
    check derivative matrix, set flag NUM ...
       and define steps for num. diff. from covariance matrix
       NUM = 0      analytical derivatives
       MUN = 1, 2   numerical derivatives
                    (=2 for comparison analytica/numerical der.)
************************************************************************/

  Int_t IJ, II, I;
  Float_t VII, SVII;

  NUM=0;
  for ( IJ=0; IJ<NX*MYF; IJ++ )  if( A[IJ] != 0.) goto L20;
//     set NUM flag for zero derivative matrix
  NUM=1;
//     force numerical derivatives for first call
 L20:
  //if (NUM==0 && IFLG==0 && IDEBUG>=0) NUM=2;
  if (NUM==0) return;
//     define steps from covariance matrix for numer. differentiation
  II=0;
  for ( I=1; I<=NX; I++ ) {
    II+=I;
    VII=fabsf(VX[II-1]);
    if ( VII != 0. ) {
      SVII = 0.5*sqrtf(VII);
      if ( ST[I-1] != 0.) {
        if ( SVII < ST[I] )  ST[I-1] = SVII;
      }
      else   ST[I-1] = SVII; 
    }
  }
}

//---------------------------------------------------------------------------
Int_t TA2KFit::SIMDER(Float_t* X, Float_t* F)
{
/************************************************************************
     derivative calculation
************************************************************************/

  static Int_t ILR=0;
  static Int_t IREDUC;
  Int_t J, IJ, JRET;
  Float_t DER;
  static Float_t XSAVE=0.;
  Bool_t LIMDEF;

//     initialize derivative calculation---------------------------------
  if (INIT==0) {
//        start with first variable
    IDR=0;
    INIT=1;
    goto L40;
  }

  JRET=-1;
//     derivative calculation -------------------------------------------
  if (IDR<0) {
//        another step required, save constraint values ...
    for ( J=0; J<MYF; J++ )  H[J]=F[J];
//        ... and set next step
    IDR=-IDR;
    X[IDR-1]=XD[1];
    return JRET;
  }

  X[IDR-1]=XSAVE;
  IJ=IDR;
  for ( J=0; J<MYF; J++ ) {
//     calculation of numerical derivative
    if (ILR==0)  DER=0.5*(H[J]-F[J])/ST[IDR-1]; // symmetric formula
    else {    //  asymmetric formula
      DER=0.5*(3.0*FC[J]+F[J]-4.0*H[J])/ST[IDR-1];
      if (ILR==2) DER=-DER;
    }
//     insert into A
    A[IJ-1]=DER;
    IJ+=NX;
  }
//     test end condition
 L30:
  if (IDR==NX) {
    JRET=INIT=0;
    if (NUM==2) NUM=0;
    return JRET;
  }
//     next variable
 L40:
  JRET=-1;
  IDR++;
  if ( ST[IDR-1] == 0.) goto L30;
  IREDUC=0;
  XSAVE=X[IDR-1];
//     central differentiation
 L50:
  ILR=0;
  LIMDEF = ( XL[IDR-1][0] != XL[IDR-1][1] );
  XD[0]=XSAVE+ST[IDR-1];
  if ( LIMDEF && ( XD[0] > XL[IDR-1][1] ) ) {
//        above upper limit
    XD[0]=XSAVE-ST[IDR-1];
    if ( LIMDEF && ( XD[0] < XL[IDR-1][0] ) ) goto L90;
    XD[1]=XSAVE-ST[IDR-1]-ST[IDR-1];
    if ( LIMDEF && ( XD[1] < XL[IDR-1][0] ) ) goto L90;
//        left differentiation
    ILR=1;
  }
  else {
    XD[1]=XSAVE-ST[IDR-1];
    if ( LIMDEF && ( XD[1] < XL[IDR-1][0] ) ) {
//           below lower limit
      XD[1]=XSAVE+ST[IDR-1]+ST[IDR-1];
      if ( LIMDEF && ( XD[1] > XL[IDR-1][1] ) ) goto L90;
//           right differentiation
      ILR=2;
    }
  }
//     first step
  X[IDR-1]=XD[0];
  IDR=-IDR;
  return JRET;
//     reduce step size
 L90:
  if ( IREDUC >= 4 ) goto L30;
  ST[IDR-1]/=3.0;
  IREDUC++;
  goto L50;
}

//---------------------------------------------------------------------------
Int_t TA2KFit::SMINVv(Float_t* V, Float_t* B, Int_t NARG, Int_t M)
{
/************************************************************************
    Obtain solution of a system of linear equations V *  X  =  B  with
    symmetric matrix V and inverse (for M =  1)  or  matrix  inversion
    only (for M = 0).

                  - - - -
       CALL SMINVv(V,B,N,M,NRANK)
                  - -     -----

          V = symmetric N-by-N matrix in symmetric storage mode
              V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
              replaced by inverse matrix
          B = N-vector   (for M = 0 use a dummy argument)
              replaced by solution vector
          M = see above


    Method of solution is by elimination selecting the  pivot  on  the
    diagonal each stage. The rank of the matrix is returned in  NRANK.
    For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
    matrix V and the corresponding elements of  B  are  set  to  zero.
    SMINVv can be used for a dimension up to 100 (see INVCDR).
************************************************************************/

  Int_t N, NI, JK, JL, LK, NRANK;
  Int_t IJ, II, I, K, J, L, JJ, KK, III;
  Float_t VKK, D, E;
  Double_t vdbl;
  const Float_t EPS=1.e-6;
  //const Float_t EPS2=1.e-16;
  const Float_t EPS2=1.e-28;

//     construct table

  N=abs(NARG);
  if (NARG>0) {
//   reset flags
    for ( I=0; I<N; I++ ) DR[I][0] = 1.;
    NI=N;
  }
  else {
//   call with negative N - matrix partially inverted
    NI=0;
    for ( I=0; I<N; I++) if ( DR[I][0] != 0.) NI++;
  }

  II=0;
  for ( I=0; I<N; I++) {
    II+=I;
    DR[I][1] = fabsf(V[II]);
  }

//     loop begin, with loop on all remaining rows/cols
  NRANK=N-NI;
  for ( I=1; I<=NI; I++) {
//      search for pivot and test for linearity and zero matrix
    K=JJ=KK=0;
    VKK=0.;
    for ( J=1; J<=N; J++ ) {    
      JJ+=J;
      if ( DR[J-1][0] != 0.) {
        if ( fabsf(V[JJ-1]) > VKK && fabsf(V[JJ-1]) > EPS*DR[J-1][1] ) {
	  VKK=fabsf(V[JJ-1]);
	  K=J;
	  KK=JJ;
	}   
      }
    }

    if ( K==0 || VKK < EPS2 ) {
//         no pivot found, clear of matrix
      IJ=0;
      for ( III=1; III<=N; III++ ) {
        if ( M == 1 && DR[III-1][0] != 0.) B[III-1]=0.;
        for ( J=1; J<=III; J++) {
	  IJ++;
          if ( DR[III-1][0]+DR[J-1][0] != 0.)  V[IJ-1]=0.;
	  V[IJ-1]=-V[IJ-1];
	} 
      }
      return NRANK;
    }
//      preparation for elimination
    NRANK++;
    DR[K-1][0]=0.;
    D=1./V[KK-1];
    if ( fabsf(D) < EPS2 ) D=0.;
    V[KK-1]=-D;
    if (M==1) B[K-1]*=D;
    JK=KK-K;
    JL=0;
//      elimination 
    for ( J=1; J<=N; J++ ) {
      if ( J==K ) {
	JK=KK;
	JL+=J;
      }
      else {
        if ( J<K ) JK++;
        else JK+=J-1;
	E = V[JK-1];
        if ( fabsf(E) < EPS2 ) E=0.;           
	vdbl = D*E;
	if ( fabs(vdbl) < 1.e-31 ) vdbl=0.;
	V[JK-1]=vdbl;
        if (M==1) B[J-1] -= B[K-1]*E;
	LK = KK-K;
        for ( L=1; L<=J; L++ ) {
	  JL++;
          if (L==K) LK=KK;
          else {
            if (L<K) LK++;
            else LK+=L-1;
	    V[JL-1]-=V[LK-1]*E;
	  }
	}
      }
    } 
  } 
//     change sign
  IJ=0;
  for ( I=1; I<=N; I++ ) {
    for ( J=1; J<=I; J++ ) {
      IJ++;
      V[IJ-1]=-V[IJ-1];
    }
  }
  return NRANK;
}

//---------------------------------------------------------------------------
Float_t Vmodf( Float_t x[], Int_t nx )
{
  Float_t vm = 0.;
  Int_t i;
  for (i=0; i<nx; i++) vm += x[i]*x[i];
  if ( vm>0.) vm = sqrtf( vm );
  return vm;
}
Double_t Vmod( Double_t x[], Int_t nx )
{
  Double_t vm = 0.;
  Int_t i;
  for (i=0; i<nx; i++) vm += x[i]*x[i];
  if ( vm>0.) vm = sqrt( vm );
  return vm;
}
void Vsubf( Float_t a[], Float_t b[], Float_t c[], Int_t n )
{
  Int_t i;
  for (i=0; i<n; i++) c[i]= a[i] - b[i];
}
void Vsub( Double_t a[], Double_t b[], Double_t c[], Int_t n )
{
  Int_t i;
  for (i=0; i<n; i++) c[i]= a[i] - b[i];
}
void Ptpxyzf( Float_t vm[3], Float_t ph3[3] )
{
  ph3[0] = vm[0] * sinf(vm[1]) * cosf(vm[2]);
  ph3[1] = vm[0] * sinf(vm[1]) * sinf(vm[2]);
  ph3[2] = vm[0] * cosf(vm[1]);
}
void Ptpxyz( Double_t vm[3], Double_t ph3[3] )
{
  ph3[0] = vm[0] * sin(vm[1]) * cos(vm[2]);
  ph3[1] = vm[0] * sin(vm[1]) * sin(vm[2]);
  ph3[2] = vm[0] * cos(vm[1]);
}
void Pxyztpf( Float_t PXYZ[3], Float_t PTP[3] )
{
  Float_t VPMOD, VPMOD2, COSTHE, COSPHI, PHI;
  const Float_t PI=3.14159265;

  VPMOD  = Vmodf( PXYZ, 3 );
  PTP[0] = VPMOD;
  if ( VPMOD < 0.000001 )  VPMOD = 0.000001;
  COSTHE = PXYZ[2] / VPMOD;
  PTP[1] = acosf(COSTHE);
  VPMOD2 = VPMOD * sqrtf( 1.- COSTHE*COSTHE );
  if ( VPMOD2 < 0.0000001 )  VPMOD2 = 0.0000001; 
  COSPHI = PXYZ[0] / VPMOD2;
  if ( COSPHI > 1.) COSPHI =  1.;          
  if ( COSPHI <-1.) COSPHI = -1.;          
  PHI = acosf( COSPHI );
  if ( PXYZ[1] < 0. ) PHI = PI*2. - PHI; 
  PTP[2] = PHI;
}
void Pxyztp( Double_t PXYZ[3], Double_t PTP[3] )
{
  Double_t VPMOD, VPMOD2, COSTHE, COSPHI, PHI;
  const Double_t PI=3.14159265;

  VPMOD  = Vmod( PXYZ, 3 );
  PTP[0] = VPMOD;
  if ( VPMOD < 0.000001 )  VPMOD = 0.000001;
  COSTHE = PXYZ[2] / VPMOD;
  PTP[1] = acos(COSTHE);
  VPMOD2 = VPMOD * sqrt( 1.- COSTHE*COSTHE );
  if ( VPMOD2 < 0.0000001 )  VPMOD2 = 0.0000001; 
  COSPHI = PXYZ[0] / VPMOD2;
  if ( COSPHI > 1.) COSPHI =  1.;          
  if ( COSPHI <-1.) COSPHI = -1.;          
  PHI = acos( COSPHI );
  if ( PXYZ[1] < 0. ) PHI = PI*2. - PHI; 
  PTP[2] = PHI;
}
/*
void plpxyz( Float_t vm[3], Float_t ph3[3] )
{
  ph3[0] = vm[0] * cosf(vm[1]) * cosf(vm[2]);
  ph3[1] = vm[0] * cosf(vm[1]) * sinf(vm[2]);
  ph3[2] = vm[0] * sinf(vm[1]);
}
void pxyzlp( Float_t PXYZ[3], Float_t PLP[3] )
{
  Float_t VPMOD, VPMOD2, SINLAM, COSPHI, PHI;
  const Float_t PI=3.14159265;

  VPMOD  = Vmod( PXYZ, 3 );
  PLP[0] = VPMOD;
  if ( VPMOD < 0.000001 )  VPMOD = 0.000001;
  SINLAM = PXYZ[2] / VPMOD;
  PLP[1] = asinf(SINLAM);
  VPMOD2 = VPMOD * sqrtf( 1.- SINLAM*SINLAM );
  if ( VPMOD2 < 0.0000001 )  VPMOD2 = 0.0000001; 
  COSPHI = PXYZ[0] / VPMOD2;
  if ( COSPHI > 1.) COSPHI =  1.;          
  if ( COSPHI <-1.) COSPHI = -1.;          
  PHI = acosf( COSPHI );
  IF ( PXYZ[1] < 0. ) PHI = PI*2. - PHI; 
  PLP[2] = PHI;
}
*/
