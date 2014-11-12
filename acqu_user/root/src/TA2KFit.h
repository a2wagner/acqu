//--Author	S Prakhov
//--Description
//                *** Acqu++ <-> Root ***
// Online/Offline Analysis of Sub-Atomic Physics Experimental Data 
//
// TA2KFit
// Kinematical fit routines
//

#ifndef __TA2KFit_h__
#define __TA2KFit_h__

#include "TLorentzVector.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
  
class TA2KFit{
private:
  static const Int_t NXMAX=50, NFMAX=20, NPMAX=20, NDMAX=10; 
  Double_t fChisq;
  Int_t   fNDF;
  Int_t   fTypDat;
  Int_t   fMCsmear;
  TLorentzVector fBeam, fBeamcm;
  TLorentzVector fParticles[NPMAX], fParticlescm[NPMAX];
  Double_t Ecm, Betacm, Velcm[4];  //lab velocities of cm frame
  
  TVector3 fVertex;
  Double_t fPulls[NXMAX];
  Double_t fDeclen;
  
  Double_t fMastg, fMass[NPMAX];
  Float_t fPin[NPMAX][4], fErrSq[NPMAX][4], fCov[10];
  Float_t fSigzvr, fPlim[NXMAX][2], fProut[NPMAX][5];
  Int_t fKind[NPMAX], fCalor[NPMAX], fLcst, fLcstCB, fLcstTAPS;
  Int_t fNpmeas, fNptall;
  void Trevcm();
  TLorentzVector Lorencm( TLorentzVector );

  Int_t NX;      // number of parameters
  Int_t MYF;     // number of transformed parameters/constraint equations
  Int_t NUM;     // flag for numerical differentiation
  Int_t IFLG;    // flag for first case (check derivative)
  Int_t INIT;
  Int_t ND;
  Float_t CHSQ;
  Float_t EPSF;
  Int_t ISTAT;
  Float_t XD[2];     // for current parameter displaced values for differ.
  Float_t XL[NXMAX][2];  // lower and upper values of parameters
  Float_t ST[NXMAX];      // step sizes for numerical differentiation
  Float_t FC[NXMAX];  // central values of parameters 
  Float_t H[NXMAX+NFMAX]; 
  Float_t A[1000];  // derivative matrix a/flags during matrix inversion/pulls
  Float_t DR[100][2];
  Int_t ITER, IDR;
//
//    CHLIM statement function for ...
//    chisquare limit for +k sigma and nd degrees of freedom is
//    approximately
  Float_t CHLIM(Int_t K, Int_t ND) {
    return 0.5*powf((Float_t)K + sqrtf((Float_t)(2*ND+1)),2.);};
  void SIMNIT(Int_t, Int_t, Float_t);
  void SMTOS(Float_t*, Int_t, Float_t*, Int_t, Int_t);
  void SIMLIM(Int_t, Float_t, Float_t);
  void SIMSTP(Int_t, Float_t);
  Int_t APLCON(Float_t*, Float_t*, Float_t*);
  void SIMMAT(Float_t*);
  Int_t SIMDER(Float_t*, Float_t*);
  Int_t SMINVv(Float_t*, Float_t*, Int_t, Int_t);
public:
  TA2KFit();
  virtual ~TA2KFit();
  Int_t Kfilbm(Double_t, Int_t, Double_t*, Double_t*, Double_t*);
  Int_t Kfilcst(Int_t, Int_t, Double_t*);
  Int_t Kfilunm(Int_t, Int_t, Double_t);
  Int_t Kinfit(Int_t, Int_t, Int_t, Int_t*, Double_t*, Int_t*,
	       Int_t*, Int_t Idecay[NDMAX][10], Int_t );
  TLorentzVector NewLVec( TVector3, Double_t, Double_t );
  void Dvpxyz( Int_t, Float_t*, Float_t*, Float_t* );
  void Trevcmf( TLorentzVector, Double_t, Double_t* );
  void Trevrest( TLorentzVector, Double_t* );
  TLorentzVector Loren( Double_t*, TLorentzVector );
  Double_t Zeta();
  void SetResol( Int_t fData, Int_t fMCres ){ fTypDat = fData;
                                              fMCsmear = fMCres; };
  Double_t Chisq(){return fChisq;};
  Int_t NDF(){return fNDF;};
  TLorentzVector Beam(){return fBeam;};
  Double_t  BeamE(){return fBeam.E();};
  TLorentzVector Beamcm(){return fBeamcm;};
  Double_t  BeamcmE(){return fBeamcm.E();};
  TLorentzVector Particle(Int_t ikind, Int_t ip)
{
  TLorentzVector  Empty(0.,0.,0.,0.);
  Int_t np=0, kind;
  for ( Int_t jp=0; jp<fNptall; jp++ ) {
    kind = abs( fKind[jp+1] );
    if ( kind > 1000 ) kind -= 1000;
    if ( kind == ikind ) np++;
    if ( np == ip ) return fParticles[jp]; 
  }
  return Empty;
};
  Double_t  ParticleE(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particle( ikind, ip );
  return Temp.E();
};
  Double_t  ParticleTheta(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particle( ikind, ip );
  return Temp.Theta();
};
  Double_t  ParticlePhi(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particle( ikind, ip );
  return Temp.Phi();
};
  TLorentzVector Particlecm(Int_t ikind, Int_t ip)
{
  TLorentzVector  Empty(0.,0.,0.,0.);
  Int_t np=0, kind;
  for ( Int_t jp=0; jp<fNptall; jp++ ) {
    kind = abs( fKind[jp+1] );
    if ( kind > 1000 ) kind -= 1000;
    if ( kind == ikind ) np++;
    if ( np == ip ) return fParticlescm[jp]; 
  }
  return Empty;
};
  Double_t  ParticlecmE(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particlecm( ikind, ip );
  return Temp.E();
};
  Double_t  ParticlecmTheta(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particlecm( ikind, ip );
  return Temp.Theta();
};
  Double_t  ParticlecmCosTheta(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particlecm( ikind, ip );
  return Temp.CosTheta();
};
  Double_t  ParticlecmPhi(Int_t ikind, Int_t ip)
{
  TLorentzVector Temp = Particlecm( ikind, ip );
  return Temp.Phi();
};
  Double_t* Pulls(){return fPulls;};
  Double_t PullEbeam(){return fPulls[0];};
  Double_t PullXvrtx(){return fPulls[1];};
  Double_t PullYvrtx(){return fPulls[2];};
  Double_t PullZvrtx(){return fPulls[3];};
  Double_t* PullsParticle(Int_t ikind, Int_t ip)
{
  static Double_t ZerPl[4] = { 4*0.};
  Int_t np=0, kind;
  for ( Int_t jp=1; jp<=fNpmeas; jp++ ) {
    kind = abs( fKind[jp] );
    if ( kind == ikind ) np++;
    if ( np == ip ) return fPulls+jp*4; 
  }
  return ZerPl;
};
  Double_t PullParticleE(Int_t ikind, Int_t ip)
{
  Double_t* Pl = PullsParticle( ikind, ip );
  return *Pl;
};
  Double_t PullParticleTheta(Int_t ikind, Int_t ip)
{
  Double_t* Pl = PullsParticle( ikind, ip );
  return *(Pl+1);
};
  Double_t PullParticlePhi(Int_t ikind, Int_t ip)
{
  Double_t* Pl = PullsParticle( ikind, ip );
  return *(Pl+2);
};
  Double_t PullParticleDepth(Int_t ikind, Int_t ip)
{
  Double_t* Pl = PullsParticle( ikind, ip );
  return *(Pl+3);
};
  TVector3 Vertex(){return fVertex;};
  Double_t VertexX(){return fVertex.X();};
  Double_t VertexY(){return fVertex.Y();};
  Double_t VertexZ(){return fVertex.Z();};
  Double_t Declen(){return fDeclen;};

//Double_t DepthGam(Double_t Eg)           // depth of e/m shower
//  {return fRped0+fRadlnai*fVhalfe*(log(Eg/fEcrit)+0.5);};

  Double_t DepthShowCB(Double_t Ecl, Int_t IPart) {  // a shower depth in CB
    Double_t dep;
    Double_t p[4] = {-3.36631,9.40334e-02,5.35372e-01,4.36397e+01}; // photon
    //Double_t p[4] = {-3.97217,7.63539e-02,4.54722e-01,4.42913e+01}; // photon
   
    if ( IPart == 1 ) {  // photon
      dep = p[0]/pow(Ecl+p[1],p[2])+p[3]; 
    }
    else if ( IPart == 14 ) {   // proton
      p[0] = 2.52512e+01; p[1] = 6.44248; p[2] = 1.96292e+02; p[3] = -1.61958e+02;
      dep = p[0]+Ecl*p[1]+Ecl*Ecl*p[2]+Ecl*Ecl*Ecl*p[3];
    }
    else dep = 45.;

    return dep;
};
  Double_t dDepthShowCB(Double_t Ecl, Int_t IPart) {  // uncertaity of shower depth in CB
    Double_t sdep;
    Double_t p[4] = {1.76634e-01,0.,6.26983e-01,2.48218}; // photon
    //Double_t p[4] = {1.29312e-01,0.,7.24244e-01,2.54382}; // photon
    if ( IPart == 1 ) {
      sdep = p[0]/pow(Ecl+p[1],p[2])+p[3];
      sdep *= 1.05; 
    }
    else if ( IPart == 14 ) {  // proton
      p[0] = 3.5783e-02; p[1] = 3.47172e-01; p[2] = 1.50307; p[3] = -4.88434e-01;
      sdep = p[0]+Ecl*p[1]+Ecl*Ecl*p[2]+Ecl*Ecl*Ecl*p[3];
      sdep *= 1.05; 
    }
    else sdep = 4.;

    return sdep;
};
  Double_t DepthShowTAPS(Double_t Ecl, Int_t IPart) {  // a shower depth in TAPS
    Double_t p[4] = {-2.99791e+01,1.75852e-03,4.99643e-02,4.14362e+01}; // photon gauss fit
    //Double_t p[4] = {-4.53249e+01,0.,3.35651e-02,5.68095e+01}; // photon gauss fit
    if ( IPart == 1 ) { // photon
      //p[0] *= 0.9425; // gp
      p[0] *= 0.978; // 0.98 // pi0
      return p[0]/pow(Ecl+p[1],p[2])+p[3];
    } 
    else if ( IPart == 14 ) {   // proton
      p[0] = -1.73216e-02; p[1] = 3.83753; p[2] = 1.54891e+02; p[3] =-1.328e+02;
      Double_t dprot = p[0]+Ecl*p[1]+Ecl*Ecl*p[2]+Ecl*Ecl*Ecl*p[3];
      //return dprot *1.025;  //gp
      return dprot *1.05;  //pi0
    }
    else return 12.5;
};

  Double_t dDepthShowTAPS(Double_t Ecl, Int_t IPart) {  // uncertaity of shower depth in TAPS
    Double_t p[4] = {2.83139,0.,1.02537e-01,-7.53507e-01}; // photon gauss sigmas fit
    //Double_t p[4] = {2.09183,0.,1.2527e-01,9.83775e-05}; // photon gauss sigmas fit
    if ( IPart == 1 ) return p[0]/pow(Ecl+p[1],p[2])+p[3]; 
    else if ( IPart == 14 ) {   // proton
      p[0] = 8.43187e-03; p[1] = 3.63264e-01; p[2] = 7.17476e-01; p[3] = 7.33715;
      return p[0]+Ecl*p[1]+Ecl*Ecl*p[2]+Ecl*Ecl*Ecl*p[3];
    }
    else return 4.;
};

  Double_t EclCorCB(Double_t Ecl, Int_t IPart, Int_t IRaw) {            // CB energy correction
    // correction for 12 MeV cluster threshold in the CB.
    //Double_t p[5] = {1.76000e-02,0.,3.76874e-01,2.43402e-03,1.10303e-02};// photon unsmeared
              //p[0] *= 0.93;   // g
              //p[3] *= 1.75;
    //p[0] *= 1.09;   //pi0
    //p[3] *= 0.7; 
    //Double_t p[5] = {1.34997e-01,0.,1.25340e-01,-1.20832e-01,1.72907e-02};// photon smeared s1
    //Double_t p[5] = {7.07721e-04,7.0e-01,1.30677e+01,4.12382e-02,-1.22663e-02};// photon smeared s2
    //Double_t p[5] = {1.14683e-04,2.48336e-01,4.657,4.01373e-02,-1.08977e-02};// photon smeared s2
    //Double_t p[5] = {4.32764e-03,3.04581e-02,9.20136e-01,2.60488e-02,1.34082e-03};// photon smeared s3
    //Double_t p[5] = {1.19028e-02,8.16827e-03,5.14487e-01,1.50047e-02,5.18801e-03};// photon smeared tasm2
    Double_t p[5] = {1.52950e-02,5.92991e-03,4.57857e-01,8.98180e-03,7.75002e-03};// photon smeared smcal11
    //p[0] *= 0.88;   //pi0
    //p[4] = 0.0086; //?
    //p[4] = 0.009; //res2
    //if ( IRaw != 0 ) {
    //p[0] = 3.97514e-05; p[1] = 0.; p[2] = 1.66502; p[3] = 5.36403e-02; p[4] = -2.62643e-02;
    //}
    if ( IPart == 14 ) {   // proton
      p[4]=0.;
      if ( IRaw !=0 ) {
	if ( Ecl<0.08 ) {p[0]=0.43; p[1]=-0.01; p[2]=0.4; p[3]=-1.18;} 
	else {p[0]=0.16; p[1] = -0.07; p[2] = 0.15; p[3] = -0.25; }
      }
      else {
	if ( Ecl<0.08 ) {p[0]=0.4; p[1]=-0.01; p[2]=0.38; p[3]=-0.9;} 
	else {p[0]=0.18; p[1] = -0.06; p[2] = 0.2; p[3] = -0.21; }
      }
      if ( Ecl<0.011 ) Ecl=0.011; 
    }
    Double_t Fcor = p[0]/pow(Ecl+p[1],p[2])+p[3]+p[4]*Ecl;
    //Fcor = 0.053;
    //Fcor *= 0.95;
    return Fcor;
};

  Double_t dEovEclCB(Double_t Ecl, Int_t IPart) {                        // CB dE/E
    Double_t Er0 = dEovEclCBInit( Ecl, IPart);
    Double_t Era = dEovEclCBAdd( Ecl, IPart);
    //Era = 0.;
    //printf("CB %lf %lf %lf %lf\n",Ecl,Er0,Era,sqrt( pow(Er0,2) + pow(Era,2) ) );
    return sqrt( pow(Er0,2) + pow(Era,2) );
};

  Double_t dEovEclCBInit(Double_t Ecl, Int_t IPart) {                        // CB dE/E
    Double_t p[5] = {5.69464e-05,1.48943e-01,3.41725,1.11244e-02,-1.77329e-03};
    /* 
    p[0] *= 1.19; //old
    p[3] *= 1.6;
    p[2] *= 1.01;
    */
    /*
    p[0] *= 0.2; // res0
    p[1] = 0.148;
    p[2] *= 1.3;
    p[3] *= 1.93;
    p[4] = 0.002;
    */
    
    p[0] = 0.014; // res1
    p[1] = 0.0025;
    p[2] = 0.35;
    p[3] = 0.;
    p[4] = 0.0032;
        
    if ( IPart == 14 ) {   // proton
      p[0] = 0.043; p[1] = 0.; p[2] = 0.43; p[3] = 0.; p[4] = 0.;
    }
    Double_t Er0 = p[0]/pow(Ecl+p[1],p[2])+p[3]+p[4]*Ecl;
    //Double_t Er0 = p[0]/pow(Ecl+p[1],p[2])+p[3]; //?
    return Er0;
};

  Double_t dEovEclCBAdd(Double_t Ecl, Int_t IPart) {                        // CB dE/E
    //Double_t Era = 0.0145/pow(Ecl,0.34);  // DMM runs 3811-3859
    //Double_t Era = 0.0145/pow(Ecl,0.28);  // res3
    //Double_t Era = 0.0145/pow(Ecl,0.1);  // res4
    Double_t Era = 0.0145/pow(Ecl-0.01,0.1);  // res1
    //Era += 0.004*Ecl; // res4
    Era += 0.0035*Ecl; // res4
    //Era += 0.0075*Ecl; // res5
    //Era += 0.018*Ecl; // res3
    if( fTypDat==2 ) Era += 0.007;  // Eta runs
    if( fTypDat==4 ) Era += 0.003;  // Eta runs 2007, pi0 2008
    //if( fTypDat==4 ) Era += 0.0065;  // Eta runs 2007, pi0 2008
    if( fTypDat==5 ) Era += 0.007;  // Eta runs April 09
    if( fTypDat==6 ) Era += 0.007;  // Butanol runs May-June 10
    //if( fTypDat==7 ) Era *= 2.0;  // EPT runs May-August 12
    //if( fTypDat==7 ) Era *= 1.8;  // EPT runs May-August 12
    if( fTypDat==7 ) Era += 0.022;  // EPT runs May-August 12
    if( fTypDat==8 ) Era += 0.009;  // December 12, March 13
    //if( fTypDat==9 ) Era += 0.0085;  // April 13 res3
    //if( fTypDat==9 ) Era += 0.02;  // April 13 res0
    //if( fTypDat==9 ) Era = 0.033+0.0045*Ecl;  // April 13 res1
    if( fTypDat==9 ) Era = 0.032+0.009*Ecl;  // April 13 res1
    //if( fTypDat==6 ) Era = 0.032+0.009*Ecl;  // Nov 13 butanol
    if( fTypDat==4 ) Era = 0.0255;  // Eta runs 2007 , pi0 2008
    if( fTypDat==6 ) Era = 0.028;  // runs 21860-21895 April 09
    //if( fTypDat==5 ) Era = 0.0305;  // Eta runs April 09
    if( fTypDat==5 ) Era = 0.0295+0.0035*Ecl;  // Eta runs April 09
    //if( fTypDat==10 ) Era = 0.032+0.009*Ecl;  // EPT 2014 temp
    //if( fTypDat==10 ) Era = 0.036+0.0035*Ecl;  // EPT 2014 res1
    //if( fTypDat==10 ) Era = 0.040+0.0025*Ecl;  // EPT 2014 res2
    //if( fTypDat==10 ) Era = 0.045+0.0015*Ecl;  // EPT 2014 res3
    if( fTypDat==10 ) Era = 0.052;  // EPT 2014 res4
    //if( fTypDat==7 ) Era = 0.040+0.0025*Ecl;  // EPT 2012 res1
    if( fTypDat==7 ) Era = 0.045+0.0015*Ecl;  // EPT 2012 res2
     if ( IPart == 14 ) {   // proton
      Era = 0.;
    }
    if ( fMCsmear==0 || fMCsmear==2 )  Era = 0.;
    return Era;
};
  Double_t ELossPID(Double_t Ecl, Double_t Thetacl, Int_t IPart) {
    Double_t p[4] = {7.71580e-02,1.02080e-01,1.96421,9.59329e-01};
    if ( IPart == 14 ) {   // proton
      p[0]=7.71580e-02; p[1]=1.02080e-01; p[2]=1.96421; p[3]=9.59329e-01;
    }
    Double_t Els = p[0]/pow(Ecl+p[1],p[2])+p[3];
    return Els/1000./sin(Thetacl);
};
  Double_t EclCorTAPS(Double_t Ecl, Int_t IPart, Int_t IRaw) {            // TAPS cluster energy correction
    // correction for 12 MeV cluster threshold.
    //Double_t p[4] = {0.090434,-0.308149,2.040689,-3.104840};  // photon smeared 0.05
    //if ( Ecl>0.32 ) { p[0] = 0.065230; p[1] = 0.124692; p[2] = -0.068446; p[3] = 0.; }
    //Double_t Fcor = p[0];
    //for ( Int_t i=1; i<4; i++ ) Fcor += pow(Ecl,i)*p[i];
    // Double_t p[4] = {2.67035e-04,1.19147,9.26359e-02,5.40767e-02}; //for err=0.05
    //Double_t p[4] = {4.16557e-02,1.72805e-01,3.27917e-02,6.19942e-02}; // for err=v2
    Double_t p[4] = {1.07370e-02,4.74417e-01,6.38448e-02,5.73663e-02}; // for err=v3
    //Double_t Fcor = p[0]/pow(Ecl,p[1])+p[2]+p[3]*Ecl;
    //Double_t Fcor = p[0]/pow(Ecl-0.009,p[1])+p[2]+p[3]*Ecl;
    Double_t Fcor = p[0]/pow(Ecl-0.005,p[1])+p[2]+p[3]*Ecl;
    if ( IPart == 14 ) {   // proton
      if ( IRaw !=0 ) {
	if ( Ecl<0.13 ) {p[0] = 0.38; p[1] = -0.012; p[2] = 0.4; p[3] =-0.78;}
	else { p[0] = 0.23; p[1] = 0.; p[2] = 0.25; p[3] = -0.275; }
      }
      else {
	if ( Ecl<0.1 ) {p[0] = 0.33; p[1] = -0.012; p[2] = 0.37; p[3] =-0.7;}
	else { p[0] = 0.23; p[1] = 0.; p[2] = 0.21; p[3] = -0.275; }
      }
      if ( Ecl<0.0125 ) Ecl=0.0125; 
      Fcor = p[0]/pow(Ecl+p[1],p[2])+p[3];
      if ( Ecl>0.4 ) Fcor = 0.;
    }
    return Fcor;
};

  Double_t dEovEclTAPS(Double_t Ecl, Int_t IPart) {                        // TAPS dE/E
    Double_t Er0 = dEovEclTAPSInit( Ecl, IPart );
    Double_t Era = dEovEclTAPSAdd( Ecl, IPart );
    //Era = 0.;
    //printf("taps %lf %lf %lf\n",Ecl,Er0,Era );
    return sqrt( pow(Er0,2) + pow(Era,2) );
};

  Double_t dEovEclTAPSInit(Double_t Ecl, Int_t IPart) {                        // TAPS dE/E
    //Double_t p[7] = {0.056157,-0.296396,2.975889,-12.064184,23.848432,-22.856702,8.504981};  // photon
    //Double_t Er0 = p[0]+Ecl*p[1]+0.005;
    //for ( Int_t i=2; i<7; i++ ) Er0 += pow(Ecl,i)*p[i];
    Double_t p[4] = {1.88319e-04,1.42657,3.96356e-02,1.52351e-02};
    //Double_t Er0 = p[0]/pow(Ecl,p[1])+p[2]+p[3]*Ecl;
    p[3] *= 1.8;
    Double_t Er0 = p[0]/pow(Ecl-0.002,p[1])+p[2]+p[3]*Ecl;
    //Double_t Er0 = p[0]/pow(Ecl-0.0025,p[1])+p[2]+p[3]*Ecl;
    if ( IPart == 14 ) {   // proton
      p[0] = 0.045; p[1] = 0.; p[2] = 0.45; p[3] = 0.;
      Er0 = p[0]/pow(Ecl+p[1],p[2])+p[3];
    }
    return Er0;
};

  Double_t dEovEclTAPSAdd(Double_t Ecl, Int_t IPart) {                        // TAPS dE/E
    Double_t Er0 = dEovEclTAPSInit( Ecl, IPart );
    //Double_t Era = 0.05;  // Eta' runs
    //Double_t Era = 0.015+0.012/pow(Ecl,0.4);;  // Eta' runs v1
    //Double_t Era = 0.015/pow(Ecl,0.485);  // Eta' runs v2
    //Double_t Era = 0.015/pow(Ecl-0.006,0.5);;  // res2
    //if( fTypDat>2 ) Era = Er0*(1.+Ecl*0.3);  // Eta' runs
    //Era += 0.015*Ecl; // res5
    //Double_t Era = 0.015/pow(Ecl+0.014,0.495)+0.012*Ecl;  // pi0 Apr'13 runs v1
    //Double_t Era = 0.011/pow(Ecl+0.014,0.495)+0.012*Ecl+0.01;  // pi0 Apr'13 runs v2
    //Double_t Era = 0.008/pow(Ecl+0.025,0.495)+0.008*Ecl+0.018;  // pi0 Apr'13 runs v3
    //Double_t Era = 0.035+0.02*Ecl;  // pi0 Apr'13    
    Double_t Era = 0.031+0.04*Ecl;  // pi0 Apr'13    
    if ( IPart == 14 ) {   // proton
      Era = 0.;
    }
    if ( fMCsmear==0 || fMCsmear==1 )  Era = 0.;
    return Era;
};

  Double_t dThetaCB(Double_t Ecl, Int_t IPart) {  // dTheta in CB
    Double_t p[4] = {7.69518e-03,4.86197e-01,1.79483,1.57948e-02}; // photon
    //Double_t p[4] = {1.70586,3.27800e-02,3.91050e-03,-1.68662}; // photon
    if ( IPart == 14 ) {   // proton
      //p[0] = 8.77387e-04; p[1] = 7.65552e-01; p[2] = 1.05023e+01; p[3] = 3.71969e-02; //angle cut
      p[0] = 1.38476e-04; p[1] = 5.30098e-01; p[2] = 7.61558; p[3] = 3.75841e-02; //all angles
      p[3] += 0.004;
    }
    return p[0]/pow(Ecl+p[1],p[2])+p[3];
};

  Double_t dTanThTAPS(Double_t Ecl, Int_t IPart) {  // error in tan(Theta) in TAPS
    Double_t p[5] = {3.28138e+02,0.,7.29002e-04,-3.27381e+02,0.}; // photon
    Double_t dtan = p[0]/pow(Ecl+p[1],p[2])+p[3];
    dtan *= 0.85;
    if ( IPart == 14 ) {   // proton
      //p[0] = 3.915273; p[1] = -19.482241; p[2] = 112.849986, p[3] = -365.957825; p[4] = 424.107806;
      //dtan = p[0]+Ecl*p[1]+Ecl*Ecl*p[2]+Ecl*Ecl*Ecl*p[3]+Ecl*Ecl*Ecl*Ecl*p[4];
      p[0] = 3.27709e+02; p[1] = 4.99670e-02; p[2] = 5.55520e-03; p[3] = -3.27819e+02;
      dtan = p[0]/pow(Ecl+p[1],p[2])+p[3];
      //dtan += 1.5;
    }
    return dtan;
};

   Double_t EKinProt(Double_t Ecl) {                  // proton E correction
          return Ecl+0.018/pow(Ecl+0.012,0.36)-0.0315;};
  //          return Ecl+0.0165/powf(Ecl+0.012,0.36)-0.0281;};
  Double_t dThetaProt(Double_t Ecl, Double_t Tcl) {        // proton dTheta
    Double_t Tclg = Tcl*TMath::RadToDeg();
    if ( Tclg<27. ) return 0.639357/pow(Ecl+1.6797,7.78758)+0.0387287;
    else return 0.00017517/pow(Ecl+0.0464427,1.58624)+0.0357436;};
  Double_t dEovEKinProt(Double_t Ecl) {                   // proton dE/Ekin
          return (0.011472+0.051288*Ecl)/EKinProt(Ecl);};
  Double_t dEProt(Double_t Ecl) {                   // proton dE 
    return (0.006/pow(Ecl-0.00733,0.588)-0.002)*Ecl;};
  //return (0.00466/pow(Ecl-0.00733,0.588)+0.00351)*Ecl;};
};

#endif
