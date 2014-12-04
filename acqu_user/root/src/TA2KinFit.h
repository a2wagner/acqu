#ifndef __TA2KinFit_h__
#define __TA2KinFit_h__

#include "KinFitter/TKinFitter.h"
#include "KinFitter/TFitParticlePThetaPhi.h"
#include "KinFitter/TFitConstraintEp.h"
#include "KinFitter/TFitConstraintM.h"

class TA2KinFit : public TKinFitter
{
	private:
	//TODO

	protected:

	public:
	TA2KinFit();
	virtual ~TA2KinFit();
	int fillMatrixDiagonal(TMatrixD* m, Double_t* e, int rows, int cols);
	int fillSquareMatrixDiagonal(TMatrixD* m, Double_t* e, int rows);

	ClassDef(TA2KinFit, 1)
};

#endif  // __TA2KinFit_h__
