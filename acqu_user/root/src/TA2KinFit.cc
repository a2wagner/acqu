/**
 * TA2KinFit
 *
 * Wrapper class serving several methods to easily use the KinFitter
 */

#include "TA2KinFit.h"

TA2KinFit::TA2KinFit()
{
}

TA2KinFit::~TA2KinFit()
{
}

int TA2KinFit::fillMatrixDiagonal(TMatrixD* m, Double_t* e, int rows, int cols)
{
	if (m->IsValid()) {
		fprintf(stderr, "Error: Matrix is not valid!\n");
		return 1;
	}
	// set matrix to zero and resize it
	m->Zero();
	m->ResizeTo(rows, cols);
	// get pointer to matrix elements, iterate over them and add array entry if diagonal element
	Double_t *ep = m->GetMatrixArray();
	int idx = 0;
	//memset(ep, 0, m->GetNoElements*sizeof(Double_t));
	memset(ep, 0, rows*cols*sizeof(Double_t));
	for (Int_t i = 0; i < rows; i++)
		for (Int_t j = 0; j < cols; j++)
			*ep++ = (i == j ? e[idx++] : 0.);

	return 0;
}

int TA2KinFit::fillSquareMatrixDiagonal(TMatrixD* m, Double_t* e, int rows)
{
	return fillMatrixDiagonal(m, e, rows, rows);
}

