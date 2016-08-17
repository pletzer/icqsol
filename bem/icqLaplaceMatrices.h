#ifndef ICQ_LAPLACE_MATRICES
#define ICQ_LAPLACE_MATRICES

#include <vtkPolyData.h>

extern "C"
void computeOffDiagonalTerms(vtkPolyData* pdata, double* gMat);

#endif // ICQ_LAPLACE_MATRICES