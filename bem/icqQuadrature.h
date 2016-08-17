/**
 * Quadrature on triangle 
 */

#ifndef ICQ_QUADRATURE
#define ICQ_QUADRATURE

#include <icqFunctor.h>
 
#include <map>
#include <vector>
#include <cmath>

struct icqQuadratureType {
    icqFunctor* func;
    std::map<int, std::vector<std::vector<double> > > gaussPtsAndWeights;
};

extern "C" 
void icqQuadratureInit(icqQuadratureType **self);

extern "C"
void icqQuadratureSetObserver(icqQuadratureType **self, const double* pObs);

extern "C"
void icqQuadratureSetFunctor(icqQuadratureType **self, icqFunctor* func);

extern "C"
int icqQuadratureGetMaxOrder(icqQuadratureType **self);

extern "C"
double icqQuadratureEvaluate(icqQuadratureType **self, int order,
                             const double* pa, const double* pb, const double* pc);

extern "C"
double icqQuadratureEvaluateDouble(icqQuadratureType **self, int order,
                                   const double* pa0, const double* pb0, const double* pc0,
                                   const double* pa1, const double* pb1, const double* pc1);
extern "C"
void icqQuadratureDel(icqQuadratureType **self);

#endif // ICQ_QUADRATURE

