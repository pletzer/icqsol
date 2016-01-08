/**
 * Quadrature on triangle 
 */

#ifndef ICQ_QUADRATURE
#define ICQ_QUADRATURE
 
#include <map>
#include <vector>

 struct icqQuadratureType {
    double (*func)(const double* pos);
    std::map<int, std::vector<std::vector<double> > > gaussPtsAndWeights;
 };

extern "C" 
void icqQuadratureInit(icqQuadratureType **self, double (*f)(const double* pos));

extern "C"
int icqQuadratureGetMaxOrder(icqQuadratureType **self);

extern "C"
double icqQuadratureEvaluate(icqQuadratureType **self, int order,
                            const double* pa, const double*pb, const double* pc);

extern "C"
void icqQuadratureDel(icqQuadratureType **self);

#endif // ICQ_QUADRATURE

