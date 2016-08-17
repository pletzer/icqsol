#include <icqLaplaceFunctor.h>
#include <cmath>

double
icqLaplaceFunctor::operator()(const double* pSrc) {
    double dx = pSrc[0] - this->pObs[0];
    double dy = pSrc[1] - this->pObs[1];
    double dz = pSrc[2] - this->pObs[2];
    return 1.0/sqrt(dx*dx + dy*dy + dz*dz)/(-4. * M_PI);
}