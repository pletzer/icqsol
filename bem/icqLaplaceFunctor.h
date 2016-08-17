#ifndef ICQ_LAPLACE_FUNCTOR
#define ICQ_LAPLACE_FUNCTOR

#include <icqFunctor.h>

class icqLaplaceFunctor : public icqFunctor {
public:
    double operator()(const double* pSrc);
};

#endif // ICQ_LAPLACE_FUNCTOR
