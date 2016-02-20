#ifndef ICQ_FUNCTOR
#define ICQ_FUNCTOR

#include <vector>

class icqFunctor {
protected:
    std::vector<double> pObs;
public:
    icqFunctor();
    virtual ~icqFunctor();
    void setObserver(const double* pObs);
    virtual double operator()(const double* pSrc) = 0;
};

#endif // ICQ_FUNCTOR