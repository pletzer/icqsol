#ifndef ICQ_FUNCTOR
#define ICQ_FUNCTOR

class icqFunctor {
protected:
    const double* pObs;
public:
    icqFunctor();
    virtual ~icqFunctor();
    void setObserver(const double* pObs);
    virtual double operator()(const double* pSrc) = 0;
};

#endif // ICQ_FUNCTOR