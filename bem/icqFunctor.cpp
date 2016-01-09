#include <icqFunctor.hpp>

icqFunctor::icqFunctor() {
	this->pObs = 0;
}

icqFunctor::~icqFunctor() {}

void
icqFunctor::setObserver(const double* pObs){
	this->pObs = pObs;
}

