#include <icqFunctor.hpp>

icqFunctor::icqFunctor() {
	this->pObs.resize(3, 0);
}

icqFunctor::~icqFunctor() {}

void
icqFunctor::setObserver(const double* pObs){
	for (size_t i = 0; i < this->pObs.size(); ++i) {
		this->pObs[i] = pObs[i];
	}
}

