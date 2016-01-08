/**
 * Quadrature on triangle 
 */

#include <cmath>
#include <icqQuadrature.hpp>

extern "C" 
void icqQuadratureInit(icqQuadratureType **self, double (*f)(const double* pos)) {

    // Allocate 
    *self = new icqQuadratureType();

    (*self)->func = f;

    // Fill in the gauss points and weights
    int order;

    {
        order = 1;
        std::vector<double> xsi(1), eta(1), wgh(1);
        xsi[0] = 0.33333333333333; eta[0] = 0.33333333333333; wgh[0] = 1.00000000000000;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 2;
        std::vector<double> xsi(3), eta(3), wgh(3);
        xsi[0] = 0.16666666666667; eta[0] = 0.16666666666667; wgh[0] = 0.33333333333333;
        xsi[1] = 0.16666666666667; eta[1] = 0.66666666666667; wgh[1] = 0.33333333333333;
        xsi[2] = 0.66666666666667; eta[2] = 0.16666666666667; wgh[2] = 0.33333333333333;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

}

extern "C"
int icqQuadratureGetMaxOrder(icqQuadratureType **self) {
    return (int) (*self)->gaussPtsAndWeights.size();
}

extern "C"
double icqQuadratureEvaluate(icqQuadratureType **self, int order,
                             const double* pa, const double*pb, const double* pc) {

    double res = 0;
    const double pb2[] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
    const double pc2[] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};

    const double p12 = pb2[1]*pc2[2];
    const double p20 = pb2[2]*pc2[0];
    const double p01 = pb2[0]*pc2[1];

    const double p21 = pb2[2]*pc2[1];
    const double p02 = pb2[0]*pc2[2];
    const double p10 = pb2[1]*pc2[0];

    const double area = sqrt((p12-p21)*(p12-p21) + 
                             (p20-p02)*(p20-p02) + 
                             (p01-p10)*(p01-p10));

    if (area == 0) return res;

    double p[] = {0., 0., 0.};

    std::map<int, std::vector< std::vector<double> > >::const_iterator 
        it = (*self)->gaussPtsAndWeights.find(order);
    if (it != (*self)->gaussPtsAndWeights.end()) {
        const std::vector<double>& xsi = it->second[0];
        const std::vector<double>& eta = it->second[1];
        const std::vector<double>& wgh = it->second[2];
        for (size_t i = 0; i < xsi.size(); ++i) {
            for (size_t j = 0; j < 3; ++j)
                p[j] = pa[j] + xsi[i]*pb2[j] + eta[i]*pc2[j];
            res += (*self)->func(p) * wgh[i];
        }
    }

    return 0.5 * area * res;
}

extern "C"
void icqQuadratureDel(icqQuadratureType **self) {
    
    delete (*self);
}

