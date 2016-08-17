/**
 * Quadrature on triangle 
 */

#include <icqQuadrature.h>
#include <icqLaplaceFunctor.h>
#include <iostream>

extern "C" 
void icqQuadratureInit(icqQuadratureType **self) {

    // Allocate 
    *self = new icqQuadratureType();

    (*self)->func = new icqLaplaceFunctor();

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

    {
        order = 3;
        std::vector<double> xsi(4), eta(4), wgh(4);
        xsi[0] = 0.33333333333333; eta[0] = 0.33333333333333; wgh[0] = -0.56250000000000;
        xsi[1] = 0.20000000000000; eta[1] = 0.20000000000000; wgh[1] = 0.52083333333333;
        xsi[2] = 0.20000000000000; eta[2] = 0.60000000000000; wgh[2] = 0.52083333333333;
        xsi[3] = 0.60000000000000; eta[3] = 0.20000000000000; wgh[3] = 0.52083333333333;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 4;
        std::vector<double> xsi(6), eta(6), wgh(6);
        xsi[0] = 0.44594849091597; eta[0] = 0.44594849091597; wgh[0] = 0.22338158967801;
        xsi[1] = 0.44594849091597; eta[1] = 0.10810301816807; wgh[1] = 0.22338158967801;
        xsi[2] = 0.10810301816807; eta[2] = 0.44594849091597; wgh[2] = 0.22338158967801;
        xsi[3] = 0.09157621350977; eta[3] = 0.09157621350977; wgh[3] = 0.10995174365532;
        xsi[4] = 0.09157621350977; eta[4] = 0.81684757298046; wgh[4] = 0.10995174365532;
        xsi[5] = 0.81684757298046; eta[5] = 0.09157621350977; wgh[5] = 0.10995174365532;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 5;
        std::vector<double> xsi(7), eta(7), wgh(7);
        xsi[0] = 0.33333333333333; eta[0] = 0.33333333333333; wgh[0] = 0.22500000000000;
        xsi[1] = 0.47014206410511; eta[1] = 0.47014206410511; wgh[1] = 0.13239415278851;
        xsi[2] = 0.47014206410511; eta[2] = 0.05971587178977; wgh[2] = 0.13239415278851;
        xsi[3] = 0.05971587178977; eta[3] = 0.47014206410511; wgh[3] = 0.13239415278851;
        xsi[4] = 0.10128650732346; eta[4] = 0.10128650732346; wgh[4] = 0.12593918054483;
        xsi[5] = 0.10128650732346; eta[5] = 0.79742698535309; wgh[5] = 0.12593918054483;
        xsi[6] = 0.79742698535309; eta[6] = 0.10128650732346; wgh[6] = 0.12593918054483;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 6;
        std::vector<double> xsi(12), eta(12), wgh(12);
        xsi[0] = 0.24928674517091; eta[0] = 0.24928674517091; wgh[0] = 0.11678627572638;
        xsi[1] = 0.24928674517091; eta[1] = 0.50142650965818; wgh[1] = 0.11678627572638;
        xsi[2] = 0.50142650965818; eta[2] = 0.24928674517091; wgh[2] = 0.11678627572638;
        xsi[3] = 0.06308901449150; eta[3] = 0.06308901449150; wgh[3] = 0.05084490637021;
        xsi[4] = 0.06308901449150; eta[4] = 0.87382197101700; wgh[4] = 0.05084490637021;
        xsi[5] = 0.87382197101700; eta[5] = 0.06308901449150; wgh[5] = 0.05084490637021;
        xsi[6] = 0.31035245103378; eta[6] = 0.63650249912140; wgh[6] = 0.08285107561837;
        xsi[7] = 0.63650249912140; eta[7] = 0.05314504984482; wgh[7] = 0.08285107561837;
        xsi[8] = 0.05314504984482; eta[8] = 0.31035245103378; wgh[8] = 0.08285107561837;
        xsi[9] = 0.63650249912140; eta[9] = 0.31035245103378; wgh[9] = 0.08285107561837;
        xsi[10] = 0.31035245103378; eta[10] = 0.05314504984482; wgh[10] = 0.08285107561837;
        xsi[11] = 0.05314504984482; eta[11] = 0.63650249912140; wgh[11] = 0.08285107561837;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 7;
        std::vector<double> xsi(13), eta(13), wgh(13);
        xsi[0] = 0.33333333333333; eta[0] = 0.33333333333333; wgh[0] = -0.14957004446768;
        xsi[1] = 0.26034596607904; eta[1] = 0.26034596607904; wgh[1] = 0.17561525743321;
        xsi[2] = 0.26034596607904; eta[2] = 0.47930806784192; wgh[2] = 0.17561525743321;
        xsi[3] = 0.47930806784192; eta[3] = 0.26034596607904; wgh[3] = 0.17561525743321;
        xsi[4] = 0.06513010290222; eta[4] = 0.06513010290222; wgh[4] = 0.05334723560884;
        xsi[5] = 0.06513010290222; eta[5] = 0.86973979419557; wgh[5] = 0.05334723560884;
        xsi[6] = 0.86973979419557; eta[6] = 0.06513010290222; wgh[6] = 0.05334723560884;
        xsi[7] = 0.31286549600487; eta[7] = 0.63844418856981; wgh[7] = 0.07711376089026;
        xsi[8] = 0.63844418856981; eta[8] = 0.04869031542532; wgh[8] = 0.07711376089026;
        xsi[9] = 0.04869031542532; eta[9] = 0.31286549600487; wgh[9] = 0.07711376089026;
        xsi[10] = 0.63844418856981; eta[10] = 0.31286549600487; wgh[10] = 0.07711376089026;
        xsi[11] = 0.31286549600487; eta[11] = 0.04869031542532; wgh[11] = 0.07711376089026;
        xsi[12] = 0.04869031542532; eta[12] = 0.63844418856981; wgh[12] = 0.07711376089026;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

    {
        order = 8;
        std::vector<double> xsi(16), eta(16), wgh(16);
        xsi[0] = 0.33333333333333; eta[0] = 0.33333333333333; wgh[0] = 0.14431560767779;
        xsi[1] = 0.45929258829272; eta[1] = 0.45929258829272; wgh[1] = 0.09509163426728;
        xsi[2] = 0.45929258829272; eta[2] = 0.08141482341455; wgh[2] = 0.09509163426728;
        xsi[3] = 0.08141482341455; eta[3] = 0.45929258829272; wgh[3] = 0.09509163426728;
        xsi[4] = 0.17056930775176; eta[4] = 0.17056930775176; wgh[4] = 0.10321737053472;
        xsi[5] = 0.17056930775176; eta[5] = 0.65886138449648; wgh[5] = 0.10321737053472;
        xsi[6] = 0.65886138449648; eta[6] = 0.17056930775176; wgh[6] = 0.10321737053472;
        xsi[7] = 0.05054722831703; eta[7] = 0.05054722831703; wgh[7] = 0.03245849762320;
        xsi[8] = 0.05054722831703; eta[8] = 0.89890554336594; wgh[8] = 0.03245849762320;
        xsi[9] = 0.89890554336594; eta[9] = 0.05054722831703; wgh[9] = 0.03245849762320;
        xsi[10] = 0.26311282963464; eta[10] = 0.72849239295540; wgh[10] = 0.02723031417443;
        xsi[11] = 0.72849239295540; eta[11] = 0.00839477740996; wgh[11] = 0.02723031417443;
        xsi[12] = 0.00839477740996; eta[12] = 0.26311282963464; wgh[12] = 0.02723031417443;
        xsi[13] = 0.72849239295540; eta[13] = 0.26311282963464; wgh[13] = 0.02723031417443;
        xsi[14] = 0.26311282963464; eta[14] = 0.00839477740996; wgh[14] = 0.02723031417443;
        xsi[15] = 0.00839477740996; eta[15] = 0.72849239295540; wgh[15] = 0.02723031417443;
        std::vector<std::vector<double> > gw(3); 
        gw[0] = xsi; gw[1] = eta; gw[2] = wgh;
        (*self)->gaussPtsAndWeights.insert(std::pair<int, std::vector< std::vector<double> > >(order, gw));
    }

}

extern "C"
void icqQuadratureSetObserver(icqQuadratureType **self, const double* pObs) {
    (*self)->func->setObserver(pObs);
}

extern "C"
void icqQuadratureSetFunctor(icqQuadratureType **self, icqFunctor* func) {
    (*self)->func = func;
}

extern "C"
int icqQuadratureGetMaxOrder(icqQuadratureType **self) {
    return (int) (*self)->gaussPtsAndWeights.size();
}

double icqQuadratureArea(const double* db, const double* dc) {
    const double p12 = db[1]*dc[2];
    const double p20 = db[2]*dc[0];
    const double p01 = db[0]*dc[1];
    const double p21 = db[2]*dc[1];
    const double p02 = db[0]*dc[2];
    const double p10 = db[1]*dc[0];
    const double area = sqrt((p12-p21)*(p12-p21) +
                             (p20-p02)*(p20-p02) +
                             (p01-p10)*(p01-p10));
    return area;
}

extern "C"
double icqQuadratureEvaluate(icqQuadratureType **self, int order,
                             const double* pa, const double*pb, const double* pc) {

    double res = 0;
    const double pb2[] = {pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
    const double pc2[] = {pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
    double area = icqQuadratureArea(pb2, pc2);

    if (area == 0) return res;

    double p[] = {0., 0., 0.};

    std::map<int, std::vector< std::vector<double> > >::const_iterator 
        it = (*self)->gaussPtsAndWeights.find(order);

    if (it != (*self)->gaussPtsAndWeights.end()) {

        const std::vector<double>& xsi = it->second[0];
        const std::vector<double>& eta = it->second[1];
        const std::vector<double>& wgh = it->second[2];
        for (size_t i = 0; i < xsi.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                p[j] = pa[j] + xsi[i]*pb2[j] + eta[i]*pc2[j];
            }
            res += (*self)->func->operator()(p) * wgh[i];
        }
    }

    return 0.5 * area * res;
}

extern "C"
double icqQuadratureEvaluateDouble(icqQuadratureType **self, int order,
                             const double* pa0, const double* pb0, const double* pc0,
                             const double* pa1, const double* pb1, const double* pc1) {
    double res = 0;
    const double db0[] = {pb0[0] - pa0[0], pb0[1] - pa0[1], pb0[2] - pa0[2]};
    const double dc0[] = {pc0[0] - pa0[0], pc0[1] - pa0[1], pc0[2] - pa0[2]};

    const double db1[] = {pb1[0] - pa1[0], pb1[1] - pa1[1], pb1[2] - pa1[2]};
    const double dc1[] = {pc1[0] - pa1[0], pc1[1] - pa1[1], pc1[2] - pa1[2]};

    //const double area0 = icqQuadratureArea(db0, dc0);
    const double area1 = icqQuadratureArea(db1, dc1);

    double p0[] = {0., 0., 0.};
    double p1[] = {0., 0., 0.};

    std::map<int, std::vector< std::vector<double> > >::const_iterator
        it = (*self)->gaussPtsAndWeights.find(order);
    if (it != (*self)->gaussPtsAndWeights.end()) {
        const std::vector<double>& xsi = it->second[0];
        const std::vector<double>& eta = it->second[1];
        const std::vector<double>& wgh = it->second[2];
        size_t n = xsi.size();
        for (size_t i0 = 0; i0 < n; ++i0) {
            for (size_t j = 0; j < 3; ++j) {
                p0[j] = pa0[j] + xsi[i0]*db0[j] + eta[i0]*dc0[j];
            }
            icqQuadratureSetObserver(self, p0);
            for (size_t i1 = 0; i1 < n; ++i1) {
                for (size_t j = 0; j < 3; ++j) {
                    p1[j] = pa1[j] + xsi[i1]*db1[j] + eta[i1]*dc1[j];
                }
                res += wgh[i0] * (*self)->func->operator()(p1) * wgh[i1];
            }
        }
    }
    return 0.5 * area1 * res;
}

extern "C"
void icqQuadratureDel(icqQuadratureType **self) {
    
    delete (*self)->func;
    delete (*self);
}

