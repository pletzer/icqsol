#include <vector>
#include <cmath>
#include <icqQuadrature.hpp>
#include <icqLaplaceMatrices.hpp>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>

double getArea(const std::vector<double>& pa, 
               const std::vector<double>& pb,
               const std::vector<double>& pc) {

    std::vector<double> db(3);
    std::vector<double> dc(3);
    for (size_t j = 0; j < 3; ++j) {
        db[j] = pb[j] - pa[j];
        dc[j] = pc[j] - pa[j];
    }

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

/**
 * Compute the off diagonal influence matrix elements
 * @param pdata vtkPolyData instance
 * @param gMat influence matrix
 */
extern "C"
void computeOffDiagonalTerms(vtkPolyData* pdata, double* gMat) {

    icqQuadratureType* self;
    icqQuadratureInit(&self);
    int maxOrder = icqQuadratureGetMaxOrder(&self);

    // Vertex coordinates for the source triangle
    std::vector<double> paSrc(3);
    std::vector<double> pbSrc(3);
    std::vector<double> pcSrc(3);

    // Observer triangle vertices
    std::vector<double> paObs(3);
    std::vector<double> pbObs(3);
    std::vector<double> pcObs(3);

    vtkCellArray* polys = pdata->GetPolys();
    vtkPoints* points = pdata->GetPoints();
    size_t numTriangles = polys->GetNumberOfCells();
    polys->InitTraversal();

    // Build the connectivity array, all triangles
    vtkIdList* ptIds = vtkIdList::New();
    std::vector<int> cells(numTriangles*3);
    for (size_t i = 0; i < numTriangles; ++i) {
        polys->GetNextCell(ptIds);
        size_t ia = ptIds->GetId(0);
        size_t ib = ptIds->GetId(1);
        size_t ic = ptIds->GetId(2);
        cells[3*i + 0] = ia;
        cells[3*i + 1] = ib;
        cells[3*i + 2] = ic; 
    }
    ptIds->Delete();

    // Iterate over the source triangles
    #pragma omp parallel for
    for (size_t jSrc = 0; jSrc < numTriangles; ++jSrc) {

        size_t ia = cells[3*jSrc + 0];
        size_t ib = cells[3*jSrc + 1];
        size_t ic = cells[3*jSrc + 2];

        points->GetPoint(ia, &paSrc[0]);
        points->GetPoint(ib, &pbSrc[0]);
        points->GetPoint(ic, &pcSrc[0]);
        double areaSrc = getArea(paSrc, pbSrc, pcSrc);

        // Iterate over the observer triangles
        for (size_t iObs = jSrc + 1; iObs < numTriangles; ++iObs) {

            points->GetPoint(cells[3*iObs + 0], &paObs[0]);
            points->GetPoint(cells[3*iObs + 1], &pbObs[0]);
            points->GetPoint(cells[3*iObs + 2], &pcObs[0]);
            double areaObs = getArea(paObs, pbObs, pcObs);

            double g = icqQuadratureEvaluateDouble(&self,  maxOrder,
                                                   &paObs[0], &pbObs[0], &pcObs[0],
                                                   &paSrc[0], &pbSrc[0], &pcSrc[0]);
            gMat[numTriangles*iObs + jSrc] = g;
            gMat[numTriangles*jSrc + iObs] = g * areaObs / areaSrc;
        }

    }

    icqQuadratureDel(&self);

}
