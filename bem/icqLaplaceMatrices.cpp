#include <vector>
#include <cmath>
#include <icqQuadrature.hpp>
#include <icqLaplaceMatrices.hpp>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>

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

    // Edges for the source triangle
    std::vector<double> pb2Src(3);
    std::vector<double> pc2Src(3);

    // Source area vector
    std::vector<double> areaSrcVec(3);

    // Observer triangle vertices
    std::vector<double> paObs(3);
    std::vector<double> pbObs(3);
    std::vector<double> pcObs(3);

    // Triangle mid-points
    std::vector<double> xObs(3);
    std::vector<double> xSrc(3);

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

        // The triangle's normal vector, mid point, and area at the center of the triangle
        for (size_t j = 0; j < 3; ++j) {
            pb2Src[j] = pbSrc[j] - paSrc[j];
            pc2Src[j] = pcSrc[j] - paSrc[j];
            xSrc[j] = (paSrc[j] + pbSrc[j] + pcSrc[j])/3.;
        }
        for (size_t j = 0; j < 3; ++j) {
            size_t j1 = (j + 1) % 3;
            size_t j2 = (j + 2) % 3;
            areaSrcVec[j] = pb2Src[j1]*pc2Src[j2] - pb2Src[j2]*pc2Src[j1];
        }
        double areaSrc = 0;
        for (size_t j = 0; j < 3; ++j) {
            areaSrc += areaSrcVec[j]*areaSrcVec[j];
        }
        areaSrc = sqrt(areaSrc);

        // Iterate over the observer triangles
        for (size_t iObs = 0; iObs < numTriangles; ++iObs) {

            if (iObs == jSrc) {
                // Singular term, treated elsewhere
                continue;
            }

            points->GetPoint(cells[3*iObs + 0], &paObs[0]);
            points->GetPoint(cells[3*iObs + 1], &pbObs[0]);
            points->GetPoint(cells[3*iObs + 2], &pcObs[0]);

            // Observer is at mid point
            for (size_t j = 0; j < 3; ++j) {
                xObs[j] = (paObs[j] + pbObs[j] + pcObs[j])/3.;
            }

            icqQuadratureSetObserver(&self, &xObs[0]);

            gMat[numTriangles*iObs + jSrc] = icqQuadratureEvaluate(&self, 
                                                                   maxOrder, 
                                                                   &paSrc[0], &pbSrc[0], &pcSrc[0]);
        }

    }

    icqQuadratureDel(&self);

}
