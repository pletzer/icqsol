#include <limits>
#include <cmath>
#include <icqInsideLocator.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <iostream>


icqInsideLocatorType::icqInsideLocatorType(vtkPolyData* pdata) {

    // A very small number
	this->eps = 1.98743174483817*std::numeric_limits<double>::epsilon(); 
	
    this->pdata = pdata;

    //  Boundaing box: (xmin,xmax, ymin,ymax, zmin,zmax)
    this->pdata->GetPoints()->ComputeBounds();
    double* bounds = this->pdata->GetPoints()->GetBounds();

    // Make the box a little bigger
    this->boxMin[0] = bounds[0] - this->eps;
    this->boxMax[0] = bounds[1] + this->eps;
    this->boxMin[1] = bounds[2] - this->eps;
    this->boxMax[1] = bounds[3] + this->eps;
    this->boxMin[2] = bounds[4] - this->eps;
    this->boxMax[2] = bounds[5] + this->eps;

    this->radius = 0;
    for (size_t k = 0; k < 3; ++k) {
        double boxLen = this->boxMax[k] - this->boxMin[k];
        this->center[k] = 0.5*(this->boxMin[k] + this->boxMax[k]);
        this->radius += boxLen*boxLen;
    }
    this->radius = 0.5*sqrt(this->radius);
}

int 
icqInsideLocatorType::isPointInside(const double* point) {

    // Quick check
    if (this->isPointInSphere(point) == 0 || 
        this->isPointInBox(point) == 0) {
        return ICQ_NO;
    }
    
    // Ray should be as short as possible
    this->setRayDirection(point);
    int numIntersections = 0;
    
    // point - pa
    double p[3];

    // pb - pa
    double b[3];

    // pc - pa
    double c[3];

    double areaVec[3];

    // The triangle's vertices
    double pa[3], pb[3], pc[3];

    // Parametric coordinates
    double xsi, eta, lam;
    
    // Iterate over the polygons
    vtkPoints* points = this->pdata->GetPoints();
    vtkCellArray* cells = this->pdata->GetPolys();
    vtkIdType numCells = cells->GetNumberOfCells();
    vtkIdList* ptIds = vtkIdList::New();
    cells->InitTraversal();
    for (vtkIdType i = 0; i < numCells; ++i) {
        cells->GetNextCell(ptIds);
        vtkIdType numPoints = ptIds->GetNumberOfIds();

        // Must have at three points to make a cell
        if (numPoints < 3) {
            continue;
        }

        vtkIdType ia = ptIds->GetId(0);
        points->GetPoint(ia, pa);
        double paDotRay = 0;
        for (size_t k = 0; k < 3; ++k) {
            p[k] = point[k] - pa[k];
            paDotRay += (pa[k] - point[k]) * this->rayDirection[k];
        }
        
        // Subdivide the polygon into triangles
        for (vtkIdType j = 1; j < numPoints - 1; ++j) {
            vtkIdType ib = ptIds->GetId(j);
            vtkIdType ic = ptIds->GetId(j + 1);
            points->GetPoint(ib, pb);
            points->GetPoint(ic, pc);
            double pbDotRay = 0;
            double pcDotRay = 0;
            for (size_t k = 0; k < 3; ++k) {
                b[k] = pb[k] - pa[k];
                c[k] = pc[k] - pa[k];
                pbDotRay += (pb[k] - point[k]) * this->rayDirection[k];
                pcDotRay += (pc[k] - point[k]) * this->rayDirection[k];
            }

            // Points cannot be degenerate
            areaVec[0] = b[1]*c[2] - b[2]*c[1];
            areaVec[1] = b[2]*c[0] - b[0]*c[2];
            areaVec[2] = b[0]*c[1] - b[1]*c[0];
            double area = 0;
            for (size_t k = 0; k < 3; ++k) {
                area += areaVec[k]*areaVec[k];
            }
            area = sqrt(area);
            if (fabs(area) < this->eps) {
                continue;
            }
            
            // At least one of the points must be in the direction of the
            // ray
            if (paDotRay > -this->eps || pbDotRay > -this->eps || pcDotRay > -this->eps) {
            	int res = this->rayIntersectsTriangle(p, b, c, &xsi, &eta, &lam);
                if (res == ICQ_YES) {
                    numIntersections++;
                }
                else if (res == ICQ_MAYBE) {

                    if (xsi > this->eps && xsi < 1. - this->eps &&
                        eta > this->eps && xsi + eta < 1. - this->eps && 
                        lam > -this->eps) {
                        // Likely point (nearly) on triangle
                        numIntersections++;
                    }
                }
        	}
        }
    }
    ptIds->Delete();
    
    int res = ICQ_NO;
    if (numIntersections % 2 == 1) {
        res = ICQ_YES;
    }
    
    return res;
}

int icqInsideLocatorType::isPointInSphere(const double* point) {

    int res = ICQ_NO;
    double radSqr = 0;
    for (size_t k = 0; k < 3; ++k) {
        radSqr += (point[k] - this->center[k])*(point[k] - this->center[k]);
    }
    if (radSqr < this->radius*this->radius) res = ICQ_YES;
    return res;
}

int icqInsideLocatorType::isPointInBox(const double* point) {

    for (size_t k = 0; k < 3; ++k) {
    	if (point[k] < this->boxMin[k]) return ICQ_NO;
    	if (point[k] > this->boxMax[k]) return ICQ_NO;
    }
    return ICQ_YES;
}

void icqInsideLocatorType::setRayDirection(const double* point) {

   // Shoot towards the box plane that is closest but avoid 
   // shooting along one of the main axes to minimize the risk
   // of hitting either an edge or a vertex.

   // Initialize to some small random directions
   this->rayDirection[0] = 0.0061246565456;
   this->rayDirection[1] = -0.0037655645623;
   this->rayDirection[2] = 0.0078962767621;

   size_t index = 0; 
   int sign = 1;
   double minD = std::numeric_limits<double>::max();

   for (size_t k = 0; k < 3; ++k) {

       double hi = this->boxMax[k] - point[k];
       double lo = point[k] - this->boxMin[k];
       double d = (lo < hi? lo: hi);
       if (d < minD) {
           index = k;
           sign = (lo < hi? -1: +1);
           minD = d;
       }
   }
   
   // Set the ray in this direction
   this->rayDirection[index] = sign;
}


int icqInsideLocatorType::rayIntersectsTriangle(const double* p,
                                                const double* b,
                                                const double* c,
                                                double* xsi, double* eta, double* lam) {

    const double* d = this->rayDirection;
    double det = b[2]*c[1]*d[0] - b[1]*c[2]*d[0] - b[2]*c[0]*d[1] + b[0]*c[2]*d[1] + b[1]*c[0]*d[2] - b[0]*c[1]*d[2];
    if (det == 0) {
        // face parallel to ray? This should never happen!!
        return ICQ_MAYBE;
    }

    *xsi = -c[2]*d[1]*p[0] + c[1]*d[2]*p[0] + c[2]*d[0]*p[1] - c[0]*d[2]*p[1] - c[1]*d[0]*p[2] + c[0]*d[1]*p[2];
    *xsi /= -det;
    *eta = b[2]*d[1]*p[0] - b[1]*d[2]*p[0] - b[2]*d[0]*p[1] + b[0]*d[2]*p[1] + b[1]*d[0]*p[2] - b[0]*d[1]*p[2];
    *eta /= -det;
    *lam = b[2]*c[1]*p[0] - b[1]*c[2]*p[0] - b[2]*c[0]*p[1] + b[0]*c[2]*p[1] + b[1]*c[0]*p[2] - b[0]*c[1]*p[2];
    *lam /= -det;

    if (*xsi < -this->eps ||
        *eta < -this->eps ||
        *xsi > 1.0 + this->eps ||
        *xsi + *eta > 1.0 + this->eps ||
        *lam < -this->eps) {
        // No intersection
        return ICQ_NO;
    }

    // A good chance to have an intersection
    if (*xsi > this->eps && *xsi < 1 - this->eps && 
        *eta > this->eps && *xsi + *eta < 1 - this->eps) {
        // Definitely an intersection
        return ICQ_YES;
    }

    // Cannot say for sure
    return ICQ_MAYBE;
}

// C interface

extern "C"
void icqInsideLocatorInit(icqInsideLocatorType** self, vtkPolyData* pdata) {
    *self = new icqInsideLocatorType(pdata);
}

extern "C"
void icqInsideLocatorDel(icqInsideLocatorType** self) {
    delete *self;
}
  
extern "C"  
int icqInsideLocatorIsPointInside(icqInsideLocatorType **self, const double* point) {
    return (*self)->isPointInside(point);
}

