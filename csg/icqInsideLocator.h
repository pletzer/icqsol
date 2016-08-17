#ifndef ICQ_INSIDE_LOCATOR
#define ICQ_INSIDE_LOCATOR

#include <vtkPolyData.h>

enum {ICQ_NO, ICQ_YES, ICQ_MAYBE};

class icqInsideLocatorType {

public:

/**
 * Constructor
 * @param pdata instance of vtkPolyData
 */

    icqInsideLocatorType(vtkPolyData* pdata);
/**
 * Destructor
 */
    ~icqInsideLocatorType(){}

/**
 * Check if a point is inside the shape
 * @param point point
 * @return ICQ_YES, ICQ_NO, or ICQ_MAYBE
 */
    int isPointInside(const double* point);

private:

    double radius;
    double eps;
    vtkPolyData* pdata;
    double boxMin[3];
    double boxMax[3];
    double center[3];
    double rayDirection[3];

    int isPointInBox(const double* point);
    int isPointInSphere(const double* point);

    int rayIntersectsTriangle(const double* p, 
                              const double* b,
                              const double* c,
                              double* xsi, double* eta, double* lam);

    void setRayDirection(const double* point);
    
};

// C interface
extern "C" {
    void icqInsideLocatorInit(icqInsideLocatorType** self, vtkPolyData* pdata);
    void icqInsideLocatorDel(icqInsideLocatorType** self);
    int icqInsideLocatorIsPointInside(icqInsideLocatorType **self, const double* point);
}

#endif // ICQ_POINT_INSIDE_VTK_POLY_DATA