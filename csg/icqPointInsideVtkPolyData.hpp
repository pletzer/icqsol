#ifndef ICQ_POINT_INSIDE_VTK_POLY_DATA
#define ICQ_POINT_INSIDE_VTK_POLY_DATA

#include <vtkPolyData.h>
#include <vector>

enum {ICQ_NO, ICQ_YES, ICQ_MAYBE};

class icqPointInsideVtkPolyData {

public:

    icqPointInsideVtkPolyData(vtkPolyData* pdata);

    ~icqPointInsideVtkPolyData(){}

    int isInShape(const double* point);

private:

    double radius;
    double eps;
    vtkPolyData* pdata;
    std::vector<double> epsOffset;
    std::vector<double> boxMin;
    std::vector<double> boxMax;
    std::vector<double> center;
    std::vector<double> rayDirection;

    int isInBox(const double* point);
    int isInSphere(const double* point);

    int rayIntersectsTriangle(const std::vector<double>& p, 
                              const std::vector<double>& b,
                              const std::vector<double>& c);

    void setRayDirection(const double* point);
    
};

#endif // ICQ_POINT_INSIDE_VTK_POLY_DATA