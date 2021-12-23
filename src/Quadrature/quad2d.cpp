// Creates quadrature rule in a cell
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//


#include <quad2d.hpp>
#include <iostream>
#include <mesh.hpp>
using namespace HArDCore2D;
QuadRuleTriangle::QuadRuleTriangle(size_t doe, bool warn)
    : _npts(dunavant_order_num(std::max(int(doe),1))),
    _xy(0),
    _w(0),
    _xyphys(0),
    area(0){
    if (doe > max_doe && warn) {
        std::cerr << "Warning: quadrature rule of degree " << doe
                  << " requested, but the maximum available is degree "
                  << max_doe << std::endl;
    }

    _w = new double[_npts];
    _xy = new double[_npts * 2];
    _xyphys = new double[_npts * 2];
    //    Input, int RULE, the index of the rule.
    //    Input, int ORDER_NUM, the order (number of points) of the rule.
    //    Output, double XY[2*ORDER_NUM], the points of the rule.
    //    Output, double W[ORDER_NUM], the weights of the rule.

    dunavant_rule(std::max(int(doe),1), _npts, _xy, _w);
}

QuadRuleTriangle::~QuadRuleTriangle(){
    delete[] _xy;
    delete[] _w;
    delete[] _xyphys;
}

size_t QuadRuleTriangle::nq() { return _npts; }
double QuadRuleTriangle::xq(size_t i) {
    if ((i * 2 + 1) >= _npts * 2) {
        throw "QuadRuleTriangle: trying to access array element out of bounds";
    }
    return _xyphys[i * 2];
}
double QuadRuleTriangle::yq(size_t i) {
    if ((i * 2 + 1) >= _npts * 2) {
        throw "QuadRuleTriangle: trying to access array element out of bounds";
    }
    return _xyphys[i * 2 + 1];
}
double QuadRuleTriangle::wq(size_t i) {
    if ((i) >= _npts) {
        throw "QuadRuleTriangle: trying to access array element out of bounds";
    }
    return _w[i] * std::abs(area);
}
void QuadRuleTriangle::setup(double xV[], double yV[]) {
    double t[6];
    for (int i = 0; i <= 2; i++) {
        t[i * 2] = xV[i];
        t[i * 2 + 1] = yV[i];
    }
    //    Input, double T[2*3], the coordinates of the vertices.
    //    The vertices are assumed to be the images of (0,0), (1,0) and
    //    (0,1) respectively.
    //    Input, int N, the number of objects to transform.
    //    Input, double REF[2*N], points in the reference triangle.
    //    Output, double PHY[2*N], corresponding points in the
    //    physical triangle.
    reference_to_physical_t3(t, _npts, _xy, _xyphys);
    area = triangle_area(t);
};
