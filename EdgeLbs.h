#ifndef EDGELBS_H
#define EDGELBS_H
/* Approximate min_lb_point by averaging weighted verts (good for Euclidean space */
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include <sys/time.h>
#include <vector>
#include "utils.h"
#include "subsimplex.h"
#include <math.h>
// Depencencies:   main  <  Algorithm  <  Simplex  <  EdgeLbs  <  Subsimplex


using namespace std;

class EdgeLbs {
    EdgeLbs(const EdgeLbs& other) {};
    EdgeLbs& operator=(const EdgeLbs& other) {};
public:
    EdgeLbs(vector<Point*> verts, vector<double> Ls, int crit_id) {
        // Note: Should pass simplex here and use store dist between verts in it
        _verts = verts;
        _L = Ls[crit_id];
        _D = _verts.size() - 1;
        _V = _verts.size();
        _crit_id = crit_id;
    };
    vector<Point*> _verts;   // Simplex vertexes
    double _L;
    int _D;
    int _V;
    int _crit_id;

    double get_lb_value(double* p) {
        /* Returns lb value at a given point */
        double cone_value;
        double dist;
        double lb_value = -numeric_limits<double>::max();;

        for (int i=0; i < _verts.size(); i++) {
            // dist = l2norm(point, verts[i]);
            dist = 0;
            //// L2norm:
            for (int j=0; j < _verts[i]->size(); j++){
                dist += pow(_verts[i]->_X[j] - p[j], 2);
            };
            dist = sqrt(dist);
            //// L1norm:
            // for (int j=0; j < _verts[i]->size(); j++){
            //     dist += fabs(_verts[i]->_X[j] - p[j]);
            // };

            cone_value = _verts[i]->_values[_crit_id] - _L*dist;
            if (lb_value < cone_value) {
                lb_value = cone_value;
            };
        };
        return lb_value;
    };

    double get_inaccurate_compromise(Point* vert, double* X, double val2) {  // Minimize edge lb
        double val1 = vert->_values[0];
        double* X1 = vert->_X;
        double X2[_D];
        for (int k=0; k < _D; k++) {
            X2[k] = X[k];
        };
        double dist = l2norm(X1, X2, _D);

        double improvement = 1.;
        while (improvement > 1e-4) {
            if (dist == 0) { return val2; };
            // Find point for Lipschitz minimum
            double c = 0.5 + (val2 - val1) / (2 * _L * dist);

            double X[_D];
            for (int k=0; k < _D; k++) {
                X[k] = c * X1[k];
                X[k] += (1. - c) * X2[k];
            };
            double lb_value = get_lb_value(X);

            if (lb_value < val2) {
                improvement = l2norm(X2, X, _D);
                val2 = lb_value;
                for (int k=0; k < _D; k++) {
                    X2[k] = X[k];
                };
            } else {
                improvement = 0;
            };

            dist = l2norm(X1, X2, _D);
        };
        return val2;
    };

    Point* minimize() {       // Returns estimated simplex_lb_min_value
        double X[_D];
        double value = numeric_limits<double>::max();

        sort(_verts.begin(), _verts.end(), Point::compare_by_value);

        // combinatorically select vert pairs
        double dist_ab, ca, cb;
        double pairs_count = 0;
        for (int a=0; a <= _D; a++) {
            for (int b=a; b <= _D; b++) {
                double lb_value = numeric_limits<double>::max();
                if (a != b) {
                    double X[_D];
                    pairs_count += 1;
                    dist_ab = l2norm(_verts[a], _verts[b]);    // Distance in variable space

                    ca = 0.5 + (_verts[b]->_values[_crit_id] - _verts[a]->_values[_crit_id]) / (2 * _L * dist_ab);

                    for (int k=0; k < _D; k++) {
                        X[k] = ca * _verts[a]->_X[k];
                        X[k] += (1. - ca) * _verts[b]->_X[k];
                    };

                    double expected_lb_value = _verts[b]->_values[_crit_id] - _L * dist_ab * ca;

                    lb_value = get_lb_value(X);

                    if (lb_value != expected_lb_value) {    // Note: should compare with precision?
                        double lb_value1 = get_inaccurate_compromise(_verts[a], X, lb_value);
                        double lb_value2 = get_inaccurate_compromise(_verts[b], X, lb_value);

                        if (lb_value1 < lb_value2) {
                            lb_value = lb_value1;
                        } else {
                            lb_value = lb_value2;
                        };
                    };
                };
                if (lb_value < value) {
                    value = lb_value;
                };
            };
        };

        Point* result = new Point(X, _D);
        result->add_value(value);
        return result;
    };

    virtual ~EdgeLbs() {
    };
};

#endif
