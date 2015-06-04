#ifndef DISIMPLV_H
#define DISIMPLV_H_ 
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h> 
#include <limits>

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "gkls.h"
#include "rnd_gen.h"

using namespace std;


class Point {
    Point(const Point& other){}
    Point& operator=(const Point& other){};
public:
    Point(){};
    Point(int *c, int argc){
        for (int i=0; i<argc ; i++){
            _X.push_back(double(c[i]));
        };
    };
    Point(double *c, int argc){
        for (int i=0; i<argc ; i++){
            _X.push_back(c[i]);
        };
    };
    Point(double c1){};
    Point(double c1, double c2){
        _X.push_back(c1);
        _X.push_back(c2);
    };

    vector<double> _X;  // Coordinates in normalised [0,1]^n space  
    vector<double> _values;
         
    void add_coord(double x) {
        _X.push_back(x);
    };
    void add_value(double value) {
        _values.push_back(value);
    };

    // double* to_double(point){    // Todo: use double x[_D] for _X instead of vector. Won't need this method when using GKLS function generator
    //     for (int i=0; i < _X.size() ; i++) {
    //         point[i] = _X[i];
    //     };
    //     return point;
    // };

    int size(){
        return _X.size();
    };

    void print(){
        cout << "       ";
        for (int i=0; i < _X.size(); i++){
            cout << _X[i] << "  \t";
        };
        for (int i=0; i < _values.size(); i++){
            if (i == 0) { cout << "->\t"; };
            cout << _values[i] << "  ";
        };
        cout << endl;
    };

    virtual ~Point(){
        vector<double>::iterator cit = _X.begin();
        while (cit != _X.end()) {
            cit = _X.erase(cit);    
        };
        vector<double>::iterator vit = _values.begin();
        while (vit != _values.end()) {
            vit = _values.erase(vit);    
        };
    };
};

class Points {  // Binary balancing tree or simply linked list
    Points(const Points& other){}
    Points& operator=(const Points& other){}
public:                
    Points(){};
    vector<Point*> _points;
   
    void add(Point* point){
        _points.push_back(point);
    };

    Point* get(double *c, int argc){
        for (int i=0; i < _points.size(); i++){
            bool matches = true;
            for (int j = 0; j < argc; j++){
                if (_points[i]->_X[j] != c[j]) {
                    matches = false;
                    break;
                };
            };
            if (matches) {
                return _points[i];
            };
        };
        return 0;
    };

    Point* get(Point* p){
        for (int i=0; i < _points.size(); i++){
            bool matches = true;
            for (int j = 0; j < p->size(); j++){
                if (_points[i]->_X[j] != p->_X[j]) {
                    matches = false;
                    break;
                };
            };
            if (matches) {
                return _points[i];
            };
        };
        return 0;
    };

    // Point* get_point(double c1, double c2){
    //     for (int i=0; i < _points.size(); i++){
    //         if (_points[i]->_X[0] == c1 && _points[i]->_X[1] == c2){
    //             return _points[i];
    //         };
    //     };
    //     return 0;
    // };

    int size() {
        return _points.size();
    };

    void print(){
        cout << "Points: " << endl;
        for (int i=0; i < _points.size(); i++) {
            _points[i]->print();
        };
    };

    virtual ~Points(){
        for (int i=0; i < _points.size(); i++) {
            delete _points[i];
        };
        _points.clear();
    };
};


class Function {
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(){
        _calls = 0;
        _f_min = numeric_limits<int>::max();
    };
    string _name;
    int _D;
    Point* _lb;
    Point* _ub;
    Point* _glob_x;  // Point where global function minimum is (should be list)
    double _glob_f;  // Predefined global function minimum
    double _L;       // Lipschitz constant

    int _calls;
    double _f_min;
    Point* _x_min;
    Points _points;


    Point* get(double *c, int argc){
        Point* cached_point = _points.get(c, argc);
        if (cached_point) {
            return cached_point;
        } else {
            Point* p = new Point(c, argc);
            double val = value(p);
            if (_f_min > val) {
                _f_min = val;
                _x_min = p;
            };
            p->add_value(val);
            _calls += 1;
            _points.add(p);
            return p;
        };
    };

    Point* get(Point* p){
        Point* cached_point = _points.get(p);
        if (cached_point) {
            return cached_point;
        } else {
            double val = value(p);
            if (_f_min > val) {
                _f_min = val;
                _x_min = p;
            };
            p->add_value(val);
            _calls += 1;
            _points.add(p);
            return p;
        };
    };

    void print(){
        cout << "\nCalls: " << _calls << "   f_min: " << _f_min << endl;
    };

    double transform(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    double pe(){
        if (_glob_f != 0) {
            return (_f_min - _glob_f) / fabs(_glob_f) * 100.; 
        };
        return _f_min * 100.; 
    };

    virtual double value(Point* point) = 0;

    virtual ~Function(){};
};

class Branin : public Function {
    Branin(const Branin& other){};
    Branin& operator=(const Branin& other){};
public:
    Branin(): Function(){
        _name = "Branin";
        _D = 2;
        _lb = new Point(-5, 0);
        _ub = new Point(10, 15);
        _glob_x = new Point(M_PI, 12.275);  // Point where global function minimum is (should be list)
        _glob_f = 0.397887;          // Predefined global function minimum
        _L = 109.94813585;
    };

    double value(Point* point) {
        double x1 = transform(point, 0); 
        double x2 = transform(point, 1); 
        double part1 = pow((x2 - 5./(4*pow(M_PI, 2))*pow(x1,2) + 5./M_PI*x1 -6), 2);
        double part2 = 10.*(1. - 1./(8*M_PI))*cos(x1) + 10.;
        return part1 + part2;
    };
};


class GKLSFunction: public Function {
    GKLSFunction(const GKLSFunction& other){};
    GKLSFunction& operator=(const GKLSFunction& other){};
public:
    GKLSFunction(int cls, int function_id): Function() {
        int _GKLS_class_D[] = {2, 2, 3, 3, 4, 4, 5, 5};
        double _GKLS_class_global_dists[] = {0.9, 0.9, 0.66, 0.9, 0.66, 0.9, 0.66, 0.66};
        double _GKLS_class_global_radiuses[] = {0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2};

        _name = "GKLSFunction";
        _D = _GKLS_class_D[cls];
        _global_dist = _GKLS_class_global_dists[cls];
        _global_radius = _GKLS_class_global_radiuses[cls]; 

        _lb = new Point();
        _ub = new Point();
        for (int i=0; i<_D; i++) {
            _lb->add_coord(-1.);
            _ub->add_coord(1.);
        };

        // _glob_x = new Point();  // Point where global function minimum is (should be list)
        _glob_f = -1.;          // Predefined global function minimum
        // _L = ;

        assert(GKLS_set_default()== GKLS_OK);     // Standartiniai nustatymai

        GKLS_dim = _D;
        GKLS_global_dist = _global_dist;
        GKLS_global_radius = _global_radius;
        GKLS_num_minima = 10;
        GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
        for (unsigned int i = 0; i < GKLS_dim; i++) {
            GKLS_domain_left[i] = -1;
            GKLS_domain_right[i] = 1;
        };
        assert(GKLS_domain_alloc() == GKLS_OK);
        assert(GKLS_parameters_check() == GKLS_OK);
        assert(GKLS_arg_generate(function_id) == GKLS_OK);
        int n = GKLS_glob.num_global_minima;
        assert(n == 1);
        int glob_idx = 1;
        assert(GKLS_minima.f[glob_idx] == GKLS_global_value);
    };

    double _global_dist;
    double _global_radius;


    double value(Point* point) {
        // for (unsigned int i = 0; i < GKLS_dim; i++) {
        //     cout << GKLS_minima.local_min[glob_idx][i] << " ";
        // }
        // cout << endl;

        double p[_D];
        for (int i; i < _D; i++) {
            p[i] = point->_X[i];
        };
        // Convert Point to double point[_D];
        return GKLS_D_func(p);
    };

    // GKLS_free();
    // GKLS_domain_free();
};




#endif
