#ifndef FUNCTIONS_H
#define FUNCTIONS_H 
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h> 
#include <limits>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include "gkls.h"
#include "rnd_gen.h"

using namespace std;

class Simplex;

class Point {
    Point(const Point& other){}
    Point& operator=(const Point& other){};
public:
    Point(int D){
        _D = D;
        _X = (double*) malloc((D)*sizeof(double));
    };
    Point(int *c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = double(c[i]);
        };
    };
    Point(double *c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = c[i];
        };
    };
    Point(double c, int argc){
        _D = argc;
        _X = (double*) malloc((argc)*sizeof(double));
        for (int i=0; i<argc ; i++){
            _X[i] = c;
        };
    };
    Point(double c1, double c2){
        _D = 2;
        _X = (double*) malloc(2*sizeof(double));
        _X[0] = c1;
        _X[1] = c2;
    };

    int _D;
    double* _X;  // Coordinates in normalised [0,1]^n space  
    vector<double> _values;
    vector<Simplex*> _simplexes;  // Simplexes, which have this point as vertex
         
    void add_value(double value) {
        _values.push_back(value);
    };

    Point* copy() {
        Point* point_copy = new Point(_X, _D);
        for (int i=0; i < _values.size(); i++) {
            point_copy->add_value(_values[i]);
        };
        return point_copy;
    };

    int size(){
        return _D;
    };

    void _neighbours_estimates_should_be_updated();
    // {
    //     for (int sid=0; sid < _simplexes.size(); sid++) {
    //         _simplexes[sid]->_should_estimates_be_updated = true;
    //     };
    // };

    void print(){
        // cout.precision(17);
        cout << "       ";
        for (int i=0; i < size(); i++){
            // cout << fixed << _X[i] << "  \t";
            cout << _X[i] << "  \t";
        };
        for (int i=0; i < _values.size(); i++){
            if (i == 0) { cout << "->\t"; };
            // cout << fixed << _values[i] << "  ";
            cout << _values[i] << "  ";
        };
        cout << endl;
    };

    static bool compare_by_value(Point* p1, Point* p2) {
        return p1->_values[0] < p2->_values[0];
    };

    friend ostream& operator<<(ostream& o, const Point& p){
        for (int i=0; i < p._D; i++) {
            o << p._X[i];
            if (i != p._D - 1) { o << ","; };
        };
        return o;
    };

    virtual ~Point(){
        free(_X);
        vector<double>::iterator vit = _values.begin();
        while (vit != _values.end()) {
            vit = _values.erase(vit);    
        };
    };
};


class PointTree;

class PointTreeNode {  // Binary balancing tree or simply linked list
    PointTreeNode(const PointTreeNode& other){}
    PointTreeNode& operator=(const PointTreeNode& other){}
public:                
    PointTreeNode(double value=numeric_limits<double>::max()){
        _height = 1;
        _value = value;
        _parent = 0;
        _left = 0;
        _right = 0;
        _subtree = 0;
        _point = 0;
    };
    int _height;
    double _value;
    PointTreeNode* _parent;
    PointTreeNode* _left;
    PointTreeNode* _right;
    PointTree* _subtree;     // Next dimension head
    Point* _point;           // Only last dimension node will have _point != 0;

    void print(){
        if (_left != 0) { cout << "l"; _left->print(); };
        cout << _value << "("<< _height << ")";
        if (_right != 0) { cout << "r"; _right->print(); };
    };

    virtual ~PointTreeNode(); // {
    //     if (_left != 0) { delete _left; };
    //     if (_right != 0) { delete _right; };
    //     if (_point != 0) { delete _point; };
    //     if (_subtree != 0) { delete _subtree; };
    // };
};

class PointTree { // Head of the tree
    PointTree(const PointTree& other){}
    PointTree& operator=(const PointTree& other){}
public:
    PointTree(){
         _tree_root = 0;
         _dim = 1;
    };
    PointTree(int dim){
         _tree_root = 0;
         _dim = dim;
    };
    PointTreeNode* _tree_root;
    int _dim;

    void update_height(PointTreeNode* node) {
        int lh = 0;
        int rh = 0;
        if (node->_left != 0) { lh = node->_left->_height; }; 
        if (node->_right != 0) { rh = node->_right->_height; };
        if (lh > rh) {
            node->_height = lh + 1;
        } else {
            node->_height = rh + 1;
        };
        // Also update all ancestors heights
        if (node->_parent != 0) {
            update_height(node->_parent);
        };
    };

    void left_right_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched_node;
        // node left right  <-  node left right left 
        diatteched_node = node->_left->_right;
        node->_left->_right = node->_left->_right->_left;
        if (node->_left->_right != 0) { node->_left->_right->_parent = node->_left; };
        // Diatteched left = node->_left
        diatteched_node->_left = node->_left;
        node->_left->_parent = diatteched_node;
        // node left  <-  node left right
        node->_left = diatteched_node;
        diatteched_node->_parent = node;
        // Update heights
        update_height(node);
        update_height(diatteched_node);
        update_height(diatteched_node->_left);
    };
    void left_left_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched;
        diatteched = node->_left;
        node->_left = node->_left->_right;
        if (node->_left != 0) { node->_left->_parent = node; };  
        diatteched->_parent = node->_parent;
        if (node->_parent != 0) {
            if (node->_parent->_left == node) {
                node->_parent->_left = diatteched;
            } else {
                node->_parent->_right = diatteched;
            };
        } else {
            _tree_root = diatteched;
        };
        diatteched->_right = node;
        node->_parent = diatteched;
        // Update heights
        update_height(node);
        update_height(diatteched);
    };
    void right_left_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched_node;
        // node left right  <-  node left right left 
        diatteched_node = node->_right->_left;
        node->_right->_left = node->_right->_left->_right;
        if (node->_right->_left != 0) { node->_right->_left->_parent = node->_right; };
        // Diatteched left = node->_left
        diatteched_node->_right = node->_right;
        node->_right->_parent = diatteched_node;
        // node left  <-  node left right
        node->_right = diatteched_node;
        diatteched_node->_parent = node;
        // Update heights
        update_height(node);
        update_height(diatteched_node);
        update_height(diatteched_node->_right);
    };
    void right_right_rebalance(PointTreeNode* node) {
        PointTreeNode* diatteched;
        diatteched = node->_right;
        node->_right = node->_right->_left;
        if (node->_right != 0) { node->_right->_parent = node; };
        diatteched->_parent = node->_parent;                       
        if (node->_parent != 0) {
            if (node->_parent->_left == node) {
                node->_parent->_left = diatteched;
            } else {
                node->_parent->_right = diatteched;
            };
        } else {
            _tree_root = diatteched;
        };
        diatteched->_left = node;
        node->_parent = diatteched;
        // Update heights
        update_height(node);
        update_height(diatteched);
    };

    void check_if_balanced(PointTreeNode* node) {
        int lh = 0;
        int rh = 0;
        int llh = 0;
        int lrh = 0;
        int rlh = 0;
        int rrh = 0;
        if (node->_left != 0) {
            lh = node->_left->_height;
            if (node->_left->_left != 0) { llh = node->_left->_left->_height; };
            if (node->_left->_right != 0) { lrh = node->_left->_right->_height; };
        };
        if (node->_right != 0) {
            rh = node->_right->_height;
            if (node->_right->_left != 0) { rlh = node->_right->_left->_height; };
            if (node->_right->_right != 0) { rrh = node->_right->_right->_height; };
        };
        if (abs(rh - lh) > 1) {
            // Not balanced, so rebalance
            if (rh > lh) {
                if (rrh > rlh) {
                    right_right_rebalance(node);
                } else {
                    right_left_rebalance(node);
                    right_right_rebalance(node);
                };
            };
            if (rh < lh) {
                if (llh > lrh) {
                    left_left_rebalance(node);
                } else {
                    left_right_rebalance(node);
                    left_left_rebalance(node);
                };
            };
        };
        if (node->_parent != 0) {
            check_if_balanced(node->_parent);
        };
    };

    Point* process_next_dimension(PointTreeNode* node, Point* point){
        // Creates next dimension tree if needed and adds point to it
        // Its last dimension
        if (point->_D == _dim) {  // Don't need next dimension
            if (node->_point == 0) {    // Save or return the point
                node->_point = point;
                return 0;
            } else {
                return node->_point;
            };
        };
        // Its not last dimension
        if (node->_subtree == 0) {  // Create subtree if it doesn't already exist
            node->_subtree = new PointTree(_dim + 1);
        };
        Point* found_point = node->_subtree->add(point);  // Get or insert point to the subtree
        if (found_point != 0) {  // We got point so return it 
            return found_point;
        } else {
            return 0;  // We inserted point
        };
    };

    Point* add(Point* point){
        // Get same point or insert given (if inserted returns 0)
        PointTreeNode* node = _tree_root;
        double value = point->_X[_dim -1];
        if (_tree_root == 0) {  // Create first tree node
            _tree_root = new PointTreeNode(value);
            node = _tree_root;
            process_next_dimension(node, point);
        } else {
            while (true) {  // Walk through tree
                if (value > node->_value) {
                    if (node->_right == 0) {
                        node->_right = new PointTreeNode(value);
                        node->_right->_parent = node;
                        update_height(node->_right);
                        process_next_dimension(node->_right, point);
                        check_if_balanced(node->_right);
                        return 0;
                    };
                    node = node->_right;
                } else if (value < node->_value) {
                    if (node->_left == 0) {
                        node->_left = new PointTreeNode(value);
                        node->_left->_parent = node;
                        update_height(node->_left);
                        process_next_dimension(node->_left, point);
                        check_if_balanced(node->_left);
                        return 0;
                    };
                    node = node->_left;
                } else {
                    // Node value matches given point value, move to next dimension.
                    return process_next_dimension(node, point);
                };
            };
        };
    };
    // setas, hsetas, steady set

    void print(){
        _tree_root->print();
        cout << endl;
    };

    virtual ~PointTree(){
        delete _tree_root;
    };
};

PointTreeNode::~PointTreeNode() {
    if (_left != 0) { delete _left; };
    if (_right != 0) { delete _right; };
    if (_point != 0) { delete _point; };
    if (_subtree != 0) { delete _subtree; };
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
        _points = new PointTree(); // new Points();
        _stopping_criteria = "x_dist";
    };
    string _name;
    string _stopping_criteria;
    int _D;
    Point* _lb;
    Point* _ub;
    Point* _glob_x;  // Point where global function minimum is (should be list)
    double _delta;   // Accuracy for stoping criteria based on distance from glob_x
    double _glob_f;  // Predefined global function minimum
    double _L;       // Global Lipschitz constant

    static bool _global_mem_allocated;

    int _calls;
    double _f_min;  // Best known function value
    Point* _x_min;  // Point where best known function value is
    double _distance_to_glob_x;  // Infinity distance to _glob_x from nearest known point
    Point* _x_nearest_to_glob_x;  // Nearest to _glob_x known point (infinity distance) 
    PointTree* _points;  // Points* _points;

    double get_distance_to_glob_x(Point* p) {
        double max_distance = 0;
        double dist = 0;
        for (int i=0; i<_D; i++) {
            dist = fabs(p->_X[i] - _glob_x->_X[i]);
            if (dist > max_distance) {
                max_distance = dist;
            };
        };
        return max_distance;
    };

    void update_meta(Point* p) {
        double val = value(p);
        if (_f_min > val) {
            _f_min = val;
            _x_min = p;
        };
        double distance_to_glob_x = get_distance_to_glob_x(p);
        if (distance_to_glob_x < _distance_to_glob_x) {
            _distance_to_glob_x = distance_to_glob_x;
            _x_nearest_to_glob_x = p;
        };

        p->add_value(val);
        _calls += 1;
    };

    void log_evaluation(Point* p) {
       ofstream log_file; 
       log_file.open("log/evaluations.txt", ios::app);
       for (int i=0; i < p->size(); i++) {
            log_file << setprecision(12) << p->_X[i] << " ";
       };
       log_file << " -> ";
       for (int i=0; i < p->_values.size(); i++) {
           log_file << setprecision(12) << p->_values[i] << " ";
       };
       log_file << endl;
       log_file.close();
    };

    Point* get(double *c, int argc){
        Point* p = new Point(c, argc);
        Point* cached_point = _points->add(p);
        if (cached_point) {
            delete p;
            log_evaluation(cached_point);
            return cached_point;
        } else {
            update_meta(p);
            log_evaluation(p);
            return p;
        };
    };

    Point* get(Point* p){  // Get Point with its function value
        Point* cached_point = _points->add(p);
        if (cached_point) {
            log_evaluation(cached_point);
            return cached_point;
        } else {
            update_meta(p);
            log_evaluation(p);
            return p;
        };
    };

    void print(){
        cout << _name << "  calls: " << _calls << "   f_min: " << _f_min << endl;
    };

    double transform(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    double transform_back(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return (point->_X[i] - _lb->_X[i]) / (_ub->_X[i]-_lb->_X[i]);
    };

    // double pe(){
    //     if (_glob_f != 0) {
    //         return (_f_min - _glob_f) / fabs(_glob_f) * 100.;
    //     };
    //     return _f_min * 100.;
    // };

    bool is_accurate_enougth() {
        double e;
        if (_stopping_criteria == "x_dist") {
            for (int i=0; i<_D; i++) {
                // if (_delta * (_ub->_X[i] - _lb->_X[i]) < fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i])) {   // Infinity norm
                //     return false;
                // };
                if (_delta * (1.0 - 0.0) < fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i])) {   // Infinity norm
                    return false;
                };
            };

            // for (int i=0; i<_D; i++) {
            //     cout << _x_nearest_to_glob_x->_X[i] << "   -  " << _glob_x->_X[i] << "  =  " << _x_nearest_to_glob_x->_X[i] - _glob_x->_X[i] << endl;
            //     // cout  << _delta * (_ub->_X[i] - _lb->_X[i]) << "   <   "<< fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i]) << endl;
            //     // if (_delta * (_ub->_X[i] - _lb->_X[i]) < fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i])) {   // Infinity norm
            //     //     return false;
            //     // };
            // };
            // cout << endl;
            // for (int i=0; i<_D; i++) {
            //     cout << _glob_x->_X[i] << "  ";
            //     // cout  << _delta * (_ub->_X[i] - _lb->_X[i]) << "   <   "<< fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i]) << endl;
            //     // if (_delta * (_ub->_X[i] - _lb->_X[i]) < fabs(_x_nearest_to_glob_x->_X[i] - _glob_x->_X[i])) {   // Infinity norm
            //     //     return false;
            //     // };
            // };
            // cout << endl;
            // cout << _delta << endl;
            return true;
        } else if (_stopping_criteria == "pe0.01") {
            e = 0.01;
        } else if (_stopping_criteria == "pe1") {
            e = 1.0;
        };
        double pe;
        if (_glob_f == 0) {
            pe = 100 * _f_min;
        } else {
            pe = 100 * (_f_min - _glob_f) / fabs(_glob_f);
        };
        if (pe < e) {
            return true;
        };
        return false;
    };

    virtual double value(Point* point) = 0;

    virtual ~Function(){
        delete _points;
    };
};
bool Function::_global_mem_allocated = false;


class Branin : public Function {
    Branin(const Branin& other){};
    Branin& operator=(const Branin& other){};
public:
    Branin(string stopping_criteria="pe1"): Function(){
        _name = "Branin";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-5., 0.);
        _ub = new Point(10., 15.);
        _glob_x = new Point(M_PI, 12.275);  // Point where global function minimum is (should be list)
        _glob_f = 0.397887;          // Predefined global function minimum
        _L = 109.94813585;
    };

    double value(Point* point) {
        double x1 = transform(point, 0); 
        double x2 = transform(point, 1); 
        // point->print();
        // cout << "    Transformed: " << x1 << ", " << x2 <<  "    lb " << _lb->_X[0] << "," << _lb->_X[1] << "  ub " << _ub->_X[0] << "," << _ub->_X[1] << endl;
        double part1 = pow((x2 - 5./(4*pow(M_PI, 2))*pow(x1,2) + 5./M_PI*x1 -6), 2);
        double part2 = 10.*(1. - 1./(8*M_PI))*cos(x1) + 10.;
        return part1 + part2;
    };
};


class GoldsteinPrice: public Function {   // http://www.sfu.ca/~ssurjano/goldpr.html
    GoldsteinPrice(const GoldsteinPrice& other){};
    GoldsteinPrice& operator=(const GoldsteinPrice& other){};
public:
    GoldsteinPrice(string stopping_criteria="pe1"): Function(){
        _name = "GoldsteinPrice";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-2., -2.);
        _ub = new Point(2., 2.);
        _glob_x = new Point(0., -1.);  // Point where global function minimum is (should be list)
        _glob_f = 3.;                 // Predefined global function minimum
        _L = 2225891.74508;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        double part1 = 1 + pow(x1 + x2 +1, 2) * (19 - 14*x1 + 3*pow(x1,2) - 14*x2 + 6*x1*x2 + 3*pow(x2,2));
        double part2 = 30 + pow(2*x1 - 3*x2, 2) * (18 - 32*x1 + 12*pow(x1,2) + 48*x2 - 36*x1*x2 + 27*pow(x2,2));
        return part1 * part2;
    };
};



class SHCamel: public Function {   // http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=44040 http://www.sfu.ca/~ssurjano/camel6.html
    SHCamel(const SHCamel& other){};
    SHCamel& operator=(const SHCamel& other){};
public:
    SHCamel(string stopping_criteria="pe1"): Function(){
        _name = "SHCamel";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-3., -2.);
        _ub = new Point(3., 2.);
        _glob_x = new Point(0.0898, -0.7126);  // Point where global function minimum is (should be list)
        _glob_f = -1.0316;                 // Predefined global function minimum
        _L = 689.202901909;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        return (4 - 2.1*pow(x1, 2) + pow(x1, 4) / 3.) * pow(x1, 2) + x1*x2 + (-4 + 4*pow(x2, 2))*pow(x2, 2);
    };
};

class Shubert: public Function {   // http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=44040 http://www.sfu.ca/~ssurjano/shubert.html
    Shubert(const Shubert& other){};
    Shubert& operator=(const Shubert& other){};
public:
    Shubert(string stopping_criteria="pe1"): Function(){
        _name = "Shubert";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-10., -10.);
        _ub = new Point(10., 10.);
        _glob_x = new Point(4.85805, 5.4828);  // Point where global function minimum is (should be list)
        _glob_f = -186.7309;                 // Predefined global function minimum
        _L = 59.4020075361;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);

        double sum1 = 0;
        double sum2 = 0;
        for (int i=1; i < 6; i++) {
            double new1 = i * cos((i+1)*x1 + i);
            double new2 = i * cos((i+1)*x2 + i);
            sum1 += new1;
            sum2 += new2;
        };
        return sum1 * sum2;
    };
};



class Alolyan: public Function {   // http://link.springer.com/article/10.1007%2Fs10898-012-0020-3
    Alolyan(const Alolyan& other){};
    Alolyan& operator=(const Alolyan& other){};
public:
    Alolyan(string stopping_criteria="pe1"): Function(){
        _name = "Alolyan";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-1., -1.);
        _ub = new Point(1., 1.);
        _glob_x = new Point(-1/3., 1.);  // Point where global function minimum is (should be list)
        _glob_f = -1.18519;                 // Predefined global function minimum
        _L = 5.65685424949;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        return x1*pow(x2, 2) + x2*pow(x1, 2) - pow(x1, 3) - pow(x2, 3);
    };
};

class Easom: public Function {   // http://link.springer.com/article/10.1007%2Fs10898-012-0020-3
    Easom(const Easom& other){};
    Easom& operator=(const Easom& other){};
public:
    Easom(string stopping_criteria="pe1"): Function(){
        _name = "Easom";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-30., -30.);
        _ub = new Point(30., 30.);
        // Leistinoji sritis nuo -100 iki 100

        _glob_x = new Point(M_PI, M_PI);  // Point where global function minimum is (should be list)
        _glob_f = -1.0;                   // Predefined global function minimum
        _L = 5.01891948878e-05;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        return -cos(x1) * cos(x2) * pow(M_E, -(pow(x1 - M_PI, 2) + pow(x2 - M_PI, 2)));
    };
};

class Rastrigin: public Function {   // http://link.springer.com/article/10.1007%2Fs10898-012-0020-3
    Rastrigin(const Rastrigin& other){};
    Rastrigin& operator=(const Rastrigin& other){};
public:
    Rastrigin(string stopping_criteria="pe1"): Function(){
        _name = "Rastrigin";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-5., -5.);
        _ub = new Point(6., 6.);
        _glob_x = new Point(0., 0.);  // Point where global function minimum is (should be list)
        _glob_f = 0.0;                   // Predefined global function minimum
        _L = 15.6204993518;
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        return 2*10. + pow(x1, 2) + pow(x2, 2) - 10*(cos(2*M_PI*x1) + cos(2*M_PI *x2));
    };
};


class Hartman3: public Function {   // http://www.sfu.ca/~ssurjano/hart3.html
    Hartman3(const Hartman3& other){};
    Hartman3& operator=(const Hartman3& other){};
public:
    Hartman3(string stopping_criteria="pe1"): Function(){
        _name = "Hartman3";
        _stopping_criteria = stopping_criteria;
        _D = 3;
        _lb = new Point(0., 3);
        _ub = new Point(1., 3);
        double X[3] = {0.114614, 0.555649, 0.852547};
        _glob_x = new Point(X, 3);  // Point where global function minimum is (should be list)
        _glob_f = -3.86278;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double alpha[4] = {1.0, 1.2, 3.0, 3.2};
        double A[4][3] = {{3.0, 10, 30}, {0.1, 10, 35}, {3.0, 10, 30}, {0.1, 10, 35}};
        double P[4][3] = {{0.3689, 0.1170, 0.2673},
                          {0.4699, 0.4387, 0.7470},
                          {0.1091, 0.8732, 0.5547},
                          {0.0381, 0.5743, 0.8828}};
        double sum1 = 0;
        for (int i=0; i < 4; i++) {
            double sum2 = 0;
            for (int j=0; j < 3; j++) {
                sum2 += A[i][j] * pow(X[j] - P[i][j], 2);
            };
            sum1 += alpha[i] * exp(-sum2);
        };
        return -sum1;
    };
};




class Shekel5: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
    Shekel5(const Shekel5& other){};
    Shekel5& operator=(const Shekel5& other){};
public:
    Shekel5(string stopping_criteria="pe1"): Function(){
        _name = "Shekel5";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(10., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.1532;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double m = 5;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;

        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};

class ReducedShekel5: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page2354.htm
    ReducedShekel5(const ReducedShekel5& other){};
    ReducedShekel5& operator=(const ReducedShekel5& other){};
public:
    ReducedShekel5(string stopping_criteria="pe1"): Function(){
        _name = "ReducedShekel5";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(6., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.1532;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double m = 5;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;

        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};


class Shekel7: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html
    Shekel7(const Shekel7& other){};
    Shekel7& operator=(const Shekel7& other){};
public:
    Shekel7(string stopping_criteria="pe1"): Function(){
        _name = "Shekel7";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(10., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.4029;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double m = 7;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;

        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};

class ReducedShekel7: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html
    ReducedShekel7(const ReducedShekel7& other){};
    ReducedShekel7& operator=(const ReducedShekel7& other){};
public:
    ReducedShekel7(string stopping_criteria="pe1"): Function(){
        _name = "ReducedShekel7";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(6., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.4029;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double m = 7;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;

        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};


class Shekel10: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html
    Shekel10(const Shekel10& other){};
    Shekel10& operator=(const Shekel10& other){};
public:
    Shekel10(string stopping_criteria="pe1"): Function(){
        _name = "Shekel10";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(10., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.5364;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        int m = 10;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;
        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};

class ReducedShekel10: public Function {   // http://www.sfu.ca/~ssurjano/shekel.html
    ReducedShekel10(const ReducedShekel10& other){};
    ReducedShekel10& operator=(const ReducedShekel10& other){};
public:
    ReducedShekel10(string stopping_criteria="pe1"): Function(){
        _name = "ReducedShekel10";
        _stopping_criteria = stopping_criteria;
        _D = 4;
        _lb = new Point(0., 4);
        _ub = new Point(6., 4);
        _glob_x = new Point(4., 4);  // Point where global function minimum is (should be list)
        _glob_f = -10.5364;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        int m = 10;
        double beta[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
        double C[4][10] = {{4.,1.,8.,6.,3.,2.,5.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,5.,1.,2.,3.6},
                           {4.,1.,8.,6.,3.,2.,3.,8.,6.,7.},
                           {4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6}};
        double sum1 = 0;
        for (int i=0; i < m; i++) {
            double sum2 = 0;
            for (int j=0; j < 4; j++) {
                sum2 += pow(X[j] - C[j][i], 2);
            };
            sum1 += 1./(sum2 + beta[i]);
        };
        return -sum1;
    };
};

class Hartman6: public Function {   // http://www.sfu.ca/~ssurjano/hart6.html
    Hartman6(const Hartman6& other){};
    Hartman6& operator=(const Hartman6& other){};
public:
    Hartman6(string stopping_criteria="pe1"): Function(){
        _name = "Hartman6";
        _stopping_criteria = stopping_criteria;
        _D = 6;
        _lb = new Point(0., _D);
        _ub = new Point(1., _D);
        double X[6] = {0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573};
        _glob_x = new Point(X, _D);  // Point where global function minimum is (should be list)
        _glob_f = -3.322;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double X[_D];
        for (int i=0; i < _D; i++) {
            X[i] = transform(point, i);
        };

        double alpha[4] = {1.0, 1.2, 3.0, 3.2};
        double A[4][10] = {{10, 3, 17, 3.5, 1.7, 8},
                           {0.05, 10, 17, 0.1, 8, 14},
                           {3.0, 3.5, 1.7, 10, 17, 8},
                           {17, 8, 0.05, 10, 0.1, 14}};
        double P[4][6] = {{0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
                          {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
                          {0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650},
                          {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}};
        double sum1 = 0;
        for (int i=0; i < 4; i++) {
            double sum2 = 0;
            for (int j=0; j < 6; j++) {
                sum2 += A[i][j] * pow(X[j] - P[i][j], 2);
            };
            sum1 += alpha[i] * exp(-sum2);
        };
        return -sum1;
    };
};

class JennrichSampson: public Function {   // http://infinity77.net/global_optimization/test_functions_nd_J.html
    JennrichSampson(const JennrichSampson& other){};
    JennrichSampson& operator=(const JennrichSampson& other){};
public:
    JennrichSampson(string stopping_criteria="pe1"): Function(){
        _name = "JennrichSampson";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(0., 0.);
        _ub = new Point(1., 1.);
        _glob_x = new Point(0.257825, 0.257825);  // Point where global function minimum is (should be list)
        _glob_f = 124.3621824;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        double sum1 = 0;
        for (int i=1; i < 11; i++ ) {
            sum1 += pow(2 + 2*i - (exp(i*x1) + exp(i*x2)), 2);
        };
        return sum1;
    };
};

class CenteredJennrichSampson: public Function {   // http://infinity77.net/global_optimization/test_functions_nd_J.html
    CenteredJennrichSampson(const CenteredJennrichSampson& other){};
    CenteredJennrichSampson& operator=(const CenteredJennrichSampson& other){};
public:
    CenteredJennrichSampson(string stopping_criteria="pe1"): Function(){
        _name = "CenteredJennrichSampson";
        _stopping_criteria = stopping_criteria;
        _D = 2;
        _lb = new Point(-0.5, -0.5);
        _ub = new Point(0.5, 0.5);
        _glob_x = new Point(0.257825, 0.257825);  // Point where global function minimum is (should be list)
        _glob_f = 124.3621824;                   // Predefined global function minimum
    };

    double value(Point* point) {
        double x1 = transform(point, 0);
        double x2 = transform(point, 1);
        double sum1 = 0;
        for (int i=1; i < 11; i++ ) {
            sum1 += pow(2 + 2*i - (exp(i*x1) + exp(i*x2)), 2);
        };
        return sum1;
    };
};

// def jennrich_sampson(X):
//     x1, x2 = X



// class RastriginShrinked: public Function {
//     RastriginShrinked(const RastriginShrinked& other){};
//     RastriginShrinked& operator=(const RastriginShrinked& other){};
// public:
//     RastriginShrinked(): Function(){
//         _name = "RastriginShrinked";
//         _D = 2;
//         _lb = new Point(-0.5, -0.5);
//         _ub = new Point(1.25, 1.25);
//         _glob_x = new Point(0., 0.);
//         _glob_f = 0;
//         _L = 15.6204993518;
//     };
//
//     double value(Point* point) {
//         double x1 = transform(point, 0);
//         double x2 = transform(point, 1);
//         return 2*10 + 4*pow(x1,2) + 4*pow(x2,2) - 10*(cos(2*M_PI*x1) + cos(2*M_PI*x2));
//     };
//
//     virtual ~RastriginShrinked(){
//         delete _lb;
//         delete _ub;
//         delete _glob_x;
//     };
// };


class GKLSFunction: public Function {
    GKLSFunction(const GKLSFunction& other){};
    GKLSFunction& operator=(const GKLSFunction& other){};
public:
    GKLSFunction(int cls, int function_id): Function() {
        int _GKLS_class_D[] = {2, 2, 3, 3, 4, 4, 5, 5};
        double _GKLS_class_global_dists[] = {0.9, 0.9, 0.66, 0.9, 0.66, 0.9, 0.66, 0.66};
        double _GKLS_class_global_radiuses[] = {0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2};
        double _GKLS_class_detlas[] = {1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-7, 1e-7};

        stringstream function_name; 
        function_name << cls << "_" << function_id;  
        _name =  function_name.str();
        _cls = cls;
        _fid = function_id;
        _stopping_criteria = "x_dist";

        cls -= 1;
        _D = _GKLS_class_D[cls];
        _global_dist = _GKLS_class_global_dists[cls];
        _global_radius = _GKLS_class_global_radiuses[cls]; 
        _delta = pow(_GKLS_class_detlas[cls], 1./_D);   // _delta = _GKLS_class_detlas[cls];

        _lb = new Point(-1., _D);
        _ub = new Point(1., _D);
        _x_nearest_to_glob_x = _lb;
        _distance_to_glob_x = numeric_limits<double>::max();

        _glob_f = -1.;          // Predefined global function minimum
        // _L = ;

        // assert(GKLS_set_default()== GKLS_OK);     // Standartiniai nustatymai

        GKLS_dim = _D;
        GKLS_global_dist = _global_dist;
        GKLS_global_radius = _global_radius;
        GKLS_num_minima = 10;
        GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
        if (Function::_global_mem_allocated == true) {
            GKLS_free();
        } else {
            assert(GKLS_domain_alloc() == GKLS_OK);
            Function::_global_mem_allocated = true;
        };
        // for (unsigned int i = 0; i < GKLS_dim; i++) {
        //     GKLS_domain_left[i] = -1;
        //     GKLS_domain_right[i] = 1;
        // };
        assert(GKLS_parameters_check() == GKLS_OK);
        assert(GKLS_arg_generate(function_id) == GKLS_OK);
        int n = GKLS_glob.num_global_minima;
        assert(n == 1);
        int glob_idx = 1;
        assert(GKLS_minima.f[glob_idx] == GKLS_global_value);
        _glob_x = new Point(_D);  // Point where global function minimum is (should be list)
        for (int i=0; i < _D; i++) {
            _glob_x->_X[i] = (GKLS_minima.local_min[glob_idx][i] - _lb->_X[i]) / (_ub->_X[i]-_lb->_X[i]);
        };
        _glob_x->add_value(_glob_f);
    };

    int _fid;
    int _cls;
    double _global_dist;
    double _global_radius;


    double transform(Point* point, int i) {     // Transforms single point coordinate from [0,1] to [l,u]
        return point->_X[i] * (_ub->_X[i]-_lb->_X[i]) + _lb->_X[i];  
    };

    double value(Point* point) {
        GKLS_free();
        assert(GKLS_arg_generate(_fid) == GKLS_OK);   // Needed for multicriteria problems. Slows down function evaluations.
        double transformed_point[_D];
        for (int i=0; i<_D; i++){
            transformed_point[i] = transform(point,i);
        };
        double value = GKLS_D_func(transformed_point);
        return value;
    };

    virtual ~GKLSFunction(){
        delete _lb;
        delete _ub;
        delete _glob_x;
        if (Function::_global_mem_allocated) {
            GKLS_free();
            GKLS_domain_free();
            Function::_global_mem_allocated = false;
        };
    };
};

Function* get_function(char* func_name, string stopping_criteria="pe1"){
    if (!strcmp(func_name, "branin")) { return new Branin(stopping_criteria); };
    if (!strcmp(func_name, "goldsteinprice")) { return new GoldsteinPrice(stopping_criteria); };
    if (!strcmp(func_name, "camel")) { return new SHCamel(stopping_criteria); };
    if (!strcmp(func_name, "shubert")) { return new Shubert(stopping_criteria); };
    if (!strcmp(func_name, "alolyan")) { return new Alolyan(stopping_criteria); };
    if (!strcmp(func_name, "easom")) { return new Easom(stopping_criteria); };
    if (!strcmp(func_name, "rastrigin")) { return new Rastrigin(stopping_criteria); };
    if (!strcmp(func_name, "hartman3")) { return new Hartman3(stopping_criteria); };
    if (!strcmp(func_name, "shekel5")) { return new Shekel5(stopping_criteria); };
    if (!strcmp(func_name, "shekel7")) { return new Shekel7(stopping_criteria); };
    if (!strcmp(func_name, "shekel10")) { return new Shekel10(stopping_criteria); };
    if (!strcmp(func_name, "hartman6")) { return new Hartman6(stopping_criteria); };
    if (!strcmp(func_name, "reducedshekel5")) { return new ReducedShekel5(stopping_criteria); };
    if (!strcmp(func_name, "reducedshekel7")) { return new ReducedShekel7(stopping_criteria); };
    if (!strcmp(func_name, "reducedshekel10")) { return new ReducedShekel10(stopping_criteria); };
    if (!strcmp(func_name, "jennrichsampson")) { return new JennrichSampson(stopping_criteria); };
    if (!strcmp(func_name, "centeredjennrichsampson")) { return new CenteredJennrichSampson(stopping_criteria); };
};

// vector<Function*> functions;
// functions.push_back(Branin);

#endif
