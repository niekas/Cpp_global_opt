/* Copyright Albertas Gimbutas 2017, all rights reserved */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>
#include <math.h>
#include <string.h>
#include <fstream>
#include <sstream>

// #include "utils.h"

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
    Point(int *c, int D){
        _D = D;
        _X = (double*) malloc((D)*sizeof(double));
        for (int i=0; i<D; i++){
            _X[i] = double(c[i]);
        };
    };
    Point(double *c, int D){
        _D = D;
        _X = (double*) malloc((D)*sizeof(double));
        for (int i=0; i<D; i++){
            _X[i] = c[i];
        };
    };

    int _D;          // Dimension of variable space
    int _C;
    double* _X;      // Coordinates in normalised [0,1]^n space
    vector<double> _values;   // Objective function value
    vector<Simplex*> _simplices;  // Simplices, which have this point as vertex

    void add_values(vector<double> values) {
        for (int i=0; i < values.size(); i++) {
            _values.push_back(values[i]);
        };
        _C = values.size();
    };

    void print(){
        cout.precision(17);
        cout << "       ";
        for (int i=0; i < size(); i++){
            cout << _X[i] << "  \t";
        };
        cout << "-> ";
        for (int i=0; i < _values.size(); i++) {
             cout << _values[i] << "  ";
        };
        // cout << endl;
    };

    int size(){
        return _D;
    };

    static bool ascending_value(Point* p1, Point* p2) {
        // Ascending by first value, if equal then by second and so on
        for (int c=0; c < p1->_values.size(); c++) {
            if (p1->_values[c] < p2->_values[c]) {
                return p1->_values[c] < p2->_values[c];
            };
            if (p1->_values[c] > p2->_values[c]) {
                return p1->_values[c] < p2->_values[c];
            };
        };
    };

    virtual ~Point(){
        free(_X);
    };
};


// Binary balancing tree (or simply linked list) for storing points
// It returns cached point if a point with the same coordinates is added
class PointTree;

class PointTreeNode {
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

    virtual ~PointTreeNode();
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

    Point* process_next_dimension(PointTreeNode* node, Point* point) {
        // Creates next dimension tree if needed and adds point to it its last dimension
        if (point->_D == _dim) {        // Don't need next dimension
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


//////////////////////////////////////////////
//// Multi-objective function definitions ////
//////////////////////////////////////////////

class Function {  // Abstract class to store information specific to a function (optimization problem)
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(int D) {};

    string _name;          // Function name
    int _D;                // Number of dimensions in variable space
    int _C;                // Number of objectives
    vector<double> _lb;    // Lower bound values
    vector<double> _ub;    // Upper bound values
    vector<double> _nadir; // Nadir point

    virtual vector<double> get_values(vector<double> X) = 0;
    virtual ~Function(){};
};


class FunctionUC {      // Abstract multi-objective function defined over a unite-cube
    FunctionUC(const FunctionUC& other){};
    FunctionUC& operator=(const FunctionUC& other){};
public:
    FunctionUC(Function* func){
        _func = func;    // Unable to use:  get_values_not_uc = func->get_values;
        _lb = func->_lb;
        _ub = func->_ub;
        _nadir = func->_nadir;
        _D = func->_D;
        _C = func->_C;
        _points = new PointTree();
        _evaluations = 0;
    };

    string _name;
    int _D;                 // Number of dimensions in variable space
    int _C;                 // Number of objectives
    vector<double> _nadir;
    vector<double> _lb;
    vector<double> _ub;
    PointTree* _points;     // Binary balancing tree to store points where objective function was evaluated
    vector<Point*> _pareto_front;
    bool _pareto_front_was_updated;
    int _evaluations;
    Function* _func;

    //// Evaluation methods
    vector<double> transform_from_uc(double* X_uc) {
        vector<double> X;
        for (int i=0; i < _D; i++){
            X.push_back(X_uc[i] * (_ub[i] - _lb[i]) + _lb[i]);
        };
        return X;
    };

    vector<double> get_values(Point* p_uc) {  // Returns values for point in UC
        return _func->get_values(transform_from_uc(p_uc->_X));
    };


    void add_values(Point* p_uc) {  // Adds obj values to point in UC
        vector<double> values = get_values(p_uc);
        _evaluations += 1;
        p_uc->add_values(values);
        if (update_pareto_front(p_uc) == true){
            _pareto_front_was_updated = true;
        };
    };

    Point* get_point_with_values(double *c_uc, int D) {
        Point* p = new Point(c_uc, D);
        Point* cached_point = _points->add(p);
        if (cached_point) {   // Value at this point is already known
            delete p;
            return cached_point;
        } else {              // Value at this point is unknown, evaluate it
            add_values(p);
            return p;
        };
    };

    Point* get_point_with_values(Point* p) {
        Point* cached_point = _points->add(p);
        if (cached_point) {   // Value at this point is already known
            return cached_point;
        } else {              // Value at this point is unknown, evaluate it
            add_values(p);
            return p;
        };
    };

    //// Pareto front methods
    int domination(Point* q, Point* p) {
        // 0 - q domiates p
        // 1 - none of them dominates
        // 2 - p dominates q
        int q_better = 0;
        int p_better = 0;
        for (int i=0; i < q->_values.size(); i++) {
            if (q->_values[i] < p->_values[i]) { q_better++; };
            if (q->_values[i] > p->_values[i]) { p_better++; };
        };
        if (q_better == q->_values.size()) { return 0; }
        if (p_better == q->_values.size()) { return 2; }
        return 1;
    };

    bool update_pareto_front(Point* p) { // Returns if point was added to pareto front
        // if any(q > p for q in S):
        //     return
        // for q in [q for q in S if p > q]:
        //     S.remove(q)
        //     S.add(p)

        vector<int> dominated;
        int relation;
        for (int i=0; i < _pareto_front.size(); i++) {
            relation = domination(_pareto_front[i], p);  // domination returns 0 > , 1 = , 2 <
            if (relation == 0) { return false; };
            if (relation == 2) { dominated.push_back(i); };
        };
        for (int i = dominated.size() - 1; i >= 0; i--) {
            _pareto_front.erase(_pareto_front.begin() + dominated[i], _pareto_front.begin() + dominated[i] +1);
        };
        _pareto_front.push_back(p);
        return true;
    };

    void show_pareto_front() {
        log_pareto_front();
        FILE* testp = popen("python log/show_pareto_front.py log/front.txt", "r");
        pclose(testp);
    };

    //// Pareto front metrics
    double values_l2norm(Point* p1, Point* p2) {
        double squared_sum = 0;
        for (int i=0; i < p1->_C; i++){
            squared_sum += pow(p1->_values[i] - p2->_values[i], 2);
        };
        return sqrt(squared_sum);
    };

    double uniformity() {
        vector<double> min_dists;
        double dist;
        for (int i=0; i < _pareto_front.size(); i++) {
            double min_dist = numeric_limits<double>::max();
            for (int j=0; j < _pareto_front.size(); j++) {
                if (i != j) {
                    dist = values_l2norm(_pareto_front[i], _pareto_front[j]);
                    if (dist < min_dist) {
                        min_dist = dist;
                    };
                };
            };
            min_dists.push_back(min_dist);
        };

        double avg_dist = 0;
        for (int i=0; i < min_dists.size(); i++) {
            avg_dist += min_dists[i];
        };
        avg_dist = avg_dist / min_dists.size();

        double uniformity = 0;
        for (int i=0; i < min_dists.size(); i++) {
           uniformity += pow(min_dists[i] - avg_dist, 2);
        };

        uniformity = sqrt(uniformity);
        return uniformity;
    };

    double hyper_volume() {
        sort(_pareto_front.begin(), _pareto_front.end(), Point::ascending_value);
        log_pareto_front();

        char buffer[50];
        FILE* fp = popen("python log/show_pareto_front.py log/front.txt -hv", "r");
        fgets(buffer, 10, fp);
        pclose(fp);
        return atof(buffer);
    };

    //// Visualization methods
    void log_pareto_front() {
        ofstream log_file;
        log_file.open("log/front.txt");
        log_file.close();
        log_file.open("log/front.txt", ios::app);
        log_file.precision(17);

        for (int i=0; i < _pareto_front.size(); i++) {
            for (int k=0; k < _pareto_front[i]->_D; k++) {
                log_file << _pareto_front[i]->_X[k] << " ";
            };
            log_file << "(";
            for (int c=0; c < _pareto_front[i]->_C; c++) {
                log_file << _pareto_front[i]->_values[c];
                if (c != _pareto_front[i]->_C -1) {
                    log_file << " ";
                };
            };
            log_file << ")" << endl;
        };

        log_file << "nadir: ";
        for (int c=0; c < _nadir.size(); c++) {
            log_file << _nadir[c] << " ";
        };
        log_file.close();
    };

    virtual ~FunctionUC(){
        delete _points;
    };
};


class EP1: public Function {             // Nonuniform Covering Method as Applied to Multicriteria Optimization Problems with Guaranteed Accuracy
    EP1(const EP1& other){};             // Yu. G. Evtushenko and M. A. Posypkin, first example
    EP1& operator=(const EP1& other){};  // Front with a break
public:
    EP1(int D) {
        _name = "EP1";
        _D = 2;
        _C = 2;
        _nadir.push_back(2.); _nadir.push_back(3.);
        _lb.push_back(0.); _lb.push_back(0.);
        _ub.push_back(2.); _ub.push_back(2.);
    };

    double value1(vector<double> X) {
        return X[0];
    };

    double value2(vector<double> X) {
        double p1 = fabs(X[0] - 1);
        double p2 = 1.5 - X[0];
        if (p1 < p2) {
            return p1 + X[1] + 1;
        };
        return p2 + X[1] + 1;
    };

    vector<double> get_values(vector<double> X) {
        vector<double> values;
        values.push_back(value1(X));
        values.push_back(value2(X));
        return values;
    };
};


class EP2: public Function {           // Nonuniform Covering Method as Applied to Multicriteria Optimization Problems with Guaranteed Accuracy
    EP2(const EP2& other){};             // Yu. G. Evtushenko and M. A. Posypkin, second example
    EP2& operator=(const EP2& other){};
public:
    EP2(int D) {
        _name = "EP2";
        _D = 2;
        _C = 2;
        _nadir.push_back(1.); _nadir.push_back(1.);
        _lb.push_back(0.); _lb.push_back(0.);
        _ub.push_back(1.); _ub.push_back(1.);
    };

    double value1(vector<double> X) {
        return (X[0] - 1) * X[1] * X[1] + 1;
    };

    double value2(vector<double> X) {
        return X[1];
    };

    vector<double> get_values(vector<double> X) {
        vector<double> values;
        values.push_back(value1(X));
        values.push_back(value2(X));
        return values;
    };
};


class GeneticFunction: public Function {
    GeneticFunction(const GeneticFunction& other){};
    GeneticFunction& operator=(const GeneticFunction& other){};
public:
    GeneticFunction(int D) {};
    vector<double> get_values(vector<double> X) {
        stringstream outs;
        outs << "python3 problems/genetic.py";
        outs << " " << _name;
        for (int i=0; i < X.size(); i++) {
            outs << " " << X[i];
        };
        FILE* p = popen(outs.str().c_str(), "r");
        char input[1500];
        if (p != NULL) {
            while(fgets(input, sizeof(input), p) != NULL) {};
        };
        string ins = string(input);
        stringstream ss(ins);

        double num;
        vector<double> values;
        while (!ss.eof()) {
            ss >> num;
            values.push_back(num);
        };
        return values;
    };
};


class ZDT1: public GeneticFunction {
    ZDT1(const ZDT1& other){};
    ZDT1& operator=(const ZDT1& other){};
public:
    ZDT1(int D) {
        _name = "zdt1";
        _C = 2;
        _D = 6;   if (D > 0) { _D = D };  // 30
        for (int c=0; c < _C; c++) {
            _nadir.push_back(11.);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class ZDT2: public GeneticFunction {
    ZDT2(const ZDT2& other){};
    ZDT2& operator=(const ZDT2& other){};
public:
    ZDT2(int D) {
        _name = "zdt2";
        _C = 2;
        _D = 6;   if (D > 0) { _D = D };  // 30
        for (int c=0; c < _C; c++) {
            _nadir.push_back(11.);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class ZDT3: public GeneticFunction {
    ZDT3(const ZDT3& other){};
    ZDT3& operator=(const ZDT3& other){};
public:
    ZDT3(int D) {
        _name = "zdt3";
        _C = 2;
        _D = 6;   if (D > 0) { _D = D };  // 30
        for (int c=0; c < _C; c++) {
            _nadir.push_back(11.);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class ZDT4: public GeneticFunction {
    ZDT4(const ZDT4& other){};
    ZDT4& operator=(const ZDT4& other){};
public:
    ZDT4(int D) {
        _name = "zdt4";
        _C = 2;
        _D = 6;   if (D > 0) { _D = D };  // 10
        for (int c=0; c < _C; c++) {
            _nadir.push_back(11.);
        };
        _lb.push_back(0.);
        _ub.push_back(1.);
        for (int i=0; i < _D -1; i++) {
            _lb.push_back(-5.);
            _ub.push_back(5.);
        };
    };
};

class ZDT6: public GeneticFunction {
    ZDT6(const ZDT6& other){};
    ZDT6& operator=(const ZDT6& other){};
public:
    ZDT6(int D) {
        _name = "zdt6";
        _C = 2;
        _D = 6;   if (D > 0) { _D = D };  // 10
        for (int c=0; c < _C; c++) {
            _nadir.push_back(11.);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};


class DTLZ1: public GeneticFunction {
    DTLZ1(const DTLZ1& other){};
    DTLZ1& operator=(const DTLZ1& other){};
public:
    DTLZ1(int D) {
        _name = "dtlz1";
        _C = 3;
        _D = 6;   if (D > 0) { _D = D };  // 30
                                       // Reference points taken from:
        for (int c=0; c < _C; c++) {   // Theory of the Hypervolume Indicator: Optimal μ-Distributions and the Choice of the Reference Point
            _nadir.push_back(1.);      // https://hal.inria.fr/inria-00430540/document
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class DTLZ2: public GeneticFunction {
    DTLZ2(const DTLZ2& other){};
    DTLZ2& operator=(const DTLZ2& other){};
public:
    DTLZ2(int D) {
        _name = "dtlz2";
        _C = 3;
        _D = 6;   if (D > 0) { _D = D };  // 30

        for (int c=0; c < _C; c++) {
            _nadir.push_back(1.5);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class DTLZ3: public GeneticFunction {
    DTLZ3(const DTLZ3& other){};
    DTLZ3& operator=(const DTLZ3& other){};
public:
    DTLZ3(int D) {
        _name = "dtlz3";
        _C = 3;
        _D = 6;   if (D > 0) { _D = D };  // 30

        for (int c=0; c < _C; c++) {
            _nadir.push_back(1.5);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };

    };
};

class DTLZ4: public GeneticFunction {
    DTLZ4(const DTLZ4& other){};
    DTLZ4& operator=(const DTLZ4& other){};
public:
    DTLZ4(int D) {
        _name = "dtlz4";
        _C = 3;
        _D = 6;   if (D > 0) { _D = D };  // 30
        for (int c=0; c < _C; c++) {
            _nadir.push_back(1.5);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class DTLZ5: public GeneticFunction {
    DTLZ5(const DTLZ5& other){};
    DTLZ5& operator=(const DTLZ5& other){};
public:
    DTLZ5(int D) {
        _name = "dtlz5";
        _C = 3;
        _D = 6;    if (D > 0) { _D = D }; // 30
        for (int c=0; c < _C; c++) {
            _nadir.push_back(1.5);   // Note: here 1.5 value is guessed
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class DTLZ6: public GeneticFunction {
    DTLZ6(const DTLZ6& other){};
    DTLZ6& operator=(const DTLZ6& other){};
public:
    DTLZ6(int D) {
        _name = "dtlz6";
        _C = 3;
        _D = 6;  if (D > 0) { _D = D };
        for (int c=0; c < _C; c++) {
            _nadir.push_back(1.5);   // Note: here 1.5 value is guessed
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};

class DTLZ7: public GeneticFunction {
    DTLZ7(const DTLZ7& other){};
    DTLZ7& operator=(const DTLZ7& other){};
public:
    DTLZ7(int D) {
        _name = "dtlz7";
        _C = 3;
        _D = 6;   if (D > 0) { _D = D };
        for (int c=0; c < _C; c++) {
            _nadir.push_back(15.);
        };
        for (int i=0; i < _D; i++) {
            _lb.push_back(0.);
            _ub.push_back(1.);

        };
    };
};


Function* get_function(char* func_name, int D=-1) {
    if (!strcmp(func_name, "ep1")) { return new EP1(D); };
    if (!strcmp(func_name, "ep2")) { return new EP2(D); };

    if (!strcmp(func_name, "zdt1")) { return new ZDT1(D); };
    if (!strcmp(func_name, "zdt2")) { return new ZDT2(D); };
    if (!strcmp(func_name, "zdt3")) { return new ZDT3(D); };
    if (!strcmp(func_name, "zdt4")) { return new ZDT4(D); };
    if (!strcmp(func_name, "zdt6")) { return new ZDT6(D); };

    if (!strcmp(func_name, "dtlz1")) { return new DTLZ1(D); };
    if (!strcmp(func_name, "dtlz2")) { return new DTLZ2(D); };
    if (!strcmp(func_name, "dtlz3")) { return new DTLZ3(D); };
    if (!strcmp(func_name, "dtlz4")) { return new DTLZ4(D); };
    if (!strcmp(func_name, "dtlz5")) { return new DTLZ5(D); };
    if (!strcmp(func_name, "dtlz6")) { return new DTLZ6(D); };
    if (!strcmp(func_name, "dtlz7")) { return new DTLZ7(D); };
};

#endif
