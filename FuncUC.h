/* Copyright Albertas Gimbutas 2017, all rights reserved */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H 
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>

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
    double* _X;      // Coordinates in normalised [0,1]^n space  
    double _value;   // Objective function value
    vector<Simplex*> _simplices;  // Simplices, which have this point as vertex
         
    void add_value(double value) {
        _value = value;
    };

    void print(){
        // cout.precision(17);
        cout << "       ";
        for (int i=0; i < size(); i++){
            cout << _X[i] << "  \t";
        };
        cout << "-> " << _value << "  ";
        // cout << endl;
    };

    int size(){
        return _D;
    };

    static bool ascending_value(Point* p1, Point* p2) {
        return p1->_value < p2->_value;
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


class Function {        // Abstract function class
    Function(const Function& other){};
    Function& operator=(const Function& other){};
public:
    Function(){
        _f_min = numeric_limits<int>::max();
        _points = new PointTree();
        _evaluations = 0;
    };

    int _evaluations;
    int _D;                     // Dimension
    double _f_min;              // Best known function value
    PointTree* _points;         // Binary balancing tree to store points where objective function was evaluated
    vector<Point*> _new_points; // Points for which stopping condition was not checked yet

    void add_value(Point* p) {
        double val = value(p);
        _evaluations += 1;
        p->add_value(val);
        if (_f_min > val) {
            _f_min = val;
        };
    };

    Point* get(double *c, int D){  
        // Returns a point with objective function value (the point may be from cache)
        Point* p = new Point(c, D);
        Point* cached_point = _points->add(p);
        if (cached_point) {   // Value at this point is already known
            delete p;
            return cached_point;
        } else {              // Value at this point is unknown, evaluate it
            add_value(p);
            _new_points.push_back(p);
            return p;
        };
    };

    Point* get(Point* p) {   
        // Returns a point with objective function value (the point may be from cache)
        Point* cached_point = _points->add(p);
        if (cached_point) {   // Value at this point is already known
            return cached_point;
        } else {              // Value at this point is unknown, evaluate it
            add_value(p);
            _new_points.push_back(p);
            return p;
        };
    };

    virtual bool is_accurate_enough() = 0;

    virtual double value(Point* point) = 0;

    virtual ~Function(){
        delete _points;
    };
};


class FuncUC: public Function {      // Function which is define over a unit-cube
    FuncUC(const FuncUC& other){};
    FuncUC& operator=(const FuncUC& other){};
public:
    FuncUC(int D, double (*get_value_uc)(vector<double>), bool (*should_stop_uc)(vector<double>)){
        _D = D;
        get_value = get_value_uc;
        should_stop = should_stop_uc;
    };

    double (*get_value) (vector<double>);   // Objective function evaluation method provided as an argument
    bool (*should_stop) (vector<double>);   // Stopping criterion method provided as an argument

    double value(Point* point) {
        // Converts a point object to a vector and evaluates objective value at it
        vector<double> point_vector;
        for (int i=0; i < point->_D; i++) {
            point_vector.push_back(point->_X[i]);
        };
        return get_value(point_vector);
    };

    bool is_accurate_enough() {
        // Checks stopping criterion for all new points
        int nr_of_new_points = _new_points.size();
        for (int j=0; j < nr_of_new_points; j++) {
            // Pop one of the new points
            Point* p = _new_points.back();      
            _new_points.pop_back();

            // Convert that point object to a vector
            vector<double> point_vector;        
            for (int i=0; i < p->_D; i++) {
                point_vector.push_back(p->_X[i]);
            };

            // Check stopping condition
            if (should_stop(point_vector) == true) {  
                return true;
            };
        };
        return false;
    };
};

#endif
