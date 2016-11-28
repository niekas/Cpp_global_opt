#ifndef SIMPLEX_H
#define SIMPLEX_H 
/* Simplex is a simplex used for solving main problem, e.g. minimizing GKLS function */
#include <list>
#include "Eigen/Dense"
#include "Elbme.h"
#include "PropConte.h"
#include "EdgeLbs.h"
#include "NelderMead.h"

using namespace std;


enum LowerBoundStrategy { MinVert, LongestEdgeLB, LowestEdgeLB, All };
const char* LBS[] = { "Min vert", "Longest edge LB", "Lowest edge LB", "All", 0 };
enum LStrategy { Self, Neighbours };
const char* LS[] = { "Self", "Neighbours", 0 };
enum DivisionStrategy { LongestHalf };
const char* DS[] = { "Longest Half", 0 };
enum SimplexGradientStrategy { FFMinVert, FFMaxVert, FFAllVertMean };
const char* SGS[] = { "FF min vert", "FF max vert", "FF all vert mean", 0 };   // "All vert max" should match "Max vert"

class SimplexTree;
class SimplexTreeNode;


class Simplex {  // Designed for outer problems
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex(LowerBoundStrategy lower_bound_strategy,
            LStrategy L_strategy,
            SimplexGradientStrategy simplex_gradient_strategy
        ) {
        _lower_bound_strategy = lower_bound_strategy;
        _L_strategy = L_strategy;
        _simplex_gradient_strategy = simplex_gradient_strategy;
        _is_in_partition = true;
        _diameter = 0;
        _tolerance = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _min_vert = 0;
        // _max_vert = 0;
        // _max_vert_value = -numeric_limits<double>::max();
        // _min_vert_value = numeric_limits<double>::max();
        _should_be_divided = false;
        _should_estimates_be_updated = true;
        // _min_lbs = 0;
        // _min_lb_value = 0;
        _D = 0;
    };

    LowerBoundStrategy _lower_bound_strategy;
    LStrategy _L_strategy;
    SimplexGradientStrategy _simplex_gradient_strategy;

    int _D;                       // Variable space dimension
    int _C;                       // Criteria (objective function) space dimension
    double _tolerance;            // Distance between _min_lb and _min_vert  or _min_lb_value and _min_vert_value

    vector<Point*> _verts;        // Simplex vertexes (points with coordinates and values)
    bool _is_in_partition;
    bool _should_be_divided;  // Should be divided in next iteration
    bool _should_estimates_be_updated;   // Should Lipschitz constant estimate and its lower bound be updated
    list<Simplex*> _neighbours;

    Point* _le_v1;      // Longest edge vertex1
    Point* _le_v2; 
    double _diameter;   // Longest edge length

    static vector<double> glob_Ls;
    static bool glob_L_was_updated;
    static double max_diameter;
    static double glob_L_coef;
    static double local_L_coef;

    // vector<double> _Ls;          // Cumulative estimates of Lipschitz constants for each criteria
    vector<double> _local_Ls;          // Cumulative estimates of Lipschitz constants for each criteria
    // vector<double> _grad_norms;  // Lipschitz constant estimate calculated by Simplex Gradient Euclidean norm.
    vector<double> _simpl_Ls;  // Lipschitz constant estimate calculated by Simplex Gradient Euclidean norm.

    Point* _min_vert;   // Pointer to vertex with lowest function value 
    // double _min_vert_value;  // _min_vert function value 
    // Point* _max_vert;
    // double _max_vert_value;
    // double _metric__vert_min_value;     // _f_min - glob_f / _diameter

    vector<Point*> _min_lbs;
    // double _min_lb_value;
    // double _metric__min_lb;       // _f_min - glob_f / _diameter
    
    void init_parameters(vector<Function*> funcs) {   // Called when all verts have been added
        _D = _verts.size() - 1;
        _C = funcs.size();
        for (int i=0; i < _C; i++) {
            _local_Ls.push_back(0.);
            _simpl_Ls.push_back(0.);
        };
        
        // Note: claculating metrics needed by algorithm here would reduce calculations

        // Sorts vertexes using first criteria values
        sort(_verts.begin(), _verts.end(), Point::compare_by_value);  // Is there calculation of intersection anywhere in the algorithm, in this case sorting by adress is needed?

        // Find longest edge length and verts
        double edge_length;  // Temporary variable
        for (int a=0; a < _verts.size(); a++) {
            // Finds _diameter
            for (int b=0; b < _verts.size(); b++){
                if (b > a) {
                    edge_length = l2norm(_verts[a], _verts[b]); 
                    if (edge_length > _diameter) {
                        _diameter = edge_length;
                        _le_v1 = _verts[a];
                        _le_v2 = _verts[b];
                    };
                };
            };
        }; 

        // Find adaptive Lipschitz constant
        // double E;
        // if (1e-4 * fabs(funcs[0]->_glob_f) > 1e-8) {
        //     E = 1e-4 * fabs(funcs[0]->_glob_f);
        // } else {
        //     E = 1e-8;
        // };


        for (int i=0; i < funcs.size(); i++) {
            _simpl_Ls[i] = find_simplex_gradient_norm(i, _simplex_gradient_strategy);      // Check in the article if global Lipschitz constant is defined

           //// Update global Ls
           if (Simplex::glob_Ls.size() < funcs.size()) {
               Simplex::glob_Ls.push_back(_simpl_Ls[i]);
           } else {
               if (Simplex::glob_Ls[i] < Simplex::glob_L_coef * _simpl_Ls[i]) {
                   Simplex::glob_Ls[i] = Simplex::glob_L_coef * _simpl_Ls[i];
                   Simplex::glob_L_was_updated = true;
               };
           };
        };

        _min_vert = _verts[0];
        for (int i=0; i < _verts.size(); i++) {
            if (_verts[i]->_values[0] < _min_vert->_values[0]) {
                _min_vert = _verts[i];
            };
        };
        // ToDo: _simpl_Ls - is not representative title, should be renamed
        // to mark that these are Lipschitz constant estimates for this simplex.
    };

    // Need a scenario where a single simplex is created and I can test with it  
    vector<Point*> find_accurate_lb_min_estimates(vector<Point*> verts, vector<double> Ls) {
        vector<Point*> estimates_of_accurate_lb_min;
        for (int i=0; i < Ls.size(); i++) {
            // Elbme* alg = new Elbme(verts, Ls, i);
            // PropConte* alg = new PropConte(verts, Ls, i, diameter);
            // NelderMead* alg = new NelderMead(verts, Ls, i, diameter);
            EdgeLbs* alg = new EdgeLbs(verts, Ls, i);
            Point* estimate_of_accurate_lb_min = alg->minimize();
            estimates_of_accurate_lb_min.push_back(estimate_of_accurate_lb_min->copy());
            delete estimate_of_accurate_lb_min;
            delete alg;
        };
        return estimates_of_accurate_lb_min;
    };

    vector<Point*> find_two_best_vert_lb_mins(Simplex* simpl, vector<double> Ls) {
        vector<Point*> lb_mins;
        for (int i=0; i < Ls.size(); i++) {    // Iterate through criterias
            // Find best vert for this criteria 
            Point* min_vert = simpl->_verts[0];
            Point* min_vert2 = simpl->_verts[0];
            Point* max_vert = simpl->_verts[0];
            for (int j=0; j < simpl->_verts.size(); j++) {
                if (simpl->_verts[j]->_values[i] < min_vert->_values[i]) {
                    min_vert = simpl->_verts[j];
                };
                if (simpl->_verts[j]->_values[i] > max_vert->_values[i]) {
                    max_vert = simpl->_verts[j];
                };
            };
            if (min_vert2 == min_vert) {
                min_vert2 = simpl->_verts[1];
            };
            for (int j=0; j < simpl->_verts.size(); j++) {
                if ((simpl->_verts[j]->_values[i] < min_vert2->_values[i]) and (min_vert != simpl->_verts[j])) {
                    min_vert2 = simpl->_verts[j];
                };
            };

            // double dist = l2norm(min_vert, min_vert2);

            double L = Ls[i];

            // Kokia formule apskaiciuoti apatinei ribai naudojant dvi virsunes?
            // double lb_value = (min_vert->_values[i] + min_vert2->_values[i] - L * simpl->_diameter) / 2.;
            double lb_value = (min_vert->_values[i] + max_vert->_values[i] - L * simpl->_diameter) / 2.;

            // cout << lb_value << " = " << min_vert->_values[i] << " - " << L << " * " << simpl->_diameter;

            Point* lb_min = new Point(simpl->_D);
            lb_min->add_value(lb_value);

            lb_mins.push_back(lb_min->copy());
            delete lb_min;
        };
        return lb_mins;
    };

    vector<Point*> find_one_vert_lb_mins(Simplex* simpl, vector<double> glob_Ls) {
        vector<Point*> lb_mins;

        for (int i=0; i < _C; i++) {    // Iterate through criterias
            // Find best vert for this criteria 
            Point* min_vert = simpl->_verts[0];
            for (int j=0; j < simpl->_verts.size(); j++) {
                // simpl->_verts[j]->print();
                if (simpl->_verts[j]->_values[i] < min_vert->_values[i]) {
                    min_vert = simpl->_verts[j];
                };
            };

            double L = simpl->_local_Ls[i];
            double lb_value = min_vert->_values[i] - L * simpl->_diameter;
            // cout << lb_value << " = " << min_vert->_values[i] << " - " << L << " * " << simpl->_diameter;

            Point* lb_min = new Point(simpl->_D);
            lb_min->add_value(lb_value);

            lb_mins.push_back(lb_min->copy());
            delete lb_min;
        };
        return lb_mins;
    };

    // static void extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth);

    static void update_estimates(vector<Simplex*> simpls, vector<Function*> funcs, vector<Point*> pareto_front, int iteration);

    double find_simplex_gradient_norm(int crit_id, SimplexGradientStrategy simplex_gradient_strategy){ 
        double L_estimate = 0;
        // Eigen::VectorXd f_diff(_D);
        // Eigen::MatrixXd x_diff(_D, _D);
        // Eigen::MatrixXd x_diff_inv_T(_D, _D);
        // Eigen::VectorXd grad(_D);
        //
        // if (simplex_gradient_strategy == FFMinVert) {  // Gradient at min vertex
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        // //     for (int i=1; i < D+1; i++) { 
        // //         f_diff(i - 1) = _verts[i]->_values[0] - _min_vert->_values[0];
        // //     }; 
        // //
        // //     for (int i=1; i < D+1; i++) {
        // //         for (int j=0; j < D; j++) {
        // //             x_diff(i-1, j) = _verts[i]->_X[j] - _min_vert->_X[j];
        // //         };
        // //     };
        // };
        // if (simplex_gradient_strategy == FFMaxVert) {  // Gradient at max vertex
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        // //     for (int i=0; i < D; i++) { 
        // //         f_diff(i) = _verts[i]->_values[0] - _max_vert->_values[0];
        // //     }; 
        // //
        // //     for (int i=0; i < D; i++) {
        // //         for (int j=0; j < D; j++) {
        // //             x_diff(i, j) = _verts[i]->_X[j] - _max_vert->_X[j];
        // //         };
        // //     };
        // };
        // if (simplex_gradient_strategy == FFMaxVert) {
        //     // throw "FFMaxVert gradient strategy not implemented yet";
        //     // FFAllVertMean
        // };
        //
        // // Find gradient at lowest point
        // x_diff_inv_T = x_diff.inverse().transpose();
        // for (int i=0; i < _D; i++) {
        //     grad[i] = x_diff_inv_T.row(i).dot(f_diff);
        // };
        //
        // // cout << grad << endl;
        //
        // // Find norm of gradient at _min_vert 
        // for (int i=0; i < grad.size(); i++){
        //     L_estimate += pow(grad[i], 2);
        // };
        // L_estimate = sqrt(L_estimate);

        // Find minimum simplex L and use it if its greater then the estimate
        double simplex_min_L = find_simplex_min_L(crit_id);
        if (simplex_min_L > L_estimate) {
            L_estimate = simplex_min_L;
        };
        return L_estimate;
    };

    double find_simplex_min_L(int crit_id) {  // Finds L using Euclidean (l2norm)
        double dist;
        double f_diff;
        double edge_L;
        double max_edge_L = -numeric_limits<double>::max();
        for (int i=0; i < _verts.size(); i++) {
            for (int j=i+1; j < _verts.size(); j++) {
                f_diff = fabs(_verts[i]->_values[crit_id] - _verts[j]->_values[crit_id]);
                dist = l2norm(_verts[i], _verts[j]);
                edge_L = f_diff / dist;   // Note: maybe dist (division by zero) protection is needed?
                if (edge_L > max_edge_L) {
                    max_edge_L = edge_L;
                };
            };
        };
        return max_edge_L;
    };

    double find_simplex_min_L_l1norm(int crit_id) {  // Finds L using City block (l1norm)
        double dist;
        double f_diff;
        double edge_L;
        double max_edge_L = -numeric_limits<double>::max();
        for (int i=0; i < _verts.size(); i++) {
            for (int j=i+1; j < _verts.size(); j++) {
                f_diff = fabs(_verts[i]->_values[crit_id] - _verts[j]->_values[crit_id]);
                dist = l1norm(_verts[i], _verts[j]);
                edge_L = f_diff / dist;   // Note: maybe dist (division by zero) protection is needed?
                if (edge_L > max_edge_L) {
                    max_edge_L = edge_L;
                };
            };
        };
        return max_edge_L;
    };



    double find_edge_lb_value(Point* _le_v1, Point* _le_v2, double _L) {   // Needs testing
        double dist = l2norm(_le_v1, _le_v2);
        double x1[2] = {0., _le_v1->_values[0]};             // (0, simplex[0][-1]['obj'][0])
        double x2[2] = {dist, _le_v1->_values[0] - _L*dist}; // (dist, simplex[0][-1]['obj'][0] - L*dist)
        double x3[2] = {dist, _le_v2->_values[0]};           // (dist, simplex[1][-1]['obj'][0])
        double x4[2] = {0, _le_v1->_values[0] - _L*dist};    // (0, simplex[1][-1]['obj'][0] - L*dist)

        // 2D line intersection based on  http://mathworld.wolfram.com/Line-LineIntersection.html
        double av[2] = {x2[0] - x1[0], x2[1] - x1[1]};
        double bv[2] = {x4[0] - x3[0], x4[1] - x3[1]};
        double cv[2] = {x3[0] - x1[0], x3[1] - x1[1]};

        // cross_product(v1, v2) = (v1.X * v2.Y) - (v1.Y * v2.X)
        double s = ((cv[0]*bv[1] - cv[1]*bv[0]) * (av[0]*bv[1] - av[1]*bv[0])/ (pow(av[0]*bv[1] - av[1]*bv[0], 2) ));

        double intersection[2] = {x1[0] + (x2[0]-x1[0])*s, x1[1] + (x2[1]-x1[1])*s};
        // X = a(simplex[0][:-1]) + s[0]/float(dist) * (a(simplex[1][:-1]) - a(simplex[0][:-1]))
        // return [list(X) + [s[1]]]
        return intersection[1];
    };

    double find_tolerance(vector<Point*> pareto_front) {
        // _min_lbs - is initialized and contains Ms for each criteria

        // Sekti pareto frontą:  naudojant tentą paprastą algoritmą.
        // rasti atstumą iki pareto fronto.
        // For now find lowest distance to vertex.
        vector<double> M;
        for (int i=0; i < _min_lbs.size(); i++) {
            M.push_back(_min_lbs[i]->_values[0]);
        };

        double min_dist = numeric_limits<double>::max();

        for (int i=0; i < pareto_front.size(); i++) {
            double dist = gtl1norm(pareto_front[i]->_values, M);
            if (dist < min_dist) {
                min_dist = dist;
            };
        };
        // If M is dominated, than do not divide this simplex
        if (min_dist == numeric_limits<double>::max()) {
            return 0;
        };
        return -min_dist;  // minus because in convex_hull we need to find max tolerance
    };

    static bool wont_be_divided(Simplex* s) {
        return !s->_should_be_divided;
    };
    static bool not_in_partition(Simplex* s) {
        return !s->_is_in_partition;
    };

    // static double ascending_min_lb_value(Simplex* s1, Simplex* s2) {
    //     return s1->_min_lb_value < s2->_min_lb_value;
    // };

    static double ascending_tolerance(Simplex* s1, Simplex* s2) {
        return s1->_tolerance < s2->_tolerance;
    };

    static double ascending_diameter(Simplex* s1, Simplex* s2) {
        return s1->_diameter < s2->_diameter;
    };

    // static double compare_metric(Simplex* s1, Simplex* s2) {
    //     return s1->_metric__vert_min_value < s2->_metric__vert_min_value;
    // };

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
        vertex->_simplexes.push_back(this);
    };

    int size() {
        return _verts.size();
    };

    void print(){
        cout << " Simplex   diam:  " << _diameter << "   tol:  " << _tolerance;
        cout << "   local-L:  ";
        for (int i=0; i < _C; i++) {
            cout << _local_Ls[i] << " ";
        };
        cout << "   simpl-L: ";
        for (int i=0; i < _C; i++) {
            cout << _simpl_Ls[i] << " ";
        };
        cout << endl;
        for (int i=0; i < _verts.size(); i++){
            _verts[i]->print();
        };
    };

    static void print(vector<Simplex*> simplexes, string label="Printing simplexes:"){
        cout << label << endl;
        for (int i=0; i < simplexes.size(); i++){
            simplexes[i]->print();
        };
    };

    static void log_front(vector<Point*> pareto_front, vector<Simplex*> simplexes_to_divide) {
       ofstream log_file; 
       log_file.open("log/front.txt");
       log_file.close();
       log_file.open("log/front.txt", ios::app);
       for (int j=0; j < pareto_front.size(); j++) {
           for (int i=0; i < pareto_front[j]->size(); i++) {
                log_file << pareto_front[j]->_X[i] << " ";
           };
           log_file << " -> ";
           for (int i=0; i < pareto_front[j]->_values.size(); i++) {
               log_file << pareto_front[j]->_values[i] << " ";
           };
           log_file << endl;
       };

       log_file << "Tolerances" << endl;
       for (int j=0; j < simplexes_to_divide.size(); j++) {
           for (int i=0; i < simplexes_to_divide[j]->_min_lbs.size(); i++) {
              log_file <<  simplexes_to_divide[j]->_min_lbs[i]->_values[0] << " ";
           };
           log_file << endl;
       };
       log_file.close();
    };


    bool is_in_simplex(Point* p, Simplex* simpl) {
//     '''Checks if point is in the triangle region using Barycentric coordinates:
//     www.farinhansford.com/dianne/teaching/cse470/materials/BarycentricCoords.pdf'''
        // Implement this in 2 dimensions
        // A = det(
        //     [simpl->_verts[0]->_X[0]]
        //         );
        // A = det(a([[t[0][0], t[1][0], t[2][0]], [t[0][1], t[1][1], t[2][1]], [1., 1., 1.]]))
        // A1 = det(a([[p[0], t[1][0], t[2][0]], [p[1], t[1][1], t[2][1]], [1, 1, 1]]))
        // A2 = det(a([[t[0][0], p[0], t[2][0]], [t[0][1], p[1], t[2][1]], [1, 1, 1]]))
        // A3 = det(a([[t[0][0], t[1][0], p[0]], [t[0][1], t[1][1], p[1]], [1, 1, 1]]))
        // u = A1 / A
        // v = A2 / A
        // w = A3 / A
        //
        // // Implement this in 3 and 4 dimensions
        // if ((u >= 0) and (v >= 0) and (w >= 0)) {
        //     return True;
        // };
        // return False;

    };

    static void log_partition(vector<Simplex*> simplexes,
                              vector<Simplex*> selected,
                              vector<Function*> funcs,
                              string label="Partition:",
                              int iteration=0) {
       ofstream log_file; 
       log_file.open("log/partition.txt");
       log_file.close();
       log_file.open("log/partition.txt", ios::app);
       log_file << label << iteration << ":" << endl;


       vector<Simplex*> simpls_near_glob_min;
       for (int i=0; i < simplexes.size(); i++) {
            // How to check weather point is in simplex?

            for (int j=0; j < simplexes[i]->_verts.size(); j++) {
                if (simplexes[i]->_verts[j] == funcs[0]->_x_nearest_to_glob_x) {
                    simpls_near_glob_min.push_back(simplexes[i]);
                };
            };
       };

       double min_dist = numeric_limits<int>::max();

       Simplex* min_dist_simpl;
       for (int i=0; i < simpls_near_glob_min.size(); i++) {
           Point* center = new Point(funcs[0]->_D);
           // Iterate through points
           Simplex* simpl = simpls_near_glob_min[i];
           for (int j=0; j < simpl->_verts.size(); j++) {
               for (int c=0; c < simpl->_D; c++) {
                    center->_X[c] += simpl->_verts[j]->_X[c];
               };
           };
           // Divide by number of verts
           for (int c=0; c < simpl->_D; c++) {
               center->_X[c] /= simpl->_verts.size();
           };
           // Find distance and save smallest
           double dist = l2norm(center, funcs[0]->_glob_x);
           if (dist < min_dist) {
               min_dist = dist;
               min_dist_simpl = simpl;
           };
       };
       // cout << "Min dist is: " << min_dist << endl;

       // cout << simpls_near_glob_min.size() << endl;


       for (int i=0; i < simplexes.size(); i++) {
           for (int j=0; j < simplexes[i]->_verts.size(); j++) {
               for (int k=0; k < simplexes[i]->_verts[j]->size(); k++){
                    log_file << simplexes[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << simplexes[i]->_verts[j]->_values[0]<<"); ";
           };
           if (simplexes[i]->_L_strategy == Self) {
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_tolerance << ")" << endl;
           };
           if (simplexes[i]->_L_strategy == Neighbours) {
                log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_tolerance << ")" << endl;
           };
       };
       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << selected[i]->_verts[j]->_values[0]<<"); ";
           };
           if (selected[i]->_L_strategy == Self) {
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_tolerance << ")" << endl;
           };
           if (selected[i]->_L_strategy == Neighbours) {
               log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_tolerance << ")" << endl;
           };
       };

       log_file << "Wanted:" << endl;
       // for (int i=0; i < simpls_near_glob_min.size(); i++) {
           // for (int j=0; j < simpls_near_glob_min[i]->_verts.size(); j++) {
           //     for (int k=0; k < simpls_near_glob_min[i]->_verts[j]->size(); k++){
           //          log_file << simpls_near_glob_min[i]->_verts[j]->_X[k] << " ";
           //     };
           //     log_file << " (" << simpls_near_glob_min[i]->_verts[j]->_values[0]<<"); ";
           // };
           // if (simpls_near_glob_min[i]->_L_strategy == Self) {
           //     log_file << " ("<< simpls_near_glob_min[i]->_diameter << "," << simpls_near_glob_min[i]->_tolerance << ")" << endl;
           // };
           // if (simpls_near_glob_min[i]->_L_strategy == Neighbours) {
           //     log_file << " ("<< simpls_near_glob_min[i]->_diameter << "," << simpls_near_glob_min[i]->_tolerance << ")" << endl;
           // };
           for (int j=0; j < min_dist_simpl->_verts.size(); j++) {
               for (int k=0; k < min_dist_simpl->_verts[j]->size(); k++){
                    log_file << min_dist_simpl->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << min_dist_simpl->_verts[j]->_values[0]<<"); ";
           };
           if (min_dist_simpl->_L_strategy == Self) {
               log_file << " ("<< min_dist_simpl->_diameter << "," << min_dist_simpl->_tolerance << ")" << endl;
           };
           if (min_dist_simpl->_L_strategy == Neighbours) {
               log_file << " ("<< min_dist_simpl->_diameter << "," << min_dist_simpl->_tolerance << ")" << endl;
           };
       // };

       log_file.close();
    };

    virtual ~Simplex(){
        for (int i=0; i < _min_lbs.size(); i++) {
            delete _min_lbs[i];
        };
        // Remove pointers from vertices to this simplex
        for (int i=0; i < _verts.size(); i++) {
            _verts[i]->_simplexes.erase(remove(_verts[i]->_simplexes.begin(), _verts[i]->_simplexes.end(), this), _verts[i]->_simplexes.end());
        };
        _min_lbs.clear();
        // delete _min_lb;
        _verts.clear();
        _neighbours.clear();
    };  
};
vector<double> Simplex::glob_Ls;
bool Simplex::glob_L_was_updated = false;
double Simplex::max_diameter = numeric_limits<double>::max();
double Simplex::glob_L_coef = 0.4;
double Simplex::local_L_coef = 0.4;


void Point::_neighbours_estimates_should_be_updated() {
    for (int sid=0; sid < _simplexes.size(); sid++) {
        _simplexes[sid]->_should_estimates_be_updated = true;
    };
};


// class SimplexTreeNode {
//     SimplexTreeNode(const SimplexTreeNode& other){}
//     SimplexTreeNode& operator=(const SimplexTreeNode& other){}
// public:                
//     SimplexTreeNode(Simplex* value){
//         _height = 1;
//         _value = value;
//         _left = 0;
//         _right = 0;
//     };
//     int _height;
//     Simplex* _value;
//     SimplexTreeNode* _parent;
//     SimplexTreeNode* _left;
//     SimplexTreeNode* _right;
//
//     void print(){
//         if (_left != 0) { cout << "l"; _left->print(); };
//         cout << _value << "("<< _height << ")";
//         if (_right != 0) { cout << "r"; _right->print(); };
//     };
//
//     virtual ~SimplexTreeNode();
// };
//
// class SimplexTree {  // Binary balancing Simplex tree for storing simplex neighbours 
//                      // Simplex adresses are stored in this tree
//     SimplexTree(const SimplexTree& other){}
//     SimplexTree& operator=(const SimplexTree& other){}
// public:
//     SimplexTree(){
//          _max_simpl_Ls = -numeric_limits<double>::max();
//          _tree_root = 0;
//     };
//     double _max_simpl_Ls;
//     SimplexTreeNode* _tree_root;
//
//     void update_height(SimplexTreeNode* node) {
//         int lh = 0;
//         int rh = 0;
//         if (node->_left != 0) { lh = node->_left->_height; }; 
//         if (node->_right != 0) { rh = node->_right->_height; };
//         if (lh > rh) {
//             node->_height = lh + 1;
//         } else {
//             node->_height = rh + 1;
//         };
//         // Also update all ancestors heights
//         if (node->_parent != 0) {
//             update_height(node->_parent);
//         };
//     };
//     void left_right_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched_node;
//         // node left right  <-  node left right left 
//         diatteched_node = node->_left->_right;
//         node->_left->_right = node->_left->_right->_left;
//         if (node->_left->_right != 0) { node->_left->_right->_parent = node->_left; };
//         // Diatteched left = node->_left
//         diatteched_node->_left = node->_left;
//         node->_left->_parent = diatteched_node;
//         // node left  <-  node left right
//         node->_left = diatteched_node;
//         diatteched_node->_parent = node;
//         // Update heights
//         update_height(node);
//         update_height(diatteched_node);
//         update_height(diatteched_node->_left);
//     };
//     void left_left_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched;
//         diatteched = node->_left;
//         node->_left = node->_left->_right;
//         if (node->_left != 0) { node->_left->_parent = node; };  
//         diatteched->_parent = node->_parent;
//         if (node->_parent != 0) {
//             if (node->_parent->_left == node) {
//                 node->_parent->_left = diatteched;
//             } else {
//                 node->_parent->_right = diatteched;
//             };
//         } else {
//             _tree_root = diatteched;
//         };
//         diatteched->_right = node;
//         node->_parent = diatteched;
//         // Update heights
//         update_height(node);
//         update_height(diatteched);
//     };
//     void right_left_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched_node;
//         // node left right  <-  node left right left 
//         diatteched_node = node->_right->_left;
//         node->_right->_left = node->_right->_left->_right;
//         if (node->_right->_left != 0) { node->_right->_left->_parent = node->_right; };
//         // Diatteched left = node->_left
//         diatteched_node->_right = node->_right;
//         node->_right->_parent = diatteched_node;
//         // node left  <-  node left right
//         node->_right = diatteched_node;
//         diatteched_node->_parent = node;
//         // Update heights
//         update_height(node);
//         update_height(diatteched_node);
//         update_height(diatteched_node->_right);
//     };
//     void right_right_rebalance(SimplexTreeNode* node) {
//         SimplexTreeNode* diatteched;
//         diatteched = node->_right;
//         node->_right = node->_right->_left;
//         if (node->_right != 0) { node->_right->_parent = node; };
//         diatteched->_parent = node->_parent;                       
//         if (node->_parent != 0) {
//             if (node->_parent->_left == node) {
//                 node->_parent->_left = diatteched;
//             } else {
//                 node->_parent->_right = diatteched;
//             };
//         } else {
//             _tree_root = diatteched;
//         };
//         diatteched->_left = node;
//         node->_parent = diatteched;
//         // Update heights
//         update_height(node);
//         update_height(diatteched);
//     };
//
//     void check_if_balanced(SimplexTreeNode* node) {  // Rebalances tree if its not balanced
//         int lh = 0;
//         int rh = 0;
//         int llh = 0;
//         int lrh = 0;
//         int rlh = 0;
//         int rrh = 0;
//         if (node->_left != 0) {
//             lh = node->_left->_height;
//             if (node->_left->_left != 0) { llh = node->_left->_left->_height; };
//             if (node->_left->_right != 0) { lrh = node->_left->_right->_height; };
//         };
//         if (node->_right != 0) {
//             rh = node->_right->_height;
//             if (node->_right->_left != 0) { rlh = node->_right->_left->_height; };
//             if (node->_right->_right != 0) { rrh = node->_right->_right->_height; };
//         };
//         if (abs(rh - lh) > 1) {
//             // Not balanced, so rebalance
//             if (rh > lh) {
//                 if (rrh > rlh) {
//                     right_right_rebalance(node);
//                 } else {
//                     right_left_rebalance(node);
//                     right_right_rebalance(node);
//                 };
//             };
//             if (rh < lh) {
//                 if (llh > lrh) {
//                     left_left_rebalance(node);
//                 } else {
//                     left_right_rebalance(node);
//                     left_left_rebalance(node);
//                 };
//             };
//         };
//         if (node->_parent != 0) {
//             check_if_balanced(node->_parent);
//         };
//     };
//
//     Simplex* add(Simplex* value) {
//         // Warning: how is this working, if simplex but not its value is used in comparisons?
//
//         // Get same point or insert given (if inserted returns 0)
//         SimplexTreeNode* node = _tree_root;
//         double grad_norm = value->_simpl_Ls;
//         if (grad_norm > _max_simpl_Ls) {
//             _max_simpl_Ls = grad_norm;
//         };
//         if (_tree_root == 0) {  // Create first tree node
//             _tree_root = new SimplexTreeNode(value);
//         } else {
//             while (true) {  // Walk through tree
//                 if (value > node->_value) {
//                     if (node->_right == 0) {
//                         node->_right = new SimplexTreeNode(value);
//                         node->_right->_parent = node;
//                         update_height(node->_right);
//                         check_if_balanced(node->_right);
//                         return 0;
//                     };
//                     node = node->_right;
//                 } else if (value < node->_value) {
//                     if (node->_left == 0) {
//                         node->_left = new SimplexTreeNode(value);
//                         node->_left->_parent = node;
//                         update_height(node->_left);
//                         check_if_balanced(node->_left);
//                         return 0;
//                     };
//                     node = node->_left;
//                 } else {
//                     // Node value matches given simplex value
//                     return value;
//                 };
//             };
//         };
//     };
//
//     void print() {
//         _tree_root->print();
//         cout << endl;
//     };
//
//     virtual ~SimplexTree();
//
// };
//
// SimplexTreeNode::~SimplexTreeNode() {
//     if (_left != 0) { delete _left; };
//     if (_right != 0) { delete _right; };
// };
//
// SimplexTree::~SimplexTree() {
//     delete _tree_root;
// };


// void Simplex::extend_region_with_vertex_neighbours(Point* vertex, SimplexTree* region, int depth) {
//     // Recursively adds vertex neighbours to region
//     for (int sid=0; sid < vertex->_simplexes.size(); sid++) {
//         Simplex* simpl = vertex->_simplexes[sid];
//         Simplex* result = region->add(simpl);
//         if (depth != 0) {
//             if (result == 0 && simpl->_is_in_partition) {
//                 for (int vid=0; vid < simpl->_verts.size(); vid++) {
//                     if (simpl->_verts[vid] != vertex) {
//                         extend_region_with_vertex_neighbours(simpl->_verts[vid], region, depth-1);
//                     };
//                 };
//             };
//         };    
//     };
// };

void Simplex::update_estimates(vector<Simplex*> simpls, vector<Function*> funcs, vector<Point*> pareto_front, int iteration) {   // Neighbours strategy - updates estimates
    for (int sid=0; sid < simpls.size(); sid++) {
        if (simpls[sid]->_should_estimates_be_updated or Simplex::glob_L_was_updated) {
            //// Use simplex's \hat{L} as initial max_simpl_Ls value
            vector<double> max_simpl_Ls;
            for (int i=0; i < simpls[sid]->_simpl_Ls.size(); i++) {
                max_simpl_Ls.push_back(simpls[sid]->_simpl_Ls[i]);
            };

            // Find max \hat{L} among neighbours
            for (list<Simplex*>::iterator it=simpls[sid]->_neighbours.begin(); it != simpls[sid]->_neighbours.end(); ++it) {
                for (int i=0; i < funcs.size(); i++) {
                    if ((*it)->_simpl_Ls[i] > max_simpl_Ls[i]) {
                    // if (((*it)->_simpl_Ls[i] > max_simpl_Ls[i]) and (simpls[sid]->_diameter <= (*it)->_diameter)) {
                        max_simpl_Ls[i] = (*it)->_simpl_Ls[i];
                    };
                };
            };

            // Update simplex's L
            for (int i=0; i < funcs.size(); i++) {
                simpls[sid]->_local_Ls[i] = max_simpl_Ls[i];
                simpls[sid]->_local_Ls[i] = Simplex::local_L_coef * simpls[sid]->_local_Ls[i];
                // simpls[sid]->_local_Ls[i] = simpls[sid]->_simpl_Ls[i];
            };

            //// Find accurate lower bound point and value estimates with given precision
            for (int i=0; i < simpls[sid]->_min_lbs.size(); i++) {
                delete simpls[sid]->_min_lbs[i];
            };

            // simpls[sid]->_min_lbs = simpls[sid]->find_accurate_lb_min_estimates(simpls[sid]->_verts, Simplex::glob_Ls);
            simpls[sid]->_min_lbs = simpls[sid]->find_one_vert_lb_mins(simpls[sid], Simplex::glob_Ls);
            // simpls[sid]->_min_lbs = simpls[sid]->find_two_best_vert_lb_mins(simpls[sid], Simplex::glob_Ls);

            simpls[sid]->_tolerance = simpls[sid]->_min_lbs[0]->_values[0];   // simpls[sid]->find_tolerance(pareto_front);

            simpls[sid]->_should_estimates_be_updated = false;
        };
    };
    Simplex::glob_L_was_updated = false;


    //// Update expected improvement
    double L = Simplex::glob_Ls[0];

    double min_min_lb = numeric_limits<double>::max();
    Simplex* s;
    for (int sid=0; sid < simpls.size(); sid++) {
        if (simpls[sid]->_min_lbs[0]->_values[0] < min_min_lb) {
            min_min_lb = simpls[sid]->_min_vert->_values[0] - L * simpls[sid]->_diameter;
            // min_min_lb = simpls[sid]->_min_lbs[0]->_values[0];
            s = simpls[sid];
        };
    };
    funcs[0]->_expected_improvement = funcs[0]->_f_min - min_min_lb;


    // Note: gali būti, kad slope apibrėžimas pas mane netinkamas atmetant
    // simpleksus su epsilon (potencialiai optimalių simpleksų parinkimo metu).
    //     simpls[sid]->_lb = simpls[sid].find_lb();
    // };
};

#endif
