#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include "FuncUC.h"
#include "Asimpl.h"
#include "gkls.h"
#include "rnd_gen.h"
#include <malloc.h>
#include <string.h>


using namespace std;

#define no_argument 0
#define required_argument 1
#define optional_argument 2



int _GKLS_class_D[] = {2, 2, 3, 3, 4, 4, 5, 5}; // Parameters of the 8 default GKLS function classes
double _GKLS_class_global_dists[] = {0.9, 0.9, 0.66, 0.9, 0.66, 0.9, 0.66, 0.66};
double _GKLS_class_global_radiuses[] = {0.2, 0.1, 0.2, 0.2, 0.2, 0.2, 0.3, 0.2};
double _GKLS_class_deltas[] = {1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-6, 1e-7, 1e-7};

int d;                              // Dimension of the optimization problem
double delta_;                      // Parameter used in stopping condition
vector<double> lb;                  // Lower bound of feasible region
vector<double> ub;                  // Upper bound of feasible region
vector< vector<double> > points;    // Points where objective value was evaluated


double get_value(vector<double> point) {  // Evaluates objective function at given point in [lb, ub]^d feasible region
    points.push_back(point);              // And memorizes the point where the evaluation was made
    return GKLS_D_func(&point[0]);
};

bool should_stop(vector<double> point) {  // Checks if stopping condition is satisfied by a given point in [lb, ub]^d 
    for (int i=0; i < d; i++) {
        if (pow(delta_, 1./d) * (ub[i] - lb[i]) < fabs(point[i] - GKLS_minima.local_min[GKLS_glob.num_global_minima][i])) {
            return false;
        };
    };
    return true;
};

// Transforms point's coordinates from [0, 1]^d unit-cube to [lb, ub]^d feasible region
vector<double> transform_from_uc(vector<double> point_uc, vector<double> _lb, vector<double> _ub) {
    vector<double> point;
    for (int i=0; i < _lb.size(); i++){
        point.push_back(point_uc[i] * (_ub[i] - _lb[i]) + _lb[i]);
    };
    return point;
};

double get_value_uc(vector<double> point_uc) {  // Evaluates objective function at a given point in [0, 1]^d feasible region
    return get_value(transform_from_uc(point_uc, lb, ub));
};

bool should_stop_uc(vector<double> point_uc) {  // Checks if stopping condition is satisfied by a given point in [0, 1]^d
    return should_stop(transform_from_uc(point_uc, lb, ub));
};


void output_surface_points(Function* func) {
   ofstream log_file; 
   log_file.open("log/surface.txt");
   log_file.close();
   log_file.open("log/surface.txt", ios::app);

   // for (int k=0; k < simplexes[i]->_verts[j]->size(); k++){
   //      log_file << simplexes[i]->_verts[j]->_X[k] << " ";
   // };
   // log_file << " (" << simplexes[i]->_verts[j]->_value << ");";

    int divide_into = 50;
    for (int i=0; i <= divide_into; i++) {
        for (int j=0; j <= divide_into; j++) {
            double coords[2] =  {i * 1./ divide_into,  j * 1./divide_into };
            Point* tmp_point = func->get(new Point(coords, 2));
            log_file << coords[0] << ", " << coords[1] << ": " << tmp_point->_value << endl;
        };
    };
    log_file.close();
};



int main(int argc, char* argv[]) {  // Minimizes 100 functions from one GKLS class; prints intermediate and summarised results
    // Parse parameters
    const struct option longopts[] = {
        {"func_cls", optional_argument, 0, 'c'},
        {"func_id", optional_argument, 0, 'f'},
        {"func_name", optional_argument, 0, 'n'},
        {"stop_crit", optional_argument, 0, 's'},
        {"task_id", required_argument, 0, 't'},
        {"callback", required_argument, 0, 'b'},
        {"max_duration", optional_argument, 0, 'd'},
        {"max_calls", optional_argument, 0, 'i'},
        {"glob_L", optional_argument, 0, 'g'},
    };
    int cls;
    int fid;
    int task_id;
    char* callback = {'\0'};
    char* func_name = {'\0'};
    char* stop_crit = {'\0'};
    int max_calls = 40000;
    int max_duration = 2*3600;
    double glob_L = numeric_limits<double>::max();

    int opt_id;
    int iarg = 0;
    while(iarg != -1) {
        iarg = getopt_long(argc, argv, "cnsftbdig", longopts, &opt_id);
        switch (iarg) {
            case 'c':
                cls = strtoul(optarg, 0, 0);
                break;
            case 'n':
                func_name = strdup(optarg);
                break;
            case 's':
                stop_crit = strdup(optarg);
                break;
            case 'f':
                fid = strtoul(optarg, 0, 0);
                break;
            case 't':
                task_id = strtoul(optarg, 0, 0);
                break;
            case 'b':
                callback = strdup(optarg);
                break;
            case 'd':
                max_duration = strtoul(optarg, 0, 0);
                break;
            case 'i':
                max_calls = strtoul(optarg, 0, 0);
                break;
            case 'g':
                glob_L = strtod(optarg, '\0');
                break;
        };
    };

    // Check which function is provided and use it


    // int number_of_funcs = 100;
    // double calls[100] = {};                         // Prepare an array to store results for each function
    // for (int i=0; i < number_of_funcs; i++) {       // Initialize the array with very big values, in order to make
    //     calls[i] = numeric_limits<double>::max();   // results bad if any of the values wasnt overwritten
    // };
    //
    // int cls = 2;   // Select GKLS class for which 100 objective functions will be minimized
    //
    // for (int func_id=1; func_id <= 1; func_id++) {  // Iterate through objective functions and minimize them


    d = _GKLS_class_D[cls-1];            // Set the dimension of the problem
    delta_ = _GKLS_class_deltas[cls-1];  // Set stopping condition Delta parameter
    for (int i=0; i < d; i++) {          // Initialize lower and upper bounds
        lb.push_back(-1.);
        ub.push_back(1.);
    };

    // Initialize GKLS function
    assert(GKLS_set_default() == GKLS_OK);
    GKLS_dim = d;
    GKLS_global_dist = _GKLS_class_global_dists[cls-1];
    GKLS_global_radius = _GKLS_class_global_radiuses[cls-1];
    GKLS_num_minima = 10;
    GKLS_global_value = GKLS_GLOBAL_MIN_VALUE;
    assert(GKLS_domain_alloc() == GKLS_OK);
    assert(GKLS_parameters_check() == GKLS_OK);
    assert(GKLS_arg_generate(fid) == GKLS_OK);
    assert(GKLS_minima.f[GKLS_glob.num_global_minima] == GKLS_global_value);

    // cout << GKLS_minima.x[1] << endl;
    // cout << "Glob min: ";
    // for (int k=0; k < 2; k++) {
    //     cout << (GKLS_minima.local_min[GKLS_glob.num_global_minima][k] - -1.) / (1. - -1.) << " " << endl;
    // };
    // cout << endl;

    // Create an object of a function, which is defined over a unit-cube [0,1]^d
    // This object can evaluate objective function values only using get_value_uc() method
    // and can check weather a point satisfies stopping condition using only should_stop_uc() method
    // These two methods are defined above and are passed as arguments to the function constructor,
    // all other global variables, methods of this file are not accessible inside the function object
    Function* func_uc = new FuncUC(d, get_value_uc, should_stop_uc);

    // output_surface_points(func_uc);
    // delete func_uc;  // Clear allocated memory
    // func_uc = new FuncUC(d, get_value_uc, should_stop_uc);


    // Minimize the function using Asimpl algorithm (it has alpha parameter set to 0.4)

    Asimpl* alg = new Asimpl(max_calls, max_duration);
    alg->minimize(func_uc);
    cout << "cls: " << cls << " fid: " << fid << " calls: " << points.size()
               << " --calls=" << points.size()
               << " --subregions=" << alg->_partition.size()
               << " --duration=" << alg->_duration
               << " --status=" << alg->_status
               << " --f_min=" << func_uc->_f_min
               << " --global_L=" << Simplex::glob_L
               << " --min_diam=" << alg->_partition[0]->_diameter
               << " --max_diam=" << alg->_partition[alg->_partition.size() - 1]->_diameter << endl;


    // Save results
    if (callback != '\0') {
        string cmd;
        stringstream cmd_ss; 
        cmd_ss.precision(10);
        cmd_ss << callback
               << " --calls=" << points.size()
               << " --subregions=" << alg->_partition.size()
               << " --duration=" << alg->_duration
               << " --task_id=" << task_id
               << " --status=" << alg->_status
               << " --f_min=" << func_uc->_f_min
               << " --global_L=" << Simplex::glob_L
               << " --min_diam=" << alg->_partition[0]->_diameter
               << " --max_diam=" << alg->_partition[alg->_partition.size() - 1]->_diameter
               << " -exe=" << argv[0] << endl;
        cmd = cmd_ss.str();
        popen(cmd.c_str(), "r");
    };

    // calls[fid-1] = points.size();  // Memorize and print the number of function evaluations
    // cout << "cls: " << cls << ", fid: " << fid << ", calls: " << points.size() << endl;

    delete func_uc;  // Clear allocated memory
    delete alg;
    points.clear();
    GKLS_free();

    // };
    //
    // // Print summarised results (for 100 functions)
    // sort(calls, calls + number_of_funcs);      // Sort numbers of function evaluations ascending
    // double calls_sum = 0;                      // Calculate total function evaluations to find the average
    // for (int i=0; i < number_of_funcs; i++) {
    //     calls_sum += calls[i];
    // };
    // cout << "Average: " << calls_sum/100. << "  Median: " << (calls[49] + calls[50])/2. << " Largest: " << calls[99] << endl;

    return 0;
};
