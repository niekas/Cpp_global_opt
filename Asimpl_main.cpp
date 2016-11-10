#include <iostream>
#include <getopt.h>
// getopt example: http://stackoverflow.com/questions/8793020/using-getopt-long-c-how-do-i-code-up-a-long-short-option-to-both-require-a
#include "Asimpl.h"
#include <stdio.h>
#include <fstream>

#define no_argument 0
#define required_argument 1
#define optional_argument 2


using namespace std;

int main(int argc, char* argv[]) {
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
    vector<Function*> funcs;
    if (func_name != '\0') {
        Function* func = get_function(func_name, stop_crit);
        funcs.push_back(func);
    } else {
        funcs.push_back(new GKLSFunction(cls, fid));
        if (stop_crit != '\0') {
            funcs[0]->_stopping_criteria = stop_crit;
        };
    };

    Asimpl* alg;
    alg = new Asimpl(max_calls, max_duration);

    if (glob_L != numeric_limits<double>::max()) {
        Simplex::glob_Ls.push_back(glob_L);
    };

    alg->minimize(funcs);

    // Print results
    double X[funcs[0]->_D];
    for (int i=0; i < funcs[0]->_D; i++) {
        X[i] = funcs[0]->transform_back(funcs[0]->_glob_x, i);
    };
    cout << "Glob value: " << funcs[0]->value(new Point(X, funcs[0]->_D)) << " == " << funcs[0]->_glob_f << endl;
    cout.precision(12);
    // cout << "Cls: " << cls << "  Fid: " << fid << endl;
    // cout << "Name: " << func_name << "  Stop_crit: " << stop_crit<< endl;
    if (alg->_status == "S") { cout << "  -->> Suspended <<--" << endl; }
    cout << "Calls: " << funcs[0]->_calls 
         << ", status: " << alg->_status  
         << ", duration: " << alg->_duration  
         << ", subregions: " << alg->_partition.size() << endl;  

    for (int i=0; i < funcs.size(); i++) {
        cout.precision(10);
        cout << "Solution for criteria " << i + 1 << ": " << funcs[i]->_f_min << endl;
        // funcs[i]->_x_nearest_to_glob_x->print();
        cout << "   Global minima for criteria " << i + 1 << ": " << funcs[i]->_glob_f << endl;
        funcs[i]->_glob_x->print();
    };


    // Save results
    if (callback != '\0') {
        string cmd;
        stringstream cmd_ss; 
        cmd_ss.precision(10);
        cmd_ss << callback
               << " --calls=" << funcs[0]->_calls
               << " --subregions=" << alg->_partition.size()
               << " --duration=" << alg->_duration
               << " --task_id=" << task_id
               << " --status=" << alg->_status
               << " --x_min=" << *funcs[0]->_x_min
               << " --f_min=" << funcs[0]->_f_min
               << " --global_L=" << Simplex::glob_Ls[0]
               << " --min_diam=" << alg->_partition[0]->_diameter
               << " --max_diam=" << alg->_partition[alg->_partition.size() - 1]->_diameter
               << " -exe=" << argv[0] << endl;
        cmd = cmd_ss.str();
        popen(cmd.c_str(), "r");
    };

    // Free memory
    delete alg;

    for (int i=0; i < funcs.size(); i++) {
        delete funcs[i];
    };
    funcs.clear();
    return 0;
};
