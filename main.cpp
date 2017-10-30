#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include "functions.h"
#include "Libre.h"
#include "gkls.h"
#include "rnd_gen.h"
#include <malloc.h>
#include <string.h>


using namespace std;

#define no_argument 0
#define required_argument 1
#define optional_argument 2


int main(int argc, char* argv[]) {
    // Parse parameters
    const struct option longopts[] = {
        {"func_name", optional_argument, 0, 'n'},
        {"alpha", optional_argument, 0, 'a'},
        {"task_id", required_argument, 0, 't'},
        {"callback", required_argument, 0, 'b'},
        {"max_calls", optional_argument, 0, 'i'},
        {"max_duration", optional_argument, 0, 'd'},
        {"d", optional_argument, 0, 'l'},
    };

    char* func_name = {'\0'};
    double alpha = numeric_limits<double>::max();
    int task_id;
    int D = 3;
    char* callback = {'\0'};
    int max_calls = 100000;
    int max_duration = 4*3600;

    int opt_id;
    int iarg = 0;
    while(iarg != -1) {
        iarg = getopt_long(argc, argv, "antbdil", longopts, &opt_id);
        switch (iarg) {
            case 'n':
                func_name = strdup(optarg);
                break;
            case 'a':
                alpha = strtod(optarg, '\0');
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
            case 'l':
                D = strtoul(optarg, 0, 0);
                break;
        };
    };

    if (alpha == numeric_limits<double>::max()) {
        alpha = 0.4;
    };

    Function::_filename = string(argv[0]);
    Function* func = get_function(func_name, D);
    FunctionUC* func_uc = new FunctionUC(func);

    Libre* alg = new Libre(max_calls, max_duration, alpha);
    alg->minimize(func_uc);

    cout << "Trials for " << func_uc->_name << ": " << func_uc->_evaluations << endl
         << "   Points in Pareto front: " << func_uc->_pareto_front.size() << endl
         << "   Hyper-volume: " << func_uc->hyper_volume() << endl
         << "   Uniformity: " << func_uc->uniformity()
         << "   Duration: " << alg->_duration
         << endl;

    // func_uc->show_pareto_front();

    //// Save results
    if (callback != '\0') {
        string cmd;
        stringstream cmd_ss;
        cmd_ss.precision(10);
        cmd_ss << callback
               << " --calls=" << func_uc->_evaluations
               << " --points_in_pareto_front=" << func_uc->_pareto_front.size()
               << " --hyper_volume=" << func_uc->hyper_volume()
               << " --uniformity=" << func_uc->uniformity()
               << " --stats_file=" << func_uc->_stats_log_filename
               << " --subregions=" << alg->_partition.size()
               << " --duration=" << alg->_duration
               << " --task_id=" << task_id
               << " --status=" << alg->_status
               << " --min_diam=" << alg->_partition[0]->_diameter
               << " --max_diam=" << alg->_partition[alg->_partition.size() - 1]->_diameter
               << " -exe=" << argv[0] << endl;
        cmd = cmd_ss.str();
        popen(cmd.c_str(), "r");
    };

    delete func_uc;  // Clear allocated memory
    delete func;
    delete alg;

    return 0;
};
