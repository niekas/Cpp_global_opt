#include <stdlib.h>
#include <iostream>
#include <getopt.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <sstream>
#include <vector>

using namespace std;


int main(int argc, char* argv[]) {
    FILE* p = popen("python problems/genetic.py zdt1 0.5 0.5 0.5 0.5 0.5 0.5", "r");
    char output[10000];
    if (p != NULL) {
        while(fgets(output, sizeof(output), p) != NULL) {
            // cout << output << endl;
        };
    };
    string output_str = string(output);
    stringstream ss(output_str);

    double num;
    vector<double> values;
    while (!ss.eof()) {
        ss >> num;
        values.push_back(num);
    };
    cout.precision(16);
    cout << values.size() << ", " << values[0] << ", "<< values[1] << endl;

    pclose(p);
    return 0;
};
