#include <iostream>
#include "Disimplv.h"

using namespace std;


void print_stats(int calls[], int subregions[], int n) {
    sort(calls, calls+n);
    sort(subregions, subregions+n);
    int calls_sum = 0;
    for (int i=0; i < n; i++) {
        calls_sum += calls[i];
    };
    cout << "Calls50: " << calls[49] << " Calls100: " << calls[99] << " Average: " << calls_sum/100. <<
        " Subregions50: " << subregions[49] << " Subregions100: " << subregions[99] << endl; 
};

int main() {
    GKLSFunction* func;
    Disimplv* alg;
    int n = 100;
    for (int cls=1; cls <= 8; cls++) {
        int calls[100];
        int subregions[100];
        for (int fid=1; fid <= n; fid++) {
            alg = new Disimplv(1.0, 1000000);
            func = new GKLSFunction(cls, fid);

            alg->minimize(func);
            cout << (*alg) << endl;

            calls[fid-1] = func->_calls;
            subregions[fid-1] = alg->_partition.size();

            delete alg;
            delete func;
            // break;
        };
        cout << "Class " << cls << " "; 
        print_stats(calls, subregions, n);
    };
    return 0;
};
