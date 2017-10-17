/* Copyright Albertas Gimbutas 2017, all rights reserved */
#ifndef LIBRE_H
#define LIBRE_H
#include <math.h>
#include "utils.h"
#include <ctime>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <fstream>
#include <sstream>


using namespace std;


class Libre {
    Libre(const Libre& other) {};
    Libre& operator=(const Libre& other) {};
public:
    Libre(int max_calls, double max_duration, double alpha) {
        _iteration = 0;
        Simplex::alpha = alpha;

        ofstream log_file;
        log_file.open("log/partition.txt");
        log_file.close();
        log_file.open("log/front.txt");
        log_file.close();

        _max_calls = max_calls;
        _max_duration = max_duration;
        _duration = 0;
    };

    vector<Simplex*> _partition;
    FunctionUC* _func;
    int _iteration;
    string _status;
    double _duration;   // Duration in seconds
    double _max_duration;   // Maximum allowed duration of minimization in seconds
    int _max_calls;

    vector<Simplex*> partition_unit_cube_into_simplices_combinatoricly(int n) {
        // Partitions n-unit-cube into simplices using combinatoric vertex triangulation algorithm
        vector<Simplex*> partition;
        int number_of_simpleces = 1;
        for (int i = 1; i <= n; i++) {
            number_of_simpleces *= i;
        };

        int teta[n];
        for (int i=0; i < n; i++){
            teta[i] = i;
        };

        int test_i = 0;

        do {
            test_i += 1;

            int triangle[n+1][n];

            for (int k = 0; k < n; k++) {
                triangle[0][k] = 0;
            };

            for (int vertex=0; vertex < n; vertex++) {
                for (int j = 0; j < n + 1; j++) {
                    triangle[vertex + 1][j] = triangle[vertex][j];
                };
                triangle[vertex + 1][teta[vertex]] = 1;
            };

            Simplex* simpl = new Simplex();
            for (int i=0; i < n + 1; i++) {
                Point* tmp_point = new Point(triangle[i], n);

                Point* point = _func->get_point_with_values(tmp_point);
                if (tmp_point != point) {
                    delete tmp_point;
                };
                simpl->add_vertex(point);
            };
            simpl->init_parameters(_func);
            partition.push_back(simpl);

        } while (next_permutation(teta, teta+n));
        return partition;
    };

    vector<Simplex*> convex_hull(vector<Simplex*> simplices) {
        int m = simplices.size() - 1;
        if (m <= 1) { return simplices; };
        int START = 0;
        int v = START;
        int w = m;
        bool flag = false;
        bool leftturn = false;
        int a, b, c;
        double det_val;
        while ((nextv(v, m) != START) or (flag == false)) {
            if (nextv(v, m) == w) {
                flag = true;
            }
            a = v;
            b = nextv(v, m);
            c = nextv(nextv(v, m), m);

            double* matrix[3];
            double line1[3] = {simplices[a]->_diameter, simplices[a]->_tolerance, 1.};
            double line2[3] = {simplices[b]->_diameter, simplices[b]->_tolerance, 1.};
            double line3[3] = {simplices[c]->_diameter, simplices[c]->_tolerance, 1.};
            matrix[0] = line1;
            matrix[1] = line2;
            matrix[2] = line3;
            det_val = Determinant(matrix, 3);

            if (det_val >= 0){
                leftturn = 1;
            } else {
                leftturn = 0;
            };
            if (leftturn) {
                v = nextv(v, m);
            } else {
                simplices.erase(simplices.begin() + nextv(v, m));
                m -= 1;
                w -= 1;
                v = predv(v, m);
            };
        };
        return simplices;
    };

    int nextv(int v, int m) {
        if (v == m) {
            return 0;
        };
        return v + 1;
    };

    int predv(int v, int m) {
        if (v == 0) {
            return m;
        };
        return v - 1;
    };

    vector<Simplex*> select_simplices_to_divide() {
        // Note: tolerance value in simplices is negative, to be able to use standard convex-hull method
        vector<Simplex*> selected_simplices;

        // Sort simplices ascending by their diameter
        vector<Simplex*> sorted_partition = _partition;
        sort(sorted_partition.begin(), sorted_partition.end(), Simplex::ascending_diameter);

        // Find simplices with best tolerance values and unique diameters
        Simplex* max_tolerance_simplex = sorted_partition[0];  // Initial value
        vector<double> diameters;
        vector<Simplex*> best_for_size;

        bool unique_diameter;
        bool found_with_same_size;
        for (int i=0; i < sorted_partition.size(); i++) {
            if (sorted_partition[i]->_tolerance < max_tolerance_simplex->_tolerance) {
                max_tolerance_simplex = sorted_partition[i];
            };
            // Saves unique diameters
            unique_diameter = true;
            for (int j=0; j < diameters.size(); j++) {
                if (diameters[j] == sorted_partition[i]->_diameter) {
                    unique_diameter = false; break;
                };
            };
            if (unique_diameter) {
                diameters.push_back(sorted_partition[i]->_diameter);
            };

            // If this simplex is better then previous with same size swap them
            found_with_same_size = false;
            for (int j=0; j < best_for_size.size(); j++) {
                if (best_for_size[j]->_diameter == sorted_partition[i]->_diameter){
                    found_with_same_size = true;
                    if (best_for_size[j]->_tolerance > sorted_partition[i]->_tolerance) {
                        best_for_size.erase(best_for_size.begin()+j);
                        best_for_size.push_back(sorted_partition[i]);
                    };
                };
            };
            if (!found_with_same_size) {
                best_for_size.push_back(sorted_partition[i]);
            };
        };

        // Find strict pareto optimal solutions using convex-hull strategy
        vector<Simplex*> selected;
        if (max_tolerance_simplex == best_for_size[best_for_size.size()-1]) {
            selected.push_back(max_tolerance_simplex);
        } else {
            if ((best_for_size.size() > 2) && (max_tolerance_simplex != best_for_size[best_for_size.size()-1])) {
                vector<Simplex*> simplices_below_line;
                double a1 = max_tolerance_simplex->_diameter;  // Should be like this based on Direct Matlab implementation
                double b1 = max_tolerance_simplex->_tolerance;
                double a2 = best_for_size[best_for_size.size()-1]->_diameter;
                double b2 = best_for_size[best_for_size.size()-1]->_tolerance;

                double slope = (b2 - b1)/(a2 - a1);
                double bias = b1 - slope * a1;

                for (int i=0; i < best_for_size.size(); i++) {
                    if (best_for_size[i]->_diameter >= a1) {  // Dont take into consideration smallel diameter simplices
                        if (best_for_size[i]->_tolerance < slope*best_for_size[i]->_diameter + bias +1e-12) {
                            simplices_below_line.push_back(best_for_size[i]);
                        };
                    };
                };
                selected = convex_hull(simplices_below_line);
            } else {
                selected = best_for_size;
            };
        };

        for (int i=0; i < selected.size(); i++) {
            selected[i]->_should_be_divided = true;
        };

        // Remove simplices which were not selected and should not be divided
        selected.erase(remove_if(selected.begin(), selected.end(), Simplex::wont_be_divided), selected.end());

        // Select all simplices which have best min_lb for its size
        for (int i=0; i < sorted_partition.size(); i++) {
            for (int j=0; j < selected.size(); j++) {
                if ((sorted_partition[i]->_diameter == selected[j]->_diameter) &&
                    (sorted_partition[i]->_tolerance == selected[j]->_tolerance)) {
                    selected_simplices.push_back(sorted_partition[i]);
                };
            };
        };
        return selected_simplices;
    };

    vector<Simplex*> divide_simplex(Simplex* simplex) {
        vector<Simplex*> divided_simplices;

        // Find longest edge middle point
        int n = _func->_D;
        double c[n];
        for (int i=0; i < n; i++) {
            c[i] = (simplex->_le_v1->_X[i] + simplex->_le_v2->_X[i]) / 2.;
        };
        Point* tmp_point = new Point(c, n);
        Point* middle_point = _func->get_point_with_values(tmp_point);
        if (tmp_point != middle_point) {
            delete tmp_point;
        };

        // Construct two new simplices using this middle point.
        Simplex* left_simplex = new Simplex();
        Simplex* right_simplex = new Simplex();

        for (int i=0; i < simplex->size(); i++) {
            if (simplex->_verts[i] != simplex->_le_v1){
                right_simplex->add_vertex(simplex->_verts[i]);
            } else {
                right_simplex->add_vertex(middle_point);
            };
            if (simplex->_verts[i] != simplex->_le_v2) {
                left_simplex->add_vertex(simplex->_verts[i]);
            } else {
                left_simplex->add_vertex(middle_point);
            };
        };

        left_simplex->init_parameters(_func);
        right_simplex->init_parameters(_func);

        simplex->_is_in_partition = false;

        divided_simplices.push_back(left_simplex);
        divided_simplices.push_back(right_simplex);
        return divided_simplices;
    };

    vector<Simplex*> divide_simplices(vector<Simplex*> simplices) {
        vector<Simplex*> new_simplices;
        for (int i=0; i < simplices.size(); i++) {
            vector<Simplex*> divided_simplices = divide_simplex(simplices[i]);

            for (int j=0; j < divided_simplices.size(); j++) {
                new_simplices.push_back(divided_simplices[j]);
            };
        };
        return new_simplices;
    };

    void show_partition(vector<Simplex*> selected) {
        Simplex::log_partition(_partition, selected, _func);
        FILE* testp = popen("python log/show_partition.py log/partition.txt", "r");
        //// FILE* testp = popen("python log/img.py", "r");
        pclose(testp);
    };

    void minimize(FunctionUC* func){
        Simplex::reset_glob_Ls(func->_C);

        timestamp_t start = get_timestamp();
        _func = func;
        _partition = partition_unit_cube_into_simplices_combinatoricly(_func->_D);

        sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);
        Simplex::update_tolerance_values(_partition, _func);

        while (_func->_evaluations < 480) {   // (!_func->is_accurate_enough()) {
                                              // Should be accurate enough:
                                              // can be defined based on tolerance.
                                              // need to define methodology for experiments comparing genetic algorithms.

            /////////////////////////////////////////////
            ////    Different stopping conditions    ////
            /////////////////////////////////////////////


            // Select simplices to divide
            vector<Simplex*> simplices_to_divide;
            if (_iteration == 0) {
                simplices_to_divide = _partition;
            } else {
                simplices_to_divide = select_simplices_to_divide();
            };

            // Divide seletected simplices
            vector<Simplex*> new_simplices = divide_simplices(simplices_to_divide);

            // Remove partitioned simplices
            _partition.erase(remove_if(_partition.begin(), _partition.end(), Simplex::not_in_partition), _partition.end());
            for (int i=0; i < simplices_to_divide.size(); i++) {
                delete simplices_to_divide[i];
            };
            simplices_to_divide.clear();

            // Add new simplices to _partition
            for (int i=0; i < new_simplices.size(); i++) {
                _partition.push_back(new_simplices[i]);
            };

            sort(_partition.begin(), _partition.end(), Simplex::ascending_diameter);
            Simplex::update_tolerance_values(_partition, _func);
            _iteration += 1;

            timestamp_t end = get_timestamp();
            _duration = (end - start) / 1000000.0L;
        };

        // Metrics should be added and tested (compare results with Fortran implementation).

        if ((_func->_evaluations <= _max_calls) && (_duration <= _max_duration)) {
            _status = "D";
        } else {
            _status = "S";
        };
    };

    virtual ~Libre(){
        for (int i=0; i < _partition.size(); i++) {
            delete _partition[i];
        };
        _partition.clear();
    };
};

#endif
