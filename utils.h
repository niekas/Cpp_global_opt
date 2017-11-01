/* Copyright Albertas Gimbutas 2017, all rights reserved */
#ifndef SIMPLEX_H
#define SIMPLEX_H
#include <fstream>
#include <sstream>
#include <sys/time.h>

using namespace std;


double l2norm(Point* p1, Point* p2) {
    // Finds Euclidean distance between two points
    double squared_sum = 0;
    for (int i=0; i < p1->size(); i++){
        squared_sum += pow(p1->_X[i] - p2->_X[i], 2);
    };
    return sqrt(squared_sum);
};

double l2norm(vector<double> p1, vector<double> p2) {
    // Finds Euclidean distance between two points
    double squared_sum = 0;
    for (int i=0; i < p1.size(); i++){
        squared_sum += pow(p1[i] - p2[i], 2);
    };
    return sqrt(squared_sum);
};

double gtl2norm(vector<double> p1, vector<double> p2) {
    // Finds Positive distance between two points
    double squared_sum = 0;
    for (int i=0; i < p1.size(); i++){
        if ((p1[i] - p2[i]) > 0) {
            squared_sum += pow(p1[i] - p2[i], 2);
        };
    };
    return sqrt(squared_sum);
};


typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp() {
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_usec + (timestamp_t) now.tv_sec * 1000000;
};

class Simplex {
    Simplex(const Simplex& other){}
    Simplex& operator=(const Simplex& other){}
public:
    Simplex() {
        _D = 0;
        _C = 0;
        _le_v1 = 0;
        _le_v2 = 0;
        _diameter = 0;
        _tolerance = -numeric_limits<double>::max();

        _is_in_partition = true;
        _should_be_divided = false;
        _should_tolerance_be_updated = true;
    };

    int _D;             // Dimension of variable space
    int _C;             // Number of criterias

    Point* _le_v1;      // First vertex of the longest edge
    Point* _le_v2;      // Second vertex of the longest edge
    double _diameter;   // Longest edge length

    vector<double> _min_Ls;      // Minimum L for this simplex
    vector<double> _min_lbs;     // Minimum of the lower bound over simplex found using vertex with the lowest value
    double _tolerance;
    bool _should_tolerance_be_updated;
    Point* _pf_tolerance_point;

    static double alpha;              // Coeficient of search globality
    static vector<double> glob_Ls;             // Globally known biggest min L

    static void reset_glob_Ls(int C) {
        glob_Ls.clear();
        // Todo: test weather clear actually removes all the elements from vector
        for (int i=0; i < C; i++) {
            glob_Ls.push_back(numeric_limits<double>::max());
            glob_Ls_was_updated.push_back(true);
        };
    };


    static vector<bool> glob_Ls_was_updated;

    vector<Point*> _verts;    // Vertices of this simplex (points with coordinates and values)
    bool _is_in_partition;    // Is this simplex in the current partition
    bool _should_be_divided;  // Should this simplex be divided in next iteration
    bool _should_lb_mins_be_updated;   // Should the minimums of Lipschitz lower bound be updated

    vector<Point*> _min_verts;     // Vertex with lowest function values

    void print() {
        cout << " Simplex   diam:  " << _diameter << "   tol:  " << _tolerance;

        cout << "   min_Ls:  ";
        for (int c=0; c < _C; c++) {
            cout << _min_Ls[c];
            if (c != _C-1) {
                cout << ", ";
            };
        };
        cout <<  endl;

        for (int i=0; i < _verts.size(); i++){
            _verts[i]->print();
        };
    };

    static void log_partition(vector<Simplex*> simplexes,
                              vector<Simplex*> selected,
                              FunctionUC* func,
                              string label="Partition:",
                              int iteration=0) {
       ofstream log_file;
       log_file.open("log/partition.txt");
       log_file.close();
       log_file.open("log/partition.txt", ios::app);
       log_file << label << iteration << ":" << endl;

       for (int i=0; i < simplexes.size(); i++) {
           for (int j=0; j < simplexes[i]->_verts.size(); j++) {
               for (int k=0; k < simplexes[i]->_verts[j]->size(); k++){
                    log_file << simplexes[i]->_verts[j]->_X[k] << " ";
               };
               log_file << "(";
               for (int c=0; c < simplexes[i]->_C; c++) {
                   log_file << simplexes[i]->_verts[j]->_values[c];
                   if (c != simplexes[i]->_C - 1) {
                       log_file << " ";
                   };
               };
               log_file << "); ";
           };
           log_file << "("<< simplexes[i]->_diameter << "," << simplexes[i]->_tolerance << ")" << endl;
       };

       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (";
               for (int c=0; c < selected[i]->_C; c++) {
                   log_file << selected[i]->_verts[j]->_values[c];
                   if (c != selected[i]->_C - 1) {
                       log_file << " ";
                   };
               };
               log_file << ");";
           };
           log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_tolerance << ")" << endl;
       };

       log_file.close();
    };


    void init_parameters(FunctionUC* func) {   // Called when all verts have been added
        _D = _verts.size() - 1;
        _C = func->_C;

        // Sorts vertexes ascending by their function value
        sort(_verts.begin(), _verts.end(), Point::ascending_value);

        // Find longest edge length (simplex diameter) and its vertices
        double edge_length;  // Temporary variable
        for (int a=0; a < _verts.size(); a++) {
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

        // Find minimum L for ths simplex
        _min_Ls = find_simplex_min_Ls();

        // Initialize or update global L if needed
        for (int i=0; i < Simplex::glob_Ls.size(); i++) {
            if (Simplex::glob_Ls[i] == numeric_limits<double>::max()) {
                Simplex::glob_Ls[i] = Simplex::alpha * _min_Ls[i];
                Simplex::glob_Ls_was_updated[i] = true;
            } else {
                if (Simplex::glob_Ls[i] < Simplex::alpha * _min_Ls[i]) {
                    Simplex::glob_Ls[i] = Simplex::alpha * _min_Ls[i];
                    Simplex::glob_Ls_was_updated[i] = true;
                };
            };
        };

        // Find vertices with minimum function values
        Point* min_vert;
        for (int c=0; c < _C; c++) {
            min_vert = _verts[0];
            for (int i=0; i < _verts.size(); i++) {
                if (_verts[i]->_values[c] < min_vert->_values[c]) {
                    min_vert = _verts[i];
                };
            };
            _min_verts.push_back(min_vert);
        };
    };

    vector<double> find_min_vert_lb_mins(Simplex* simpl, vector<double> Ls) {
        vector<double> min_lbs;
        for (int c=0; c < simpl->_C; c++) {
            min_lbs.push_back(_min_verts[c]->_values[c] - Ls[c] * simpl->_diameter);
        };
        // Finds minimum of lower bound, which is constructed from the vertex with lowest function value
        return min_lbs;  // _min_vert->_value - L * simpl->_diameter;
    };

    vector<double> find_simplex_min_Ls() {
        // Finds minimum L for this simplex by finding min L for each simplex edge

        vector<double> max_edge_Ls;
        double dist;
        double f_diff;
        double edge_L;

        for (int c=0; c < _C; c++) {
            max_edge_Ls.push_back(-numeric_limits<double>::max());

            for (int i=0; i < _verts.size(); i++) {
                for (int j=i+1; j < _verts.size(); j++) {
                    f_diff = fabs(_verts[i]->_values[c] - _verts[j]->_values[c]);
                    dist = l2norm(_verts[i], _verts[j]);
                    edge_L = f_diff / dist;         // Note: maybe dist (division by zero) protection is needed?
                    if (edge_L > max_edge_Ls[c]) {  //       Practically this case does not occur.
                        max_edge_Ls[c] = edge_L;
                    };
                };
            };

        };

        return max_edge_Ls;
    };

    bool dominates(vector<double> y, vector<double> r) {
        // Does y dominate r;  True if all y values are better than r
        bool all_better = true;
        for (int i=0; i < r.size(); i++) {
            if (y[i] >= r[i]) {
                all_better = false;
                break;
            };
        };
        return all_better;
    };

    double find_tolerance(vector<Point*> pareto_front) {
        // Find R
        vector<double> R;
        for (int k=0; k < _C; k++) {
            R.push_back(_min_lbs[k]);
        };

        // Is R dominated?
        bool r_dominated = false;
        for (int i=0; i < pareto_front.size(); i++) {
            if (dominates(pareto_front[i]->_values, R) == true) {
                r_dominated = true;
            };
        };

        double min_dist = numeric_limits<double>::max();
        double dist;
        if (r_dominated == true) {
            // Find smallest gtl2norm, when R is dominated.
            for (int i=0; i < pareto_front.size(); i++) {
                dist = gtl2norm(R, pareto_front[i]->_values);
                if (dist < min_dist) {
                    min_dist = dist;
                    _pf_tolerance_point = pareto_front[i];
                };
            };
            min_dist = -1 * min_dist;
        } else {
            // Find smallest gtl2norm, when R is not dominated.
            for (int i=0; i < pareto_front.size(); i++) {
                dist = gtl2norm(pareto_front[i]->_values, R);
                if (dist < min_dist) {
                    min_dist = dist;
                    _pf_tolerance_point = pareto_front[i];
                };
            };
        };
        return -min_dist;    // minus because in convex_hull we need to find max tolerance
    };


    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
        vertex->_simplices.push_back(this);
    };

    int size() {
        return _verts.size();
    };

    // static void update_min_lb_values(vector<Simplex*> simpls, FunctionUC* func);
    static void update_tolerance_values(vector<Simplex*> simpls, FunctionUC* func);

    static bool wont_be_divided(Simplex* s) {
        return !s->_should_be_divided;
    };

    static bool not_in_partition(Simplex* s) {
        return !s->_is_in_partition;
    };

    static double ascending_diameter(Simplex* s1, Simplex* s2) {
        return s1->_diameter < s2->_diameter;
    };

    virtual ~Simplex(){
        // for (int i=0; i < _verts.size(); i++) {
        //     _verts[i]->_simplices.erase(remove(_verts[i]->_simplices.begin(), _verts[i]->_simplices.end(), this), _verts[i]->_simplices.end());
        // };
        _verts.clear();
    };
};
// bool Simplex::glob_L_was_updated = false;
// double Simplex::glob_L = numeric_limits<double>::max();
// bool Simplex::glob_Ls_was_updated
vector<bool> glob_Ls_was_updated_vector;
vector<bool> Simplex::glob_Ls_was_updated = glob_Ls_was_updated_vector;
vector<double> glob_Ls_vector;
vector<double> Simplex::glob_Ls = glob_Ls_vector;
double Simplex::alpha = numeric_limits<double>::max();


// void Simplex::update_min_lb_values(vector<Simplex*> simpls, FunctionUC* func) {
//     // For each criteria update tolerances if neeeded.            Should only update Ls and tolerances for criterias simplices which were
//
//     for (int sid=0; sid < simpls.size(); sid++) {
//         if (simpls[sid]->_should_lb_mins_be_updated or Simplex::glob_L_was_updated) {
//             simpls[sid]->_min_lb = simpls[sid]->find_min_vert_lb_min(simpls[sid], Simplex::glob_L);
//             simpls[sid]->_should_lb_mins_be_updated = false;
//         };
//     };
//     Simplex::glob_L_was_updated = false;
// };



    // Todo: tolerance must be updated every time pareto_front changes, should track pareto_front updates.
    // Pareto front should be assigned to FunctionUC and parameter wheather it was updated should be also attached to FunctionUC.




void Simplex::update_tolerance_values(vector<Simplex*> simpls, FunctionUC* func) {     // Simpls are all simplices in the partition
    for (int sid=0; sid < simpls.size(); sid++) { // each simplex
        for (int c=0; c < func->_C; c++) {  // each criteria
            if (simpls[sid]->_should_tolerance_be_updated or Simplex::glob_Ls_was_updated[c]) {
                simpls[sid]->_min_lbs = simpls[sid]->find_min_vert_lb_mins(simpls[sid], Simplex::glob_Ls);
            };
            if (func->_pareto_front_was_updated or simpls[sid]->_should_tolerance_be_updated or Simplex::glob_Ls_was_updated[c]) {
                simpls[sid]->_tolerance = simpls[sid]->find_tolerance(func->_pareto_front);;
                simpls[sid]->_should_tolerance_be_updated = false;
            };
        };
    };

    // Reset flags
    for (int c=0; c < func->_C; c++) {
        Simplex::glob_Ls_was_updated[c] = false;
    };
    func->_pareto_front_was_updated = false;
};


double Determinant(double **a, int n) {
    // Based on http://paulbourke.net/miscellaneous/determinant/
    int i, j, j1, j2;
    double det = 0;
    double **m = NULL;

    if (n < 1) { /* Error */ cout << "Determinant cannot be calculated for empty matrix" << endl;
    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = (double**) malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = (double*) malloc((n-1)*sizeof(double));
            for (i=1; i<n; i++) {
                j2 = 0;
                for (j=0; j<n; j++) {
                    if (j == j1) continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * Determinant(m,n-1);
            for (i=0;i<n-1;i++) free(m[i]);
            free(m);
        }
    }
    return(det);
};

#endif
