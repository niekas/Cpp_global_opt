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
        _le_v1 = 0;
        _le_v2 = 0;
        _diameter = 0;
        _min_lb = numeric_limits<double>::max();

        _is_in_partition = true;
        _should_be_divided = false;
        _should_lb_mins_be_updated = true;
    };

    int _D;             // Dimension of variable space

    Point* _le_v1;      // First vertex of the longest edge
    Point* _le_v2;      // Second vertex of the longest edge
    double _diameter;   // Longest edge length

    double _min_L;      // Minimum L for this simplex
    double _min_lb;     // Minimum of the lower bound over simplex found using vertex with the lowest value

    static double alpha;              // Coeficient of search globality
    static double glob_L;             // Globally known biggest min L
    static bool glob_L_was_updated;

    vector<Point*> _verts;    // Vertices of this simplex (points with coordinates and values)
    bool _is_in_partition;    // Is this simplex in the current partition
    bool _should_be_divided;  // Should this simplex be divided in next iteration
    bool _should_lb_mins_be_updated;   // Should the minimums of Lipschitz lower bound be updated

    Point* _min_vert;   // Vertex with lowest function value
    Point* _max_vert;   // Vertex with lowest function value

    void print(){
        cout << " Simplex   diam:  " << _diameter << "   tol:  " << _min_lb;
        cout << "   min_L:  " << _min_L <<  endl;
        for (int i=0; i < _verts.size(); i++){
            _verts[i]->print();
        };
    };

    static void log_partition(vector<Simplex*> simplexes,
                              vector<Simplex*> selected,
                              Function* funcs,
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
               log_file << " (" << simplexes[i]->_verts[j]->_value << ");";
           };
           log_file << " ("<< simplexes[i]->_diameter << "," << simplexes[i]->_min_lb << ")" << endl;
       };


       log_file << "Selected:" << endl;
       for (int i=0; i < selected.size(); i++) {
           for (int j=0; j < selected[i]->_verts.size(); j++) {
               for (int k=0; k < selected[i]->_verts[j]->size(); k++){
                    log_file << selected[i]->_verts[j]->_X[k] << " ";
               };
               log_file << " (" << selected[i]->_verts[j]->_value<<");";
           };
           log_file << " ("<< selected[i]->_diameter << "," << selected[i]->_min_lb << ")" << endl;
       };


       log_file.close();
    };


    void init_parameters(Function* func) {   // Called when all verts have been added
        _D = _verts.size() - 1;

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
        _min_L = find_simplex_min_L();

        // Initialize or update global L if needed
        if (Simplex::glob_L == numeric_limits<double>::max()) {
            Simplex::glob_L = Simplex::alpha * _min_L;
        } else {
            if (Simplex::glob_L < Simplex::alpha * _min_L) {
                Simplex::glob_L = Simplex::alpha * _min_L;
                Simplex::glob_L_was_updated = true;
            };
        };

        // Find vertex with minimum function value
        _min_vert = _verts[0];
        _max_vert = _verts[_verts.size()-1];
        for (int i=0; i < _verts.size(); i++) {
            if (_verts[i]->_value < _min_vert->_value) {
                _min_vert = _verts[i];
            };
        };
    };

    double find_min_vert_lb_min(Simplex* simpl, double L) {
        // Finds minimum of lower bound, which is constructed from the vertex with lowest function value
        // return _min_vert->_value - L * simpl->_diameter;
        return _max_vert->_value - L * simpl->_diameter;
    };

    double find_simplex_min_L() {
        // Finds minimum L for this simplex by finding min L for each simplex edge
        double dist;
        double f_diff;
        double edge_L;
        double max_edge_L = -numeric_limits<double>::max();
        for (int i=0; i < _verts.size(); i++) {
            for (int j=i+1; j < _verts.size(); j++) {
                f_diff = fabs(_verts[i]->_value - _verts[j]->_value);
                dist = l2norm(_verts[i], _verts[j]);
                edge_L = f_diff / dist;     // Note: maybe dist (division by zero) protection is needed?
                if (edge_L > max_edge_L) {  //       Practically this case does not occur.
                    max_edge_L = edge_L;
                };
            };
        };
        return max_edge_L;
    };

    void add_vertex(Point* vertex){
        _verts.push_back(vertex);
        vertex->_simplices.push_back(this);
    };

    int size() {
        return _verts.size();
    };

    static void update_min_lb_values(vector<Simplex*> simpls, Function* func);

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
        for (int i=0; i < _verts.size(); i++) {
            _verts[i]->_simplices.erase(remove(_verts[i]->_simplices.begin(), _verts[i]->_simplices.end(), this), _verts[i]->_simplices.end());
        };
        _verts.clear();
    };
};
bool Simplex::glob_L_was_updated = false;
double Simplex::glob_L = numeric_limits<double>::max();
double Simplex::alpha = numeric_limits<double>::max();


void Simplex::update_min_lb_values(vector<Simplex*> simpls, Function* func) {
    for (int sid=0; sid < simpls.size(); sid++) {
        if (simpls[sid]->_should_lb_mins_be_updated or Simplex::glob_L_was_updated) {
            simpls[sid]->_min_lb = simpls[sid]->find_min_vert_lb_min(simpls[sid], Simplex::glob_L);
            simpls[sid]->_should_lb_mins_be_updated = false;
        };
    };
    Simplex::glob_L_was_updated = false;
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
