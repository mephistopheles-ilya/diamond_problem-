#include <vector>
#include <cmath>
#include <random>

#include <iostream>
#include <cmath>

using namespace std;

void findMinEigenvector(double a[9], double eigenvector[3]) {
    const double eps = 1e-12;

    double trace_A = a[0] + a[4] + a[8];
    
    double trace_A2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] 
                    + a[3]*a[3] + a[4]*a[4] + a[5]*a[5]
                    + a[6]*a[6] + a[7]*a[7] + a[8]*a[8];
    
    double coefficient_b = (trace_A * trace_A - trace_A2) / 2.0;
    
    double det_A = a[0]*(a[4]*a[8] - a[5]*a[5]) 
                 - a[1]*(a[1]*a[8] - a[5]*a[2]) 
                 + a[2]*(a[1]*a[5] - a[4]*a[2]);
    
    double a_coeff = -trace_A;
    double b_coeff = coefficient_b;
    double c_coeff = -det_A;
    
    double p = b_coeff - (a_coeff * a_coeff) / 3.0;
    double q = (2.0 * a_coeff * a_coeff * a_coeff) / 27.0 
               - (a_coeff * b_coeff) / 3.0 + c_coeff;
    
    double lambda[3];
    if (fabs(p) < eps) {
        double lambda_single = -a_coeff / 3.0;
        lambda[0] = lambda[1] = lambda[2] = lambda_single;
    } else {
        double denominator = sqrt(-p * p * p / 27.0);
        double theta = acos(fmin(fmax(-q / (2.0 * denominator), -1.0), 1.0)) / 3.0;
        double sqrt_neg_p_over3 = sqrt(-p / 3.0);
        lambda[0] = 2.0 * sqrt_neg_p_over3 * cos(theta) - a_coeff / 3.0;
        lambda[1] = 2.0 * sqrt_neg_p_over3 * cos(theta + 2.0 * M_PI / 3.0) - a_coeff / 3.0;
        lambda[2] = 2.0 * sqrt_neg_p_over3 * cos(theta + 4.0 * M_PI / 3.0) - a_coeff / 3.0;
    }
    
    double min_lambda = lambda[0];
    if (lambda[1] < min_lambda) min_lambda = lambda[1];
    if (lambda[2] < min_lambda) min_lambda = lambda[2];
    
    double m[3][3] = {
        {a[0] - min_lambda, a[1], a[2]},
        {a[1], a[4] - min_lambda, a[5]},
        {a[2], a[5], a[8] - min_lambda}
    };
    
    double v[3];
    v[0] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
    v[1] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
    v[2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];
    
    double norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    
    if (norm < eps) {
        v[0] = m[0][1]*m[2][2] - m[0][2]*m[2][1];
        v[1] = m[0][2]*m[2][0] - m[0][0]*m[2][2];
        v[2] = m[0][0]*m[2][1] - m[0][1]*m[2][0];
        norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    
    if (norm < eps) {
        v[0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
        v[1] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
        v[2] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
        norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    
    if (norm < eps) {
        if (fabs(m[0][0]) < eps && fabs(m[1][1]) < eps && fabs(m[2][2]) < eps) {
            eigenvector[0] = 1.0;
            eigenvector[1] = eigenvector[2] = 0.0;
        } else {
            if (fabs(a[0] - min_lambda) < eps) {
                eigenvector[0] = 1.0;
                eigenvector[1] = eigenvector[2] = 0.0;
            } else if (fabs(a[4] - min_lambda) < eps) {
                eigenvector[1] = 1.0;
                eigenvector[0] = eigenvector[2] = 0.0;
            } else {
                eigenvector[2] = 1.0;
                eigenvector[0] = eigenvector[1] = 0.0;
            }
        }
    } else {
        eigenvector[0] = v[0] / norm;
        eigenvector[1] = v[1] / norm;
        eigenvector[2] = v[2] / norm;
    }
}

int main() {
    double a[9] = {
        4.0, 2.0, 1.0,
        2.0, 5.0, 3.0,
        1.0, 3.0, 6.0
    };
    
    double eigenvector[3];
    findMinEigenvector(a, eigenvector);
    
    cout << "Eigenvector: [";
    cout << eigenvector[0] << ", " << eigenvector[1] << ", " << eigenvector[2] << "]" << endl;
    
    return 0;
}

struct Vector3d {
    double x, y, z;
    Vector3d() : x(0), y(0), z(0) {}
    Vector3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

struct Matrix3x3 {
    double m[3][3]; // m[row][col]
};

Vector3d compute_centroid(const std::vector<Vector3d>& points) {
    Vector3d centroid;
    for (const auto& p : points) {
        centroid.x += p.x;
        centroid.y += p.y;
        centroid.z += p.z;
    }
    size_t n = points.size();
    if (n == 0) return centroid;
    centroid.x /= n;
    centroid.y /= n;
    centroid.z /= n;
    return centroid;
}

Matrix3x3 compute_covariance(const std::vector<Vector3d>& points, const Vector3d& centroid) {
    Matrix3x3 cov = {};
    size_t n = points.size();
    if (n == 0) return cov;

    for (const auto& p : points) {
        double dx = p.x - centroid.x;
        double dy = p.y - centroid.y;
        double dz = p.z - centroid.z;

        cov.m[0][0] += dx * dx;
        cov.m[0][1] += dx * dy;
        cov.m[0][2] += dx * dz;

        cov.m[1][0] += dy * dx;
        cov.m[1][1] += dy * dy;
        cov.m[1][2] += dy * dz;

        cov.m[2][0] += dz * dx;
        cov.m[2][1] += dz * dy;
        cov.m[2][2] += dz * dz;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cov.m[i][j] /= n;
        }
    }

    return cov;
}

Vector3d matrix_multiply(const Matrix3x3& m, const Vector3d& v) {
    Vector3d result;
    result.x = m.m[0][0] * v.x + m.m[0][1] * v.y + m.m[0][2] * v.z;
    result.y = m.m[1][0] * v.x + m.m[1][1] * v.y + m.m[1][2] * v.z;
    result.z = m.m[2][0] * v.x + m.m[2][1] * v.y + m.m[2][2] * v.z;
    return result;
}

double dot_product(const Vector3d& a, const Vector3d& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3d normalize(const Vector3d& v) {
    double len = sqrt(dot_product(v, v));
    if (len == 0.0) return v;
    return Vector3d(v.x / len, v.y / len, v.z / len);
}

Vector3d subtract_component(const Vector3d& v, const Vector3d& direction) {
    double scalar = dot_product(v, direction);
    return Vector3d(v.x - scalar * direction.x,
                    v.y - scalar * direction.y,
                    v.z - scalar * direction.z);
}

Vector3d cross_product(const Vector3d& a, const Vector3d& b) {
    return Vector3d(a.y * b.z - a.z * b.y,
                    a.z * b.x - a.x * b.z,
                    a.x * b.y - a.y * b.x);
}

Vector3d power_method(const Matrix3x3& matrix, const Vector3d& initial_vector,
                      int iterations, const Vector3d* orthogonal_to) {
    Vector3d v = normalize(initial_vector);
    for (int i = 0; i < iterations; ++i) {
        v = matrix_multiply(matrix, v);
        if (orthogonal_to) {
            v = subtract_component(v, *orthogonal_to);
        }
        v = normalize(v);
    }
    return v;
}

Vector3d random_unit_vector() {
    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<double> dist(-1.0, 1.0);

    Vector3d v(dist(rng), dist(rng), dist(rng));
    return normalize(v);
}

void least_squares_plane(const std::vector<Vector3d>& points,
                         Vector3d& centroid, Vector3d& normal) {
    centroid = compute_centroid(points);
    Matrix3x3 cov = compute_covariance(points, centroid);

    // Find the first principal component
    Vector3d v1 = power_method(cov, random_unit_vector(), 100, nullptr);

    // Find the second principal component, orthogonal to v1
    Vector3d v2_initial = random_unit_vector();
    v2_initial = subtract_component(v2_initial, v1);
    Vector3d v2 = power_method(cov, normalize(v2_initial), 100, &v1);

    // Compute the normal as the cross product of v1 and v2
    normal = cross_product(v1, v2);
    normal = normalize(normal);

    // Adjust normal direction based on majority of points
    double sum = 0.0;
    for (const auto& p : points) {
        Vector3d vec = {p.x - centroid.x, p.y - centroid.y, p.z - centroid.z};
        sum += dot_product(normal, vec);
    }
    if (sum < 0.0) {
        normal.x = -normal.x;
        normal.y = -normal.y;
        normal.z = -normal.z;
    }
}
