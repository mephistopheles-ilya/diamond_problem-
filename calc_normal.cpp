#include <vector>
#include <cmath>
#include <random>

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
