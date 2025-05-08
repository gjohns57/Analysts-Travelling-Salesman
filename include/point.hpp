#pragma once
#include <math.h>
template <int N>
class point {
    double coords[N];
public:


    point() {
        for (int i = 0; i < N; i++)
            coords[i] = 0;
    }

    double norm2() const {
        double sum = 0;
        for (int i = 0; i < N; i++)
            sum += coords[i] * coords[i];
        return sqrt(sum);
    }

    double& operator[](int i) {
        return coords[i];
    }

    point<N> operator+(const point<N>& other) const {
        point<N> result;
        for (int i = 0; i < N; i++)
            result.coords[i] = this->coords[i] + other.coords[i];
        return result;
    }

    point<N> operator-(const point<N>& other) const {
        point<N> result;
        for (int i = 0; i < N; i++)
            result.coords[i] = this->coords[i] - other.coords[i];
        return result;
    }

    double operator*(const point<N>& other) const {
        double result = 0;
        for (int i = 0; i < N; i++)
            result += this->coords[i] * other.coords[i];
        return result;
    }

    point<N> operator*(double scalar) const {
        point<N> result;
        for (int i = 0; i < N; i++)
            result.coords[i] = this->coords[i] * scalar;
        return result;
    }

    point<N> operator/(double scalar) const {
        point<N> result;
        for (int i = 0; i < N; i++)
            result.coords[i] = this->coords[i] / scalar;
        return result;
    }
};
