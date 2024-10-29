
template <int N>
class point {
    double coords[N];
public:
    point(const std::array<double, N>& arr) {
        std::copy(arr.begin(), arr.end(), coords);
    }

    double norm2() const {
        double sum = 0;
        for (int i = 0; i < N; i++)
            sum += coords[i] * coords[i];
        return sqrt(sum);
    }
    double operator[](int i) const {
        return coords[i];
    }
};