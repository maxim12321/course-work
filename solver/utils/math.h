#pragma once

template <typename T>
T pow(const T &value, int n) {
    if (n == 0) {
        return value / value;
    }
    if (n == 1) {
        return value;
    }

    auto part = pow(value, n / 2);
    part *= part;

    return n % 2 == 0 ? part : part * value;
}

const long double kEpsilon = 0.0000000001;

inline bool is_zero(const double &value) {
    return std::abs(value) < kEpsilon;
}   

inline bool is_zero(const long double &value) {
    return std::abs(value) < kEpsilon;
}
