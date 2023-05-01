#ifndef ORBITAL_MANEUVERS_DENSE_H
#define ORBITAL_MANEUVERS_DENSE_H


#include <vector>
#include <iostream>
#include <cmath>

template<typename T>
class Dense {
private:
    std::vector<T> matrix_;
    int h_, w_;
public:
    Dense(int h_, int w_, const std::vector<T> input_) : h_(h_), w_(w_), matrix_(input_) {};

    Dense(int h_, int w_) : h_(h_), w_(w_) {
        matrix_.resize(w_ * h_);
    };

    T operator()(int i, int j) const { return matrix_[i * w_ + j]; }

    T &operator()(int i, int j) { return matrix_[i * w_ + j]; }

    int get_height() const { return h_; }

    int get_width() const { return w_; }

    Dense<T> operator*(const Dense<T> &mult_) {
        std::vector<T> result_;
        result_.reserve(h_ * mult_.w_);
        for (int i = 0; i < h_; i++) {
            for (int k = 0; k < mult_.w_; k++) {
                T sum = static_cast<T>(0);
                for (int j = 0; j < w_; j++) sum += matrix_[i * w_ + j] * mult_(j, k);
                result_.push_back(sum);
            }
        }
        return Dense<T>{h_, mult_.w_, result_};
    }

    std::vector<T>  operator*(const std::vector<T> &mult_) const {
        std::vector<T> result_;
        result_.reserve(mult_.size());
        for (int i = 0; i < h_; i++) {
            T sum = static_cast<T>(0);
            for (int k = 0; k < mult_.size(); k++) {

                sum += matrix_[i * w_ + k] * mult_[k];

            }
            result_.push_back(sum);
        }
        return result_;
    }

    Dense<T> operator*(T mult_) {
        std::vector<T> result_;
        result_.reserve(h_ * w_);
        for (int i = 0; i < h_; i++) for (int j = 0; j < w_; j++) result_.push_back(matrix_[i * w_ + j] * mult_);
        return Dense<T>{h_, w_, result_};
    }

    Dense<T> operator+(const Dense<T> &add_) {
        std::vector<T> result_;
        result_.reserve(w_ * h_);
        for (int i = 0; i < h_; i++) {
            for (int j = 0; j < w_; j++) result_.push_back(matrix_[i * w_ + j] + add_(i, j));
        }
        return Dense<T>{h_, w_, result_};
    }

    Dense<T> operator+=(const Dense<T> &add_) {
        for (int i = 0; i < h_; i++) {
            for (int j = 0; j < w_; j++) matrix_[i * w_ + j] += add_(i, j);
        }
        return *this;
    }

    Dense<T> Transpose() {
        std::vector<T> trans_;
        trans_.reserve(h_ * w_);
        for (int i = 0; i < w_; i++) {
            for (int j = 0; j < h_; j++) {
                trans_.push_back(matrix_[j * w_ + i]);
            }
        }
        return Dense<T>{w_, h_, trans_};
    }

    Dense<T> get_column(int j) const { // getting the j column from matrix
        std::vector<T> col_;
        col_.reserve(h_);
        for (int i = 0; i < h_; i++) col_.push_back(matrix_[i * w_ + j]);
        return Dense<T>{h_, 1, col_};
    }

    T norm() const // only for vectors
    {
        T sum = static_cast<T>(0);
        for (int i = 0; i < h_; i++) sum += matrix_[i] * matrix_[i];
        return sqrt(static_cast<double>(sum));
    }

    T scalar(const Dense<T> &mult_) // only for vectors
    {
        T sum = static_cast<T>(0);
        for (int i = 0; i < h_; i++) sum += matrix_[i] * mult_(i, 0);
        return sum;
    }

    Dense<T> eval_norm_vec(
            int i) { // evaluation of the normal vector, taking into account that the first i coordinates are equal to 0
        Dense<T> res_{h_, 1, matrix_};
        for (int j = 0; j < i; j++) res_(j, 0) = 0;
        if (res_(i, 0) >= 0) {
            res_(i, 0) += res_.norm();
        } else res_(i, 0) -= res_.norm();
        double cur_norm = res_.norm();
        for (int j = i; j < h_; j++) res_(j, 0) /= cur_norm;
        return res_;
    }

    friend std::ostream &operator<<(std::ostream &out, const Dense<T> &output_) {
        for (int i = 0; i < output_.h_; i++) {
            for (int j = 0; j < output_.w_; j++) out << output_(i, j) << " ";
            out << "\n";
        }
        return out;
    }
};


#endif //ORBITAL_MANEUVERS_DENSE_H
