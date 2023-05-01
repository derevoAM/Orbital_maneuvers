#ifndef ORBITAL_MANEUVERS_VECTOR_H
#define ORBITAL_MANEUVERS_VECTOR_H


#include <iostream>
#include <vector>
#include <cmath>


template<typename T>
std::vector<T> operator+(const std::vector<T> &vec_1_, const std::vector<T> &vec_2_) {
    std::vector<T> res_(3);
    for (int i = 0; i < vec_1_.size(); i++) res_[i] = vec_1_[i] + vec_2_[i];
    return res_;
}

template<typename T>
std::vector<T> operator*(const std::vector<T> &vec_1_, T mult_) {
    std::vector<T> res_(3);
    for (int i = 0; i < vec_1_.size(); i++) res_[i] = vec_1_[i] * mult_;
    return res_;
}

template<typename T>
std::vector<T> operator/(const std::vector<T> &vec_1_, T mult_) {
    std::vector<T> res_;
    res_.reserve(vec_1_.size());
    for (int i = 0; i < vec_1_.size(); i++) res_.push_back(vec_1_[i] / mult_);
    return res_;
}

template<typename T>
T scalar(const std::vector<T> &mult_1, const std::vector<T> &mult_2)
{
    T sum = 0;
    for(int i = 0; i < mult_1.size(); i ++) sum += mult_1[i] * mult_2[i];
    return sum;
}

template<typename T>
std::vector<T> cross_product(const std::vector<T> &mult_1, const std::vector<T> &mult_2)
{
    std::vector<T> res_(3);
    res_[0] = mult_1[1] * mult_2[2] - mult_1[2] * mult_2[1];
    res_[1] = mult_1[2] * mult_2[0] - mult_1[0] * mult_2[2];
    res_[2] = mult_1[0] * mult_2[1] - mult_1[1] * mult_2[0];
    return res_;
}

template<typename T>
std::vector<T> operator-(const std::vector<T> &vec_1_, const std::vector<T> &vec_2_) {
    std::vector<T> res_(3);
    for (int i = 0; i < vec_1_.size(); i++) res_[i] = vec_1_[i] - vec_2_[i];
    return res_;
}

template<typename T>
T norm(const std::vector<T> &vec_) {
    T res_ = 0;
    for (int i = 0; i < vec_.size(); i++) res_ += vec_[i] * vec_[i];
    return std::sqrt(res_);
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec_) {

    for (int j = 0; j < vec_.size(); j++) out << vec_[j] << " ";
    return out;
}

#endif //ORBITAL_MANEUVERS_VECTOR_H
