#pragma once
#include <complex>
#include <Eigen/Dense>

namespace ecg1d {

using Cd = std::complex<double>;
using MatrixXcd = Eigen::MatrixXcd;
using VectorXcd = Eigen::VectorXcd;
using MatrixXi = Eigen::MatrixXi;
using VectorXi = Eigen::VectorXi;

struct CosineTerm {
    double C = 0.0;
    double A = 0.0;
    double b0 = 0.0;
};

} // namespace ecg1d
