#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
typedef Eigen::Matrix<double, 2, Eigen::Dynamic> MatrixX2;
typedef Eigen::Matrix<double, 3, Eigen::Dynamic> MatrixX3;
typedef Eigen::Index Index;
typedef Eigen::SparseMatrix<double, 0, int> SparseMatrixX;
template <typename T> using VectorS = Eigen::Matrix<T, Eigen::Dynamic, 1>;
typedef VectorS<double> VectorX;
template <typename T> using FunctorAD = std::function<T(const VectorS<T> &)>;
