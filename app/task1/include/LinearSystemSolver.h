#ifndef LINEARSYSTEMSOLVER_H
#define LINEARSYSTEMSOLVER_H

#include "Matrix.h"

#include <iostream>
#include <ranges>

template <typename T>
auto in_range(T first, T last)
{
    assert(first <= last);
    return std::views::iota(first, last);
}

template <typename T>
auto in_range_reverse(T first, T last)
{
    assert(first >= last);
    std::swap(first, last);
    return std::views::reverse(in_range(first, last));
}

enum SelectingPolicy
{
    FirstNotZeroInCol = 0,
    AbsMaxInCol = 1
};

class LinearSystemSolver
{
public:

    template <SelectingPolicy Policy = FirstNotZeroInCol>
    static Matrix gauss(const Matrix& A, const Matrix& B)
    {
        Matrix System = Matrix::Zero(A.rows(), A.cols() + B.cols());
        System << A, B;

        const Eigen::Index n = System.rows();
        Matrix X = Matrix::Zero(n, B.cols());

        directMotion<Policy>(System);
        reverseMotion<Policy>(System, X);

        return X;
    }

private:

    template <SelectingPolicy Policy>
    static double selectMainElement(Matrix& System, Eigen::Index i)
    {
        return 0.0;
    };

    template <SelectingPolicy Policy>
    static void restoreOrder(){};

    template <SelectingPolicy Policy>
    static void directMotion(Matrix& System)
    {
        const Eigen::Index n = System.rows();
        for (auto i : in_range(0L, n))
        {
            const double mainElement = selectMainElement<Policy>(System, i);
            auto row = System.row(i);
            row /= mainElement;
            for (auto j : in_range(i + 1, n))
            {
                System.row(j) -= row * System(j, i);
            }
        }
    }

    template <SelectingPolicy Policy>
    static void reverseMotion(const Matrix& System, Matrix& X)
    {
        const Eigen::Index n = System.rows();

        for (auto i : in_range(0L, X.cols()))
        {
            X(n - 1, i) = System(n - 1, n + i);
            for (auto j : in_range_reverse(n - 1, 0L))
            {
                auto aSlice = Eigen::seq(j + 1, n - 1);
                auto xSlice = Eigen::seq(j + 1, n - 1);

                auto a = System(j, aSlice);
                auto x = X(xSlice, i);

                X(j, i) = System(j, n + i) - a.dot(x);
            }
        }

        restoreOrder<Policy>();
    }

public:

    static Matrix lu(const Matrix& A, const Matrix& B)
    {
        const auto n = A.rows();
        const auto m = B.cols();

        Matrix L;
        Matrix U;
        LUDecomposition(A, L, U);

        Eigen::MatrixXd X(n, m);

        for (int j = 0; j < m; ++j) {
            Eigen::VectorXd b = B.col(j);

            Eigen::VectorXd y(n);
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int k = 0; k < i; ++k) {
                    sum += L(i, k) * y(k);
                }
                y(i) = (b(i) - sum) / L(i, i);
            }

            Eigen::VectorXd x(n);
            for (int i = n - 1; i >= 0; --i) {
                double sum = 0.0;
                for (int k = i + 1; k < n; ++k) {
                    sum += U(i, k) * x(k);
                }
                x(i) = (y(i) - sum) / U(i, i);
            }

            X.col(j) = x;
        }

        return X;
    };

private:

    static void LUDecomposition(const Matrix& A, Matrix& L, Matrix& U)
    {
        const auto n = A.rows();

        L = Matrix::Identity(n, n);
        U = A;

        for (int k = 0; k < n - 1; ++k) {
            for (int i = k + 1; i < n; ++i) {
                if (U(k, k) != 0.0) {
                    L(i, k) = U(i, k) / U(k, k);
                    U.row(i) -= L(i, k) * U.row(k);
                }
            }
        }

        std::cout << "\n\n" << L << "\n\n" << U << "\n\n\n";
    }

};

#endif // LINEARSYSTEMSOLVER_H
