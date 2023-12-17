#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>

class Matrix : public Eigen::MatrixXd
{
    using Eigen::MatrixXd::MatrixXd;

    struct ColSwapper{};
    struct RowSwapper{};

public:

    static ColSwapper Column;
    static RowSwapper Row;

    void swap(Eigen::Index i, Eigen::Index j, ColSwapper);
    void swap(Eigen::Index i, Eigen::Index j, RowSwapper);

    [[nodiscard]] Matrix inverse() const;
    Matrix& inverseInPlace();
};

#endif // MATRIX_H
