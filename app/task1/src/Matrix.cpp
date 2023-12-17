#include "Matrix.h"
#include "LinearSystemSolver.h"

Matrix& Matrix::inverseInPlace()
{
    Eigen::Index n = rows();
    const Matrix I = Identity(n, n);

    *this = LinearSystemSolver::gauss<AbsMaxInCol>(*this, I);
    return *this;
}

Matrix Matrix::inverse() const
{
    Matrix M = *this;
    return M.inverseInPlace();
}

void Matrix::swap(Eigen::Index i, Eigen::Index j, ColSwapper)
{
    col(i).swap(col(j));
}

void Matrix::swap(Eigen::Index i, Eigen::Index j, RowSwapper)
{
    row(i).swap(row(j));
}