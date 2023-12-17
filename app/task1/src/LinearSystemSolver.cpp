#include "LinearSystemSolver.h"

template<>
double LinearSystemSolver::selectMainElement<FirstNotZeroInCol>(Matrix& System, Eigen::Index i)
{
    const auto notZeroFilter = [](double element) { return element != 0.0; };

    auto column = System.col(i);
    const auto first = column.begin() + i;
    const auto last = column.end();

    const auto elementIter = std::ranges::find_if(first, last, notZeroFilter);
    if (elementIter == last) throw std::logic_error("The matrix has linearly dependent columnds");

    const double mainElement = *elementIter;
    const auto idx = std::distance(column.begin(), elementIter);
    System.swap(i, idx, Matrix::Row);

    return mainElement;
}

template<>
double LinearSystemSolver::selectMainElement<AbsMaxInCol>(Matrix& System, Eigen::Index i)
{
    const auto absComp = [](double lhs, double rhs) { return std::abs(lhs) < std::abs(rhs); };

    auto column = System.col(i);
    const auto first = column.begin() + i;
    const auto last = column.end();

    const auto elementIter = std::ranges::max_element(first, last, absComp);
    const double mainElement = *elementIter;
    if (mainElement == 0.0) throw std::logic_error("The matrix has linearly dependent columnds");

    const auto idx = std::distance(column.begin(), elementIter);
    System.swap(i, idx, Matrix::Row);

    return mainElement;
}