#include <iostream>

#include "LinearSystemSolver.h"

int main()
{
    const Matrix I = Matrix::Identity(3, 3);
    Matrix A(3, 3);
    A << 5, 2, 3, 4, 5, 6, 7, 8, 9;
    std::cout << A << "\n";

    try
    {
        //std::cout << A.inverse();
        std:: cout << LinearSystemSolver::lu(A, I);
    }
    catch (const std::logic_error& e)
    {
        std::cout << e.what() << "\nExit...";
    }

    return 0;
}