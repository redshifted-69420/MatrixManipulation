```c++
    Mat2D values = {
        {2.0f, 2.0f, 3.0f},
        {4.0f, 5.0f, 6.0f},
        {7.0f, 8.0f, 9.0f}
    };
    Matrix mat1(values);
    mat1.print();

    Matrix mat2(3, 3, Identity);
    mat2.print();

    Matrix mat3(4, 4, Zeros);
    mat3.print();

    Mat2D values1 = {
        {1.0f, 2.0f},
        {3.0f, 4.0f}
    };
    Mat2D values2 = {
        {5.0f, 6.0f},
        {7.0f, 8.0f}
    };

    Matrix mat5(values1);
    Matrix mat4(values2);

    Matrix result = mat5 * mat4;

    result.print();

    Mat2D values6 = {
        {4.0f, 3.0f, 2.0f},
        {2.0f, 1.0f, 5.0f},
        {2.0f, 1.0f, 9.0f}
    };

    Matrix mat(values6);
    std::cout << "Matrix:" << std::endl;
    mat.print();

    float det = mat.determinant();
    std::cout << "Determinant: " << det << std::endl;

    Mat2D matrix = {
        {-1.0, 2.0, 3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0, 9.0, 10.0},
        {11.0, 12.1, 13.4, 14.4, 15.8},
        {16.9, 17.8, 18.0, 19.7, 20.6},
        {21.3, 22.2, 23.6, 24.1, 25.0}
    };

    Matrix mat9(matrix);
    std::cout << "Matrix:" << std::endl;
    mat9.print();

    float det2 = mat9.determinant();
    std::cout << "Determinant: " << det2 << std::endl;

    std::cout << "LaTeX representation:\n"
              << mat9.toLatex() << std::endl;

    if (mat9.generatePDF("matrix", true)) {
        std::cout << "PDF generated successfully as 'matrix.pdf'!\n\n";
    } else {
        std::cout << "Error generating PDF\n";
    }

    Matrix mat10(values6);
    std::cout << "Original Matrix:" << std::endl;
    mat10.print();

    try {
        Matrix inv = mat10.inverse();
        std::cout << "Inverse Matrix:" << std::endl;
        inv.print();

        Matrix identity = mat10 * inv;
        std::cout << "Verification (should be identity matrix):" << std::endl;
        identity.print();
    } catch (const std::runtime_error &e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    std::cout << "LaTeX report generated. Compile with:\n"
              << "pdflatex matrix_report.tex\n";

    std::cout << "Trace of Matrix is: " << mat10.trace() << std::endl;

    try
    {
        std::vector<float> eigenvalues = mat10.calculateEigenvalues();

        std::cout << "Eigenvalues:" << std::endl;
        for (float eigenvalue : eigenvalues)
        {
            std::cout << eigenvalue << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cout << "Error: " << e.what() << std::endl;
    }
```
