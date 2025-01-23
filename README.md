# MatrixManipulation

to run this 
```bash
cd build && cmake ..
cmake --build .
```
or just run it on clion or vs code 
<img width="460" alt="Screenshot 2025-01-23 at 23 06 31" src="https://github.com/user-attachments/assets/d3e32c36-6462-4505-a0c0-0e0dbaa421a6" />
This code defines the implementation of a Matrix class in C++ that provides a variety of matrix operations and functionalities. Here is a brief overview of the main components and methods:

Constructors:

Matrix(const Mat2D &mat_vals): Initializes a matrix with given 2D vector values.
Matrix(const size_t rows, const size_t cols, const MatrixType type): Initializes a matrix with given dimensions and type (e.g., identity matrix).
Methods:

void print() const: Prints the matrix to the console.
float &at(const size_t i, const size_t j): Accesses an element at (i, j) with bounds checking.
const float &at(const size_t i, const size_t j) const: Accesses an element at (i, j) with bounds checking (const version).
Matrix multiply(const Matrix &other) const: Multiplies this matrix by another matrix.
float determinant() const: Calculates the determinant of the matrix.
std::string toLatex() const: Converts the matrix to a LaTeX formatted string.
bool compileToPDF(const std::string &filename, const bool DarkMode) const: Compiles the matrix to a PDF file using LaTeX.
bool generatePDF(const std::string &filename, const bool DarkMode) const: Generates a PDF file of the matrix.
Matrix getMinor(const size_t row, const size_t col) const: Gets the minor matrix after removing a specified row and column.
float cofactor(const size_t row, const size_t col) const: Calculates the cofactor of an element.
Matrix getCofactorMatrix() const: Gets the cofactor matrix.
Matrix getAdjugateMatrix() const: Gets the adjugate matrix.
Matrix inverse() const: Calculates the inverse of the matrix.
std::string matrixToLatex() const: Converts the matrix to a LaTeX formatted string.
float trace() const: Calculates the trace of the matrix.
Matrix elementWiseMultiply(const Matrix &other) const: Performs element-wise multiplication (Hadamard product).
std::vector rowMeans() const: Calculates the mean of each row.
Matrix transpose() const: Transposes the matrix.
float vectorNorm(const std::vector &vec): Calculates the norm of a vector.
bool isUpperTriangular(const float tolerance) const: Checks if the matrix is upper triangular.
std::vector getDiagonal() const: Gets the diagonal elements of the matrix.
std::pair<Matrix, Matrix> qrDecomposition() const: Performs QR decomposition.
std::vector calculateEigenvalues(int maxIterations) const: Calculates the eigenvalues.
Matrix pow(const int n) const: Raises the matrix to the power of n.
Matrix kroneckerProduct(const Matrix &other) const: Calculates the Kronecker product.
std::vector colMeans() const: Calculates the mean of each column.
float sum() const: Calculates the sum of all elements.
float mean() const: Calculates the mean of all elements.
std::pair<Matrix, Matrix> luDecomposition() const: Performs LU decomposition.
std::pair<Matrix, Matrix> choleskyDecomposition() const: Performs Cholesky decomposition.
std::string generateHistogram(const std::vector &values): Generates a histogram using LaTeX/TikZ.
std::string generateHeatmap() const: Generates a heatmap using LaTeX/TikZ.
float frobeniusNorm() const: Calculates the Frobenius norm.
float maxNorm() const: Calculates the maximum norm.
bool isSymmetric(const float tolerance) const: Checks if the matrix is symmetric.
bool isOrthogonal(const float tolerance) const: Checks if the matrix is orthogonal.
bool isPositiveDefinite() const: Checks if the matrix is positive definite.
std::pair<size_t, size_t> rank() const: Calculates the rank of the matrix.
float variance() const: Calculates the variance.
void generateLatexReport(const std::string &filename, const bool DarkMode, const bool includeAdvancedStats, const bool includeDecompositions, const bool includeVisualization) const: Generates a comprehensive LaTeX report of the matrix including various properties and visualizations.
This class provides a comprehensive suite of matrix operations, statistical analysis, and LaTeX report generation, making it a versatile tool for matrix manipulation and analysis.
