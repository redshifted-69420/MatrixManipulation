#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <numeric>
#include <vector>

using Mat2D = std::vector<std::vector<float>>;
using Mat1D = std::vector<float>;

enum MatrixType { Identity, Zeros };

class Matrix {
private:
  Mat2D data_;
  size_t rows_;
  size_t cols_;
  static const float TOLERANCE;
  [[nodiscard]] static std::string
  generateHistogram(const std::vector<float> &values);
  [[nodiscard]] std::string generateHeatmap() const;

public:
  explicit Matrix(const Mat2D &mat_vals);
  Matrix(size_t rows, size_t cols, MatrixType type = Zeros);
  void print() const;
  [[nodiscard]] size_t getRows() const { return rows_; }
  [[nodiscard]] size_t getCols() const { return cols_; }
  float &at(size_t i, size_t j);
  [[nodiscard]] const float &at(size_t i, size_t j) const;
  [[nodiscard]] Matrix multiply(const Matrix &other) const;
  Matrix operator*(const Matrix &other) const { return multiply(other); }
  [[nodiscard]] float determinant() const;
  [[nodiscard]] std::string toLatex() const;
  [[nodiscard]] bool compileToPDF(const std::string &filename,
                                  bool DarkMode = false) const;
  [[nodiscard]] bool generatePDF(const std::string &filename = "matrix",
                                 bool DarkMode = false) const;
  [[nodiscard]] Matrix getCofactorMatrix() const;
  [[nodiscard]] Matrix getAdjugateMatrix() const;
  [[nodiscard]] Matrix inverse() const;
  [[nodiscard]] Matrix getMinor(size_t row, size_t col) const;
  [[nodiscard]] float cofactor(size_t row, size_t col) const;
  [[nodiscard]] std::string toLatexString() const;
  void generateLatexReport(const std::string &filename, bool DarkMode = false,
                           bool includeAdvancedStats = true,
                           bool includeDecompositions = true,
                           bool includeVisualization = true) const;
  [[nodiscard]] std::string matrixToLatex() const;
  [[nodiscard]] float trace() const;
  [[nodiscard]] std::vector<float>
  calculateEigenvalues(int maxIterations = 100) const;
  [[nodiscard]] std::pair<Matrix, Matrix> qrDecomposition() const;
  [[nodiscard]] bool isUpperTriangular(float tolerance = TOLERANCE) const;
  [[nodiscard]] Matrix transpose() const;
  [[nodiscard]] static float vectorNorm(const std::vector<float> &vec);
  [[nodiscard]] std::vector<float> getDiagonal() const;
  [[nodiscard]] float frobeniusNorm() const;
  [[nodiscard]] float maxNorm() const;
  [[nodiscard]] bool isSymmetric(float tolerance = TOLERANCE) const;
  [[nodiscard]] bool isOrthogonal(float tolerance = TOLERANCE) const;
  [[nodiscard]] bool isPositiveDefinite() const;
  [[nodiscard]] std::pair<size_t, size_t> rank() const;
  [[nodiscard]] Matrix pow(int n) const;
  [[nodiscard]] Matrix elementWiseMultiply(const Matrix &other) const;
  [[nodiscard]] Matrix kroneckerProduct(const Matrix &other) const;
  [[nodiscard]] std::vector<float> rowMeans() const;
  [[nodiscard]] std::vector<float> colMeans() const;
  [[nodiscard]] float sum() const;
  [[nodiscard]] float mean() const;
  [[nodiscard]] float variance() const;
  [[nodiscard]] std::pair<Matrix, Matrix> luDecomposition() const;
  [[nodiscard]] std::pair<Matrix, Matrix> choleskyDecomposition() const;
  [[nodiscard]] std::vector<std::vector<float>> findNullVectors() const;
};

constexpr float Matrix::TOLERANCE = 1e-10;

#endif // MATRIX_H
