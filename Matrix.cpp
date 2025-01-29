#include "Matrix.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

Matrix::Matrix(const Mat2D &mat_vals) {
  if (mat_vals.empty() || mat_vals[0].empty()) {
    throw std::invalid_argument("Matrix cannot be empty");
  }
  rows_ = mat_vals.size();
  cols_ = mat_vals[0].size();
  for (const auto &row : mat_vals) {
    if (row.size() != cols_) {
      throw std::invalid_argument("All rows must have the same length");
    }
  }
  data_ = mat_vals;
}

Matrix::Matrix(const size_t rows, const size_t cols, const MatrixType type)
    : rows_(rows), cols_(cols) {
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("Matrix dimensions must be positive");
  }

  if (type == Identity && rows != cols) {
    throw std::invalid_argument("Identity matrix must be square");
  }

  data_ = Mat2D(rows, std::vector<float>(cols, 0.0L));

  if (type == Identity) {
    for (size_t i = 0; i < rows; ++i) {
      data_[i][i] = 1.0L;
    }
  }
}

void Matrix::print() const {
  for (size_t i = 0; i < rows_; ++i) {
    std::cout << " [";
    for (size_t j = 0; j < cols_; ++j) {
      std::cout << "\033[32m" << data_[i][j] << "\033[0m" << " ";
    }
    std::cout << "\b]";
    std::cout << '\n';
  }
  std::cout << '\n';
}

float &Matrix::at(const size_t i, const size_t j) {
  if (i >= rows_ || j >= cols_) {
    throw std::out_of_range("Matrix indices out of range");
  }
  return data_[i][j];
}

const float &Matrix::at(const size_t i, const size_t j) const {
  if (i >= rows_ || j >= cols_) {
    throw std::out_of_range("Matrix indices out of range");
  }
  return data_[i][j];
}

Matrix Matrix::multiply(const Matrix &other) const {
  if (cols_ != other.getRows()) {
    throw std::invalid_argument(
        "Matrix dimensions incompatible for multiplication");
  }

  Mat2D result(rows_, std::vector(other.getCols(), 0.0f));

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < other.getCols(); ++j) {
      for (size_t k = 0; k < cols_; ++k) {
        result[i][j] += data_[i][k] * other.at(k, j);
      }
    }
  }
  return Matrix(result);
}

float Matrix::determinant() const {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Determinant calculation is only implemented for square matrices");
  }

  if (rows_ == 2) {
    return (data_[0][0] * data_[1][1]) - (data_[0][1] * data_[1][0]);
  }

  float det = 0.0L;

  for (size_t col = 0; col < cols_; col++) {
    Mat2D minor;
    for (size_t i = 1; i < rows_; ++i) {
      std::vector<float> row;
      for (size_t j = 0; j < cols_; ++j) {
        if (j != col) {
          row.push_back(data_[i][j]);
        }
      }
      minor.push_back(row);
    }
    Matrix minorMatrix(minor);
    const float sign = (col % 2 == 0) ? 1.0 : -1.0;
    det += sign * data_[0][col] * minorMatrix.determinant();
  }
  return det;
}

std::string Matrix::toLatex() const {
  std::stringstream latex;
  latex << "\\begin{bmatrix}\n";
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      latex << std::fixed << std::setprecision(2) << data_[i][j];
      if (j < cols_ - 1) {
        latex << " & ";
      }
    }
    if (i < rows_ - 1) {
      latex << " \\\\\n";
    } else {
      latex << "\n";
    }
  }
  latex << "\\end{bmatrix}\n";
  return latex.str();
}

bool Matrix::compileToPDF(const std::string &filename,
                          const bool DarkMode) const {
  try {

    generateLatexReport(filename + ".tex", DarkMode);
    const std::string command = "pdflatex -interaction=nonstopmode " +
                                filename + ".tex > /dev/null 2>&1";
    const int result = std::system(command.c_str());

    if (result == 0) {
      const std::string cleanupCommand = "rm -f " + filename + ".aux " +
                                         filename + ".log " + filename +
                                         ".fdb_latexmk " + filename + ".fls";
      std::system(cleanupCommand.c_str());
      return true;
    }

    std::cerr << "pdflatex failed with exit code: " << result << std::endl;
    return false;
  } catch (const std::exception &e) {
    std::cerr << "Error in compileToPDF: " << e.what() << std::endl;
    return false;
  }
}

bool Matrix::generatePDF(const std::string &filename,
                         const bool DarkMode) const {
  try {
    if (compileToPDF(filename, DarkMode)) {
      return true;
    }
    std::cerr << "Error: Failed to compile PDF for file: " << filename
              << std::endl;
    return false;
  } catch (const std::exception &e) {
    std::cerr << "Error generating PDF: " << e.what() << std::endl;
    return false;
  }
}

Matrix Matrix::getMinor(const size_t row, const size_t col) const {
  Mat2D minorData;

  for (size_t i = 0; i < rows_; ++i) {
    if (i != row) {
      std::vector<float> minorRow;
      for (size_t j = 0; j < cols_; ++j) {
        if (j != col) {
          minorRow.push_back(data_[i][j]);
        }
      }
      minorData.push_back(minorRow);
    }
  }
  return Matrix(minorData);
}

float Matrix::cofactor(const size_t row, const size_t col) const {
  const Matrix minor = getMinor(row, col);
  const float sign = ((row + col) % 2 == 0) ? 1.0L : -1.0L;
  return sign * minor.determinant();
}

Matrix Matrix::getCofactorMatrix() const {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Matrix must be square to compute cofactor matrix");
  }

  Mat2D cofactorData(rows_, std::vector<float>(cols_));

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      cofactorData[i][j] = cofactor(i, j);
    }
  }
  return Matrix(cofactorData);
}

Matrix Matrix::getAdjugateMatrix() const {
  const Matrix cofactorMatrix = getCofactorMatrix();
  Mat2D adjugateData(cols_, std::vector<float>(rows_));

  // Transpose of cofactor matrix
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      adjugateData[j][i] = cofactorMatrix.data_[i][j];
    }
  }
  return Matrix(adjugateData);
}

Matrix Matrix::inverse() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("Only square matrices can be inverted");
  }

  const float det = determinant();
  if (std::abs(det) < TOLERANCE) {
    throw std::runtime_error("Matrix is singular (non-invertible)");
  }

  const Matrix adjugate = getAdjugateMatrix();
  Mat2D inverseData(rows_, std::vector<float>(cols_));

  // Multiply adjugate by 1/determinant
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      inverseData[i][j] = adjugate.data_[i][j] / det;
    }
  }
  return Matrix(inverseData);
}

std::string Matrix::matrixToLatex() const {
  std::stringstream ss;
  ss << "\\begin{bmatrix}\n";
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      ss << std::fixed << std::setprecision(4) << data_[i][j];
      if (j < cols_ - 1)
        ss << " & ";
    }
    if (i < rows_ - 1)
      ss << " \\\\\n";
  }
  ss << "\n\\end{bmatrix}\n";
  return ss.str();
}

float Matrix::trace() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("Trace is only defined for square matrices");
  }

  float sum = 0.0f;
  for (size_t i = 0; i < rows_; ++i) {
    sum += data_[i][i];
  }
  return sum;
}

Matrix Matrix::elementWiseMultiply(const Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument(
        "Matrices must have same dimensions for Hadamard product");
  }

  std::vector result(rows_, std::vector<float>(cols_));
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      result[i][j] = data_[i][j] * other.data_[i][j];
    }
  }
  return Matrix(result);
}

std::vector<float> Matrix::rowMeans() const {
  std::vector<float> means(rows_);
  for (size_t i = 0; i < rows_; ++i) {
    float sum = 0.0f;
    for (size_t j = 0; j < cols_; ++j) {
      sum += data_[i][j];
    }
    means[i] = sum / static_cast<float>(cols_);
  }
  return means;
}

Matrix Matrix::transpose() const {
  std::vector transposed(cols_, std::vector<float>(rows_));
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      transposed[j][i] = data_[i][j];
    }
  }
  return Matrix(transposed);
}

float Matrix::vectorNorm(const std::vector<float> &vec) {
  float sum = 0.0f;
  for (const float val : vec) {
    sum += val * val;
  }
  return std::sqrt(sum);
}

bool Matrix::isUpperTriangular(const float tolerance) const {
  if (rows_ != cols_)
    return false;

  for (size_t i = 1; i < rows_; ++i) {
    for (size_t j = 0; j < i; ++j) {
      if (std::abs(data_[i][j]) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

std::vector<float> Matrix::getDiagonal() const {
  size_t minDim = std::min(rows_, cols_);
  std::vector<float> diagonal(minDim);
  for (size_t i = 0; i < minDim; ++i) {
    diagonal[i] = data_[i][i];
  }
  return diagonal;
}

std::pair<Matrix, Matrix> Matrix::qrDecomposition() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("QR decomposition requires a square matrix");
  }

  size_t n = rows_;
  Matrix Q = *this;
  Matrix R(n, n, Zeros);

  for (size_t j = 0; j < n; ++j) {
    std::vector<float> v(n);
    for (size_t i = 0; i < n; ++i) {
      v[i] = Q.data_[i][j];
    }

    R.data_[j][j] = vectorNorm(v);

    if (R.data_[j][j] > TOLERANCE) {
      for (size_t i = 0; i < n; ++i) {
        Q.data_[i][j] /= R.data_[j][j];
      }
    }

    for (size_t k = j + 1; k < n; ++k) {
      float sum = 0.0f;
      for (size_t i = 0; i < n; ++i) {
        sum += Q.data_[i][j] * Q.data_[i][k];
      }
      R.data_[j][k] = sum;

      for (size_t i = 0; i < n; ++i) {
        Q.data_[i][k] -= R.data_[j][k] * Q.data_[i][j];
      }
    }
  }

  return std::make_pair(Q, R);
}

std::vector<float> Matrix::calculateEigenvalues(int maxIterations) const {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Eigenvalue calculation requires a square matrix");
  }

  Matrix A = *this;
  int iterations = 0;

  while (!A.isUpperTriangular(TOLERANCE) && iterations < maxIterations) {
    auto [Q, R] = A.qrDecomposition();
    A = R * Q;
    iterations++;
  }

  return A.getDiagonal();
}

Matrix Matrix::pow(const int n) const {
  if (rows_ != cols_) {
    throw std::invalid_argument("Matrix power requires square matrix");
  }
  if (n < 0) {
    throw std::invalid_argument("Negative powers not implemented");
  }
  if (n == 0) {
    return {rows_, cols_, Identity};
  }

  Matrix result = *this;
  for (int i = 1; i < n; ++i) {
    result = result.multiply(*this);
  }
  return result;
}

Matrix Matrix::kroneckerProduct(const Matrix &other) const {
  const size_t resultRows = rows_ * other.rows_;
  const size_t resultCols = cols_ * other.cols_;
  std::vector result(resultRows, std::vector<float>(resultCols));

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      for (size_t k = 0; k < other.rows_; ++k) {
        for (size_t l = 0; l < other.cols_; ++l) {
          result[i * other.rows_ + k][j * other.cols_ + l] =
              data_[i][j] * other.data_[k][l];
        }
      }
    }
  }
  return Matrix(result);
}

std::vector<float> Matrix::colMeans() const {
  std::vector<float> means(cols_, 0.0f);
  for (size_t j = 0; j < cols_; ++j) {
    for (size_t i = 0; i < rows_; ++i) {
      means[j] += data_[i][j];
    }
    means[j] /= static_cast<float>(rows_);
  }
  return means;
}

float Matrix::sum() const {
  float total = 0.0f;
  for (const auto &row : data_) {
    for (const float val : row) {
      total += val;
    }
  }
  return total;
}

float Matrix::mean() const { return sum() / static_cast<float>(rows_ * cols_); }

std::pair<Matrix, Matrix> Matrix::luDecomposition() const {
  if (rows_ != cols_) {
    throw std::invalid_argument("LU decomposition requires square matrix");
  }

  const size_t n = rows_;
  Matrix L(n, n, Zeros);
  Matrix U = *this;

  for (size_t i = 0; i < n; ++i) {
    L.data_[i][i] = 1.0f; // Diagonal elements of L are 1

    for (size_t j = i + 1; j < n; ++j) {
      if (std::abs(U.data_[i][i]) < TOLERANCE) {
        throw std::runtime_error("Zero pivot encountered in LU decomposition");
      }

      const float factor = U.data_[j][i] / U.data_[i][i];
      L.data_[j][i] = factor;

      for (size_t k = i; k < n; ++k) {
        U.data_[j][k] -= factor * U.data_[i][k];
      }
    }
  }

  return std::make_pair(L, U);
}

std::pair<Matrix, Matrix> Matrix::choleskyDecomposition() const {
  if (!isSymmetric() || !isPositiveDefinite()) {
    throw std::invalid_argument(
        "Cholesky decomposition requires symmetric positive definite matrix");
  }

  const size_t n = rows_;
  Matrix L(n, n, Zeros);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      float sum = 0.0f;

      if (j == i) {
        for (size_t k = 0; k < j; ++k) {
          sum += L.data_[j][k] * L.data_[j][k];
        }
        L.data_[j][j] = std::sqrt(data_[j][j] - sum);
      } else {
        for (size_t k = 0; k < j; ++k) {
          sum += L.data_[i][k] * L.data_[j][k];
        }
        L.data_[i][j] = (data_[i][j] - sum) / L.data_[j][j];
      }
    }
  }

  return std::make_pair(L, L.transpose());
}

std::string Matrix::generateHistogram(const std::vector<float> &values) {
  if (values.empty()) {
    return "\\text{\\textbf{No data available for histogram}}";
  }

  // Find min and max values
  const float min_val = *std::ranges::min_element(values);
  float max_val = *std::ranges::max_element(values);

  // Handle case where all values are the same
  if (std::abs(max_val - min_val) < TOLERANCE) {
    max_val = min_val + 1.0f; // Add range to avoid division by zero
  }

  // Create 10 bins for the histogram
  constexpr int num_bins = 10;
  std::vector bins(num_bins, 0);
  const float bin_width = (max_val - min_val) / num_bins;

  // Count values in each bin
  for (const float value : values) {
    const int bin =
        std::min(static_cast<int>((value - min_val) / bin_width), num_bins - 1);
    bins[bin]++;
  }

  // Find maximum bin height for scaling
  const int max_height = *std::ranges::max_element(bins);
  if (max_height == 0) {
    return "\\text{\\textbf{Unable to generate histogram: all bins are empty}}";
  }

  // Generate TikZ picture
  std::stringstream ss;
  ss << "\\begin{tikzpicture}[scale=0.7]\n";

  // Draw axes
  ss << "\\draw[->] (0,0) -- (12,0) node[right] {Value};\n";
  ss << "\\draw[->] (0,0) -- (0,7) node[above] {Frequency};\n";

  // Draw bars
  constexpr float bar_width = 10.0 / num_bins;
  for (int i = 0; i < num_bins; ++i) {
    const float x = static_cast<float>(i) * bar_width;
    const float height =
        bins[i] > 0 ? std::max(0.1f, (static_cast<float>(bins[i]) * 6.0f) /
                                         static_cast<float>(max_height))
                    : 0;

    ss << "\\fill[blue!40] (" << x << ",0) rectangle (" << (x + bar_width)
       << "," << height << ");\n";

    // Add value labels on x-axis
    const float value = min_val + (static_cast<float>(i) * bin_width);
    if (i % 2 == 0) { // Label every other bin for clarity
      ss << "\\node[below] at (" << (x + bar_width / 2) << ",0) {" << std::fixed
         << std::setprecision(1) << value << "};\n";
    }
  }

  ss << "\\end{tikzpicture}";
  return ss.str();
}

std::string Matrix::generateHeatmap() const {
  // Check for empty matrix
  if (rows_ == 0 || cols_ == 0) {
    return "\\text{\\textbf{No data available for heatmap: Matrix is empty}}";
  }

  try {
    // Find min and max values for color scaling
    float min_val = data_.at(0).at(0);
    float max_val = data_.at(0).at(0);

    for (const auto &row : data_) {
      if (row.empty()) {
        return "\\text{\\textbf{Invalid matrix structure: Empty row detected}}";
      }
      for (float val : row) {
        min_val = std::min(min_val, val);
        max_val = std::max(max_val, val);
      }
    }

    // Handle case where all values are the same
    if (std::abs(max_val - min_val) < TOLERANCE) {
      max_val = min_val + 1.0f; // Add range to avoid division by zero
    }

    std::stringstream ss;
    ss << "\\begin{tikzpicture}[scale=0.7]\n";

    // Define color gradient
    ss << "\\definecolor{minColor}{RGB}{240,240,255}\n";
    ss << "\\definecolor{maxColor}{RGB}{0,0,255}\n";

    // Draw cells
    constexpr float cell_size = 1.0;
    for (size_t i = 0; i < rows_; ++i) {
      for (size_t j = 0; j < cols_; ++j) {
        const float value = data_.at(i).at(j); // Use at() for bounds checking
        const float normalized_value = (value - min_val) / (max_val - min_val);

        ss << "\\fill[minColor!" << 100 * normalized_value << "!maxColor] ("
           << static_cast<float>(j) * cell_size << ","
           << static_cast<float>(rows_ - 1 - i) * cell_size << ") rectangle ("
           << static_cast<float>(j + 1) * cell_size << ","
           << static_cast<float>(rows_ - i) * cell_size << ");\n";

        // Add value labels if matrix is small enough
        if (rows_ <= 8 && cols_ <= 8) {
          ss << "\\node at ("
             << static_cast<float>(j) * cell_size + cell_size / 2 << ","
             << static_cast<float>(rows_ - 1 - i) * cell_size + cell_size / 2
             << ") {\\tiny " << std::fixed << std::setprecision(1) << value
             << "};\n";
        }
      }
    }

    // Add row labels
    for (size_t i = 0; i < rows_; ++i) {
      ss << "\\node[left] at (0,"
         << (static_cast<float>(rows_ - i) - 0.5) * cell_size << ") {" << i
         << "};\n";
    }

    // Add column labels
    for (size_t j = 0; j < cols_; ++j) {
      ss << "\\node[below] at ("
         << (static_cast<float>(j) * cell_size + cell_size / 2) << ",0) {" << j
         << "};\n";
    }

    // Add colorbar
    constexpr float colorbar_width = 0.3;
    const float colorbar_height = static_cast<float>(rows_) * cell_size;
    const float colorbar_x = static_cast<float>(cols_ + 1) * cell_size;

    for (int i = 0; i < 100; ++i) {
      const float y = static_cast<float>(i) * colorbar_height / 100.0f;
      ss << "\\fill[minColor!" << i << "!maxColor] (" << colorbar_x << "," << y
         << ") rectangle (" << (colorbar_x + colorbar_width) << ","
         << (y + colorbar_height / 100.0) << ");\n";
    }

    // Add colorbar labels
    ss << "\\node[right] at (" << (colorbar_x + colorbar_width) << ","
       << colorbar_height << ") {" << std::fixed << std::setprecision(1)
       << max_val << "};\n";
    ss << "\\node[right] at (" << (colorbar_x + colorbar_width) << ",0) {"
       << min_val << "};\n";

    ss << "\\end{tikzpicture}";
    return ss.str();
  } catch (const std::exception &e) {
    return "\\text{\\textbf{Error generating heatmap: " +
           std::string(e.what()) + "}}";
  }
}

float Matrix::frobeniusNorm() const {
  float sum = 0.0f;
  for (const auto &row : data_) {
    for (const float val : row) {
      sum += val * val;
    }
  }
  return std::sqrt(sum);
}

float Matrix::maxNorm() const {
  float maxVal = 0.0f;
  for (const auto &row : data_) {
    for (const float val : row) {
      maxVal = std::max(maxVal, std::abs(val));
    }
  }
  return maxVal;
}

bool Matrix::isSymmetric(const float tolerance) const {
  if (rows_ != cols_)
    return false;

  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = i + 1; j < cols_; ++j) {
      if (std::abs(data_[i][j] - data_[j][i]) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

bool Matrix::isOrthogonal(const float tolerance) const {
  if (rows_ != cols_)
    return false;

  // Compute A * A^T
  const Matrix product = multiply(transpose());

  // Check if it's close to identity matrix
  for (size_t i = 0; i < rows_; ++i) {
    for (size_t j = 0; j < cols_; ++j) {
      if (const float expected = (i == j) ? 1.0f : 0.0f;
          std::abs(product.data_[i][j] - expected) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

bool Matrix::isPositiveDefinite() const {
  if (rows_ != cols_)
    return false;

  try {
    // Attempt Cholesky decomposition
    auto [L, _] = choleskyDecomposition();
    return true;
  } catch ([[maybe_unused]] const std::exception &e) {
    return false;
  }
}

std::pair<size_t, size_t> Matrix::rank() const {
  // Simple rank calculation using Gaussian elimination
  Matrix reduced = *this;
  size_t r = 0;

  for (size_t col = 0; col < cols_; ++col) {
    // Find pivot
    bool pivotFound = false;
    for (size_t row = r; row < rows_; ++row) {
      if (std::abs(reduced.data_[row][col]) > 1e-10) {
        // Swap rows
        if (row != r) {
          std::swap(reduced.data_[row], reduced.data_[r]);
        }
        pivotFound = true;
        break;
      }
    }

    if (pivotFound) {
      // Eliminate below
      for (size_t row = r + 1; row < rows_; ++row) {
        const float factor = reduced.data_[row][col] / reduced.data_[r][col];
        for (size_t j = col; j < cols_; ++j) {
          reduced.data_[row][j] -= factor * reduced.data_[r][j];
        }
      }
      r++;
    }
  }

  return {r, cols_ - r};
}

float Matrix::variance() const {
  const float avg = mean();
  float sumSquaredDiff = 0.0f;

  for (const auto &row : data_) {
    for (float val : row) {
      float diff = val - avg;
      sumSquaredDiff += diff * diff;
    }
  }

  return sumSquaredDiff / static_cast<float>(rows_ * cols_);
}

void Matrix::generateLatexReport(const std::string &filename,
                                 const bool DarkMode,
                                 const bool includeAdvancedStats,
                                 const bool includeDecompositions,
                                 const bool includeVisualization) const {
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Could not create LaTeX file");
  }

  std::string mode = "\\pagecolor{bg}\n\\color{text}\n";
  std::string defineColour = "\\definecolor{bg}{RGB}{51,51,51}\n\\definecolor{"
                             "text}{RGB}{220,220,220}\n";
  std::string HeaderColour = DarkMode
                                 ? "\\definecolor{Othertext}{RGB}{251,251,251}"
                                 : "\\definecolor{Othertext}{RGB}{51,51,51}";

  file << "\\documentclass{article}\n"
       << "\\usepackage{color}\n"
       << "\\usepackage{amsmath}\n"
       << "\\usepackage{float}\n"
       << "\\usepackage{tikz}\n"
       << "\\usepackage{pgfplots}\n"
       << "\\pgfplotsset{compat=1.18}\n"
       << "\\usepackage[margin=1in]{geometry}\n"
       << HeaderColour << "\n"
       << (DarkMode ? defineColour : "") << (DarkMode ? mode : "")
       << "\\title{\\color{Othertext}Matrix Analysis Report}\n"
       << "\\author{\\color{Othertext}Generated by Matrix Class}\n"
       << "\\begin{document}\n"
       << "\\date{\\color{Othertext}\\today}\n"
       << "\\maketitle\n"
       << "\n\n";

  file << "\\section{Original Matrix}\n"
       << "\\[\n"
       << matrixToLatex() << "\\]\n\n";

  // Check if matrix is square before attempting square-matrix operations
  const bool isSquare = (rows_ == cols_);

  if (isSquare) {
    // Square matrix operations
    try {
      file << "\\section{Determinant}\n"
           << "The determinant of the matrix is:\n"
           << "\\[\n"
           << "det(A) = " << std::fixed << std::setprecision(3) << determinant()
           << "\\]\n\n";

      file << "\\section{Trace}\n"
           << "The trace of the matrix is:\n"
           << "\\[\n"
           << "Tr(A) = " << std::fixed << std::setprecision(3) << trace()
           << "\\]\n\n";

      file << "\\section{Inverse Matrix}\n"
           << "The inverse of the matrix is:\n"
           << "\\[\n";
      try {
        Matrix inv = inverse();
        file << inv.matrixToLatex();
      } catch (const std::exception &e) {
        file << "\\text{Matrix is not invertible: " << e.what() << "}\n";
      }
      file << "\\]\n\n";

      file << "\\section{Verification}\n"
           << "Multiplying the original matrix with its inverse ($A \\cdot "
              "A^{-1}$) should give the identity matrix:\n"
           << "\\[\n";
      try {
        Matrix inv = inverse();
        Matrix verification = multiply(inv);
        file << verification.matrixToLatex();
      } catch (const std::exception &e) {
        file << "\\text{Verification not possible: " << e.what() << "}\n";
      }
      file << "\\]\n\n";

      file << "\\section{Eigenvalues}\n"
           << "The eigenvalues of the matrix are:\n"
           << "\\[\n";
      try {
        std::vector<float> eigenvals = calculateEigenvalues();
        file << "\\begin{aligned}\n";
        for (size_t i = 0; i < eigenvals.size(); ++i) {
          file << "\\lambda_" << (i + 1) << " &= " << std::fixed
               << std::setprecision(3) << eigenvals[i];
          if (i < eigenvals.size() - 1) {
            file << " \\\\\n";
          }
        }
        file << "\n\\end{aligned}\n";
      } catch (const std::exception &e) {
        file << "\\text{Eigenvalue calculation error: " << e.what() << "}\n";
      }
      file << "\\]\n";
    } catch (const std::exception &e) {
      file << "\\text{Error in square matrix calculations: " << e.what()
           << "}\n";
    }
  } else {
    // Non-square matrix notice
    file << "\\section{Matrix Properties}\n"
         << "This is a " << rows_ << " × " << cols_ << " matrix.\n"
         << "Properties such as determinant, trace, inverse, and eigenvalues "
         << "are only defined for square matrices.\n\n";
  }

  file << "\\section{Null Space Vectors}\n"
       << "The basis vectors of the null space are:\n"
       << "\\[\n";
  try {
    std::vector<std::vector<float>> nullVecs = findNullVectors();
    if (nullVecs.empty()) {
      file << "\\text{The null space contains only the zero vector}\n";
    } else {
      file << "\\begin{aligned}\n";
      for (size_t i = 0; i < nullVecs.size(); ++i) {
        file << "\\vec{v}_" << (i + 1) << " &= \\begin{bmatrix}\n";
        for (size_t j = 0; j < nullVecs[i].size(); ++j) {
          file << std::fixed << std::setprecision(3) << nullVecs[i][j];
          if (j < nullVecs[i].size() - 1) {
            file << " \\\\\n";
          }
        }
        file << "\\end{bmatrix}";
        if (i < nullVecs.size() - 1) {
          file << " \\\\\n";
        }
      }
      file << "\n\\end{aligned}\n";
    }
  } catch (const std::exception &e) {
    file << "\\text{Error calculating null space vectors: " << e.what()
         << "}\n";
  }
  file << "\\]\n\n";

  if (includeAdvancedStats) {
    file << "\\section{Matrix Norms}\n"
         << "\\begin{align*}\n"
         << "\\text{Frobenius Norm} &= " << std::fixed << std::setprecision(3)
         << frobeniusNorm() << "\\\\\n"
         << "\\text{Maximum Norm} &= " << maxNorm() << "\n"
         << "\\end{align*}\n\n";

    file << "\\section{Matrix Properties}\n"
         << "\\begin{itemize}\n"
         << "\\item Symmetric: " << (isSymmetric() ? "Yes" : "No") << "\n"
         << "\\item Orthogonal: " << (isOrthogonal() ? "Yes" : "No") << "\n"
         << "\\item Positive Definite: "
         << (isPositiveDefinite() ? "Yes" : "No") << "\n";

    auto [rank_val, nullity] = rank();
    file << "\\item Rank: " << rank_val << " (Nullity: " << nullity << ")\n"
         << "\\end{itemize}\n\n";

    file << "\\section{Statistical Analysis}\n"
         << "\\subsection{Basic Statistics}\n"
         << "\\begin{align*}\n"
         << "\\text{Mean} &= " << mean() << "\\\\\n"
         << "\\text{Variance} &= " << variance() << "\\\\\n"
         << "\\text{Sum} &= " << sum() << "\n"
         << "\\end{align*}\n\n";
  }

  if (includeDecompositions) {
    file << "\\section{Matrix Decompositions}\n"
         << "\\subsection{LU Decomposition}\n";
    try {
      auto [L, U] = luDecomposition();
      file << "\\[L\n = " << L.matrixToLatex() << "\\]\n\n"
           << "\\[U\n = " << U.matrixToLatex() << "\\]\n\n";
    } catch (const std::exception &e) {
      file << "\\[LU\n decomposition not possible: " << e.what() << "\\]\n\n";
    }

    if (isSymmetric() && isPositiveDefinite()) {
      file << "\\subsection{Cholesky Decomposition}\n";
      try {
        auto [L, Lt] = choleskyDecomposition();
        file << "L = " << L.matrixToLatex() << "\n\n";
      } catch (const std::exception &e) {
        file << "Cholesky decomposition not possible: " << e.what() << "\n\n";
      }
    }
  }

  if (includeVisualization) {
    file << "\\section{Matrix Visualization}\n"
         << "\\subsection{Heatmap}\n"
         << generateHeatmap() << "\n\n"
         << "\\subsection{Row Means Distribution}\n"
         << generateHistogram(rowMeans()) << "\n\n";
  }

  file << "\\end{document}\n";
  file.close();
}

// std::vector<std::vector<float>> Matrix::findNullVectors() const {
//   // Use Singular Value Decomposition (SVD) approach for better numerical
//   // stability This is more reliable than RREF for ill-conditioned matrices

//   // Helper function for Householder reflection
//   auto householderReflection = [](const std::vector<float> &x, size_t start)
//   {
//     float norm = 0.0f;
//     for (size_t i = start; i < x.size(); ++i) {
//       norm += x[i] * x[i];
//     }
//     norm = std::sqrt(norm);

//     if (norm < TOLERANCE)
//       return std::vector<float>(x.size(), 0.0f);

//     std::vector<float> v(x);
//     v[start] += (v[start] > 0 ? norm : -norm);

//     float scale = 0.0f;
//     for (size_t i = start; i < v.size(); ++i) {
//       scale += v[i] * v[i];
//     }
//     scale = 2.0f / scale;

//     return v;
//   };

//   // Create working copy with increased precision
//   Matrix A(rows_, cols_);
//   for (size_t i = 0; i < rows_; ++i) {
//     for (size_t j = 0; j < cols_; ++j) {
//       A.data_[i][j] = data_[i][j];
//     }
//   }

//   // Compute QR decomposition using Householder reflections
//   std::vector<std::vector<float>> reflectors;
//   const size_t min_dim = std::min(rows_, cols_);

//   for (size_t j = 0; j < min_dim; ++j) {
//     std::vector<float> column(rows_);
//     for (size_t i = j; i < rows_; ++i) {
//       column[i] = A.data_[i][j];
//     }

//     auto v = householderReflection(column, j);
//     if (v[j] == 0.0f)
//       continue;

//     reflectors.push_back(v);

//     // Apply reflection to remaining columns
//     for (size_t k = j; k < cols_; ++k) {
//       float dot = 0.0f;
//       for (size_t i = j; i < rows_; ++i) {
//         dot += v[i] * A.data_[i][k];
//       }
//       dot *= 2.0f / (vectorNorm(v) * vectorNorm(v));

//       for (size_t i = j; i < rows_; ++i) {
//         A.data_[i][k] -= dot * v[i];
//       }
//     }
//   }

//   // Find numerical rank using SVD-like tolerance
//   float max_diag = 0.0f;
//   for (size_t i = 0; i < min_dim; ++i) {
//     max_diag = std::max(max_diag, std::abs(A.data_[i][i]));
//   }

//   size_t rank = 0;
//   const float rank_tolerance = max_diag * std::max(rows_, cols_) * TOLERANCE;

//   for (size_t i = 0; i < min_dim; ++i) {
//     if (std::abs(A.data_[i][i]) > rank_tolerance) {
//       rank++;
//     }
//   }

//   // Initialize null space basis
//   std::vector<std::vector<float>> nullVectors;
//   if (rank >= cols_)
//     return nullVectors; // No null space

//   // For each free variable (cols - rank of them)
//   for (size_t j = rank; j < cols_; ++j) {
//     std::vector<float> nullVector(cols_, 0.0f);
//     nullVector[j] = 1.0f;

//     // Back-substitution with improved numerical stability
//     for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
//       float sum = 0.0f;
//       for (size_t k = i + 1; k < cols_; ++k) {
//         sum += A.data_[i][k] * nullVector[k];
//       }

//       if (std::abs(A.data_[i][i]) > rank_tolerance) {
//         nullVector[i] = -sum / A.data_[i][i];
//       } else {
//         nullVector[i] = 0.0f;
//       }
//     }

//     // Apply Gram-Schmidt orthogonalization to improve numerical stability
//     for (const auto &existing : nullVectors) {
//       float dot = 0.0f;
//       for (size_t k = 0; k < cols_; ++k) {
//         dot += nullVector[k] * existing[k];
//       }
//       for (size_t k = 0; k < cols_; ++k) {
//         nullVector[k] -= dot * existing[k];
//       }
//     }

//     // Normalize the vector
//     float norm = vectorNorm(nullVector);
//     if (norm > TOLERANCE) {
//       for (float &component : nullVector) {
//         component /= norm;
//       }
//       nullVectors.push_back(nullVector);
//     }
//   }

//   // Additional validation and refinement
//   if (!nullVectors.empty()) {
//     // Verify orthogonality to row space
//     for (auto &nullVector : nullVectors) {
//       std::vector<float> residual(rows_, 0.0f);
//       for (size_t i = 0; i < rows_; ++i) {
//         for (size_t j = 0; j < cols_; ++j) {
//           residual[i] += data_[i][j] * nullVector[j];
//         }
//       }

//       // Iterative refinement if needed
//       float residualNorm = vectorNorm(residual);
//       if (residualNorm > TOLERANCE * 10.0f) {
//         std::vector<float> correction(cols_, 0.0f);
//         for (size_t j = 0; j < cols_; ++j) {
//           for (size_t i = 0; i < rows_; ++i) {
//             correction[j] += data_[i][j] * residual[i];
//           }
//         }

//         for (size_t j = 0; j < cols_; ++j) {
//           nullVector[j] -= correction[j] / (residualNorm * residualNorm);
//         }

//         // Renormalize after refinement
//         float norm = vectorNorm(nullVector);
//         if (norm > TOLERANCE) {
//           for (float &component : nullVector) {
//             component /= norm;
//           }
//         }
//       }
//     }
//   }

//   return nullVectors;
// }

std::vector<std::vector<float>> Matrix::findNullVectors() const {
  //--------------------------------------------------------------------------
  // 1) Check matrix dimensions.
  //--------------------------------------------------------------------------
  size_t rows = getRows();
  size_t cols = getCols();
  if (rows == 0 || cols == 0) {
    return {};
  }

  //--------------------------------------------------------------------------
  // 2) Make a local copy for in-place pivoted QR factorization.
  //    We'll store data_ in a nested std::vector<float> for convenience.
  //--------------------------------------------------------------------------
  std::vector<std::vector<float>> M(rows, std::vector<float>(cols));
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      M[i][j] = at(i, j);
    }
  }

  //--------------------------------------------------------------------------
  // 3) Set up a pivot array (to track column permutations) and column norms.
  //--------------------------------------------------------------------------
  std::vector<size_t> pivots(cols);
  for (size_t j = 0; j < cols; ++j) {
    pivots[j] = j;
  }

  std::vector<float> colNorm(cols, 0.0f);
  for (size_t j = 0; j < cols; ++j) {
    double sumSq = 0.0;
    for (size_t i = 0; i < rows; ++i) {
      sumSq += double(M[i][j]) * double(M[i][j]);
    }
    colNorm[j] = static_cast<float>(std::sqrt(sumSq));
  }

  //--------------------------------------------------------------------------
  // Helper to compute Euclidean norm of M[k..end][col]
  //--------------------------------------------------------------------------
  auto columnSegmentNorm = [&](size_t k, size_t col) {
    double s = 0.0;
    for (size_t i = k; i < rows; ++i) {
      s += double(M[i][col]) * double(M[i][col]);
    }
    return static_cast<float>(std::sqrt(s));
  };

  //--------------------------------------------------------------------------
  // 4) Perform pivoted Householder QR factorization in-place.
  //--------------------------------------------------------------------------
  size_t minRC = std::min(rows, cols);
  constexpr float TOL = TOLERANCE; // usage of the class-defined tolerance

  for (size_t k = 0; k < minRC; ++k) {
    // 4a) Pivot: select the column j >= k with the largest norm
    size_t pivotCol = k;
    float maxNorm = colNorm[k];
    for (size_t j = k + 1; j < cols; ++j) {
      if (colNorm[j] > maxNorm) {
        maxNorm = colNorm[j];
        pivotCol = j;
      }
    }

    // 4b) Swap the chosen pivot column with column k
    if (pivotCol != k) {
      std::swap(pivots[pivotCol], pivots[k]);
      std::swap(colNorm[pivotCol], colNorm[k]);
      for (size_t i = 0; i < rows; ++i) {
        std::swap(M[i][pivotCol], M[i][k]);
      }
    }

    // 4c) Form the Householder vector to eliminate sub-diagonal elements in
    // column k
    float normX = columnSegmentNorm(k, k);
    if (normX < TOL) {
      // If the column is too small, skip reflection
      continue;
    }
    float alpha = (M[k][k] > 0.f) ? -normX : normX;
    float v0 = M[k][k] - alpha;
    M[k][k] = alpha; // on diagonal

    // 4d) Apply reflection to the trailing columns [k+1..cols-1]
    //     Using the simplified Householder formula:
    //     H = I - (v * v^T) / (v^T v), where v = [v0 M[k+1][k] ...
    //     M[rows-1][k]]^T We'll compute dot = v^T * nextCol, then scale, then
    //     subtract.
    for (size_t j = k + 1; j < cols; ++j) {
      double dot = double(v0) * double(M[k][j]);
      for (size_t i = k + 1; i < rows; ++i) {
        dot += double(M[i][k]) * double(M[i][j]);
      }
      double denom = double(v0) * double(v0);
      for (size_t i = k + 1; i < rows; ++i) {
        denom += double(M[i][k]) * double(M[i][k]);
      }
      if (std::fabs(denom) < 1e-32) {
        continue; // avoid dividing by zero
      }
      double scale = dot / denom;
      M[k][j] -= static_cast<float>(scale * double(v0));
      for (size_t i = k + 1; i < rows; ++i) {
        M[i][j] -= static_cast<float>(scale * double(M[i][k]));
      }
    }

    // 4e) Update the norms for columns k+1..cols-1
    for (size_t j = k + 1; j < cols; ++j) {
      double s = 0.0;
      for (size_t i = k; i < rows; ++i) {
        s += double(M[i][j]) * double(M[i][j]);
      }
      colNorm[j] = static_cast<float>(std::sqrt(s));
    }
  }

  //--------------------------------------------------------------------------
  // 5) Determine rank from the diagonal of R
  //--------------------------------------------------------------------------
  float maxDiag = 0.0f;
  for (size_t i = 0; i < minRC; ++i) {
    float val = std::fabs(M[i][i]);
    if (val > maxDiag) {
      maxDiag = val;
    }
  }
  float rankThreshold = maxDiag * float(std::max(rows, cols)) * TOL;
  size_t rank = 0;
  for (size_t i = 0; i < minRC; ++i) {
    if (std::fabs(M[i][i]) > rankThreshold) {
      rank++;
    } else {
      break;
    }
  }

  //--------------------------------------------------------------------------
  // 6) If rank >= cols, no null space
  //--------------------------------------------------------------------------
  if (rank >= cols) {
    return {};
  }

  //--------------------------------------------------------------------------
  // 7) Back-substitution to find each null-space vector
  //--------------------------------------------------------------------------
  std::vector<std::vector<float>> nullSpace;
  nullSpace.reserve(cols - rank);

  for (size_t freeIdx = rank; freeIdx < cols; ++freeIdx) {
    // We'll solve for xPrime in the pivoted system, then un-permute to original
    std::vector<float> xPrime(cols, 0.0f);
    xPrime[freeIdx] = 1.0f;

    // R(0..rank-1, 0..rank-1) * x(0..rank-1) = -R(0..rank-1, freeIdx) =>
    // back-substitute
    for (int i = int(rank) - 1; i >= 0; --i) {
      float rhs = -M[size_t(i)][freeIdx];
      for (size_t c = size_t(i) + 1; c < rank; ++c) {
        rhs -= M[size_t(i)][c] * xPrime[c];
      }
      float diagVal = M[size_t(i)][size_t(i)];
      if (std::fabs(diagVal) > TOL) {
        xPrime[size_t(i)] = rhs / diagVal;
      } else {
        xPrime[size_t(i)] = 0.0f;
      }
    }

    // Un-permute xPrime back to the original column ordering
    std::vector<float> basis(cols, 0.0f);
    for (size_t c = 0; c < cols; ++c) {
      basis[pivots[c]] = xPrime[c];
    }

    nullSpace.push_back(std::move(basis));
  }

  //--------------------------------------------------------------------------
  // 8) Optional: normalize each basis vector
  //--------------------------------------------------------------------------
  for (auto &v : nullSpace) {
    double sumSq = 0.0;
    for (float val : v) {
      sumSq += double(val) * double(val);
    }
    float nrm = static_cast<float>(std::sqrt(sumSq));
    if (nrm > TOL) {
      // Use the class's built-in vectorNorm if you want:
      // float nrm = vectorNorm(v);
      for (float &val : v) {
        val /= nrm;
      }
    }
  }

  return nullSpace;
}
