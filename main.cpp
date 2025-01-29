#include "Matrix.hpp"
#include <iostream>

int main() {
  try {
    Mat2D values = {{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}, {7.0f, 8.0f, 9.0f}};
    Mat2D vals = {{2, 4, 6, 8}, {1, 3, 0, 5}, {1, 1, 6, 3}};
    const Matrix testMatrix(vals);

    if (bool success = testMatrix.generatePDF("matrix")) {
      std::cout << "PDF generated successfully!" << std::endl;
    } else {
      std::cout << "Failed to generate PDF." << std::endl;
    }
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}
