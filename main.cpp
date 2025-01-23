#include "Matrix.hpp"
#include <iostream>

int main() {
  try {
    Mat2D matrix = {{-1.0, 2.0, 3.0, 4.0, 5.0},
                    {6.0, 7.0, 8.0, 9.0, 10.0},
                    {11.0, 12.1, 13.4, 14.4, 15.8},
                    {16.9, 17.8, 18.0, 19.7, 20.6},
                    {21.3, 22.2, 23.6, 24.1, 25.0}};

    const Matrix testMatrix(matrix);

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
