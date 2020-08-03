#include "common_tools.h"
#include <random>

using namespace std;


std::complex<double> randomComplex()
{
  const double real_part = randomReal();
  const double imag_part = randomReal();
  const std::complex<double> result(real_part, imag_part);
  return result;
}

double randomReal()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> distrib(0, 1);
  const int num = distrib(gen);
  const double result = static_cast<double>(num) - 0.5;
  return result;
}

