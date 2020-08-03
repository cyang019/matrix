#include "cxdbl_mat_op.h"
#include <iostream>

using namespace std;

int main()
{
  mat_mat_mul(4, 50000);
  mat_mat_mul(4, 50000);
  mat_mat_mul(4, 50000);
  cout << "..." << endl;
  mat_mat_mul(8, 10000);
  mat_mat_mul(8, 10000);
  mat_mat_mul(8, 10000);
  cout << "..." << endl;
  mat_mat_mul(16, 10000);
  mat_mat_mul(16, 10000);
  mat_mat_mul(16, 10000);
  cout << "..." << endl;
  mat_mat_mul(64, 5000);
  mat_mat_mul(64, 5000);
  mat_mat_mul(64, 5000);
  cout << "..." << endl;
  mat_mat_mul(256, 1000);
  mat_mat_mul(256, 1000);
  mat_mat_mul(256, 1000);

  return 0;
}
