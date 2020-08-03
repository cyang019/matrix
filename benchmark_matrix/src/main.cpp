#include "cxdbl_mat_op.h"
#include <iostream>
#include <cmath>

using namespace std;

int main()
{
  double t1 = mat_mat_mul(4, 80000);
  double t2 = mat_mat_mul(4, 80000);
  double t3 = mat_mat_mul(4, 80000);
  double t4 = mat_mat_mul(4, 80000);
  double t5 = mat_mat_mul(4, 80000);
  double total_time = t1 + t2 + t3 + t4 + t5;
  double mean_time = total_time / 5.0;
  double var_time = 
    (t1 - mean_time) * (t1 - mean_time) +
    (t2 - mean_time) * (t2 - mean_time) +
    (t3 - mean_time) * (t3 - mean_time) +
    (t4 - mean_time) * (t4 - mean_time) +
    (t5 - mean_time) * (t5 - mean_time);
  double std_time = sqrt(var_time/4);
  cout << "avg: " << mean_time << endl;
  cout << "std: " << std_time << endl;

  cout << "..." << endl;
  t1 = mat_mat_mul(8, 40000);
  t2 = mat_mat_mul(8, 40000);
  t3 = mat_mat_mul(8, 40000);
  t4 = mat_mat_mul(8, 40000);
  t5 = mat_mat_mul(8, 40000);
  total_time = t1 + t2 + t3 + t4 + t5;
  mean_time = total_time / 5.0;
  var_time = 
    (t1 - mean_time) * (t1 - mean_time) +
    (t2 - mean_time) * (t2 - mean_time) +
    (t3 - mean_time) * (t3 - mean_time) +
    (t4 - mean_time) * (t4 - mean_time) +
    (t5 - mean_time) * (t5 - mean_time);
  std_time = sqrt(var_time/4);
  cout << "avg: " << mean_time << endl;
  cout << "std: " << std_time << endl;

  cout << "..." << endl;
  t1 = mat_mat_mul(16, 20000);
  t2 = mat_mat_mul(16, 20000);
  t3 = mat_mat_mul(16, 20000);
  t4 = mat_mat_mul(16, 20000);
  t5 = mat_mat_mul(16, 20000);
  total_time = t1 + t2 + t3 + t4 + t5;
  mean_time = total_time / 5.0;
  var_time = 
    (t1 - mean_time) * (t1 - mean_time) +
    (t2 - mean_time) * (t2 - mean_time) +
    (t3 - mean_time) * (t3 - mean_time) +
    (t4 - mean_time) * (t4 - mean_time) +
    (t5 - mean_time) * (t5 - mean_time);
  std_time = sqrt(var_time/4);
  cout << "avg: " << mean_time << endl;
  cout << "std: " << std_time << endl;

  cout << "..." << endl;
  t1 = mat_mat_mul(64, 5000);
  t2 = mat_mat_mul(64, 5000);
  t3 = mat_mat_mul(64, 5000);
  t4 = mat_mat_mul(64, 5000);
  t5 = mat_mat_mul(64, 5000);
  total_time = t1 + t2 + t3 + t4 + t5;
  mean_time = total_time / 5.0;
  var_time = 
    (t1 - mean_time) * (t1 - mean_time) +
    (t2 - mean_time) * (t2 - mean_time) +
    (t3 - mean_time) * (t3 - mean_time) +
    (t4 - mean_time) * (t4 - mean_time) +
    (t5 - mean_time) * (t5 - mean_time);
  std_time = sqrt(var_time/4);
  cout << "avg: " << mean_time << endl;
  cout << "std: " << std_time << endl;
  cout << "..." << endl;
  t1 = mat_mat_mul(256, 1000);
  t2 = mat_mat_mul(256, 1000);
  t3 = mat_mat_mul(256, 1000);
  t4 = mat_mat_mul(256, 1000);
  t5 = mat_mat_mul(256, 1000);
  total_time = t1 + t2 + t3 + t4 + t5;
  mean_time = total_time / 5.0;
  var_time = 
    (t1 - mean_time) * (t1 - mean_time) +
    (t2 - mean_time) * (t2 - mean_time) +
    (t3 - mean_time) * (t3 - mean_time) +
    (t4 - mean_time) * (t4 - mean_time) +
    (t5 - mean_time) * (t5 - mean_time);
  std_time = sqrt(var_time/4);
  cout << "avg: " << mean_time << endl;
  cout << "std: " << std_time << endl;

  return 0;
}
