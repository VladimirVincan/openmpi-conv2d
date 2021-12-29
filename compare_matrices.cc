#include <fstream>
#include <iostream>

#include "colors.hpp"

using namespace std;

int main() {
  const char distributed_file_name[] = "outputs/distributed_matrix.txt";
  const char serial_file_name[] = "outputs/serial_matrix.txt";
  const char regular_distributed_file_name[] =
      "outputs/regular_distributed_matrix.txt";
  const char regular_serial_file_name[] = "outputs/regular_serial_matrix.txt";
  const char python_file_name[] = "outputs/python_matrix.txt";

  ifstream D(distributed_file_name);
  ifstream S(serial_file_name);
  // ifstream RD(regular_distributed_file_name);
  // ifstream RS(regular_serial_file_name);
  ifstream P(python_file_name);

  int i = 0;
  double d, s, p, rd, rs;

  while (true) {
    D >> d;
    S >> s;
    // RD >> rd;
    // RS >> rs;
    P >> p;
    if (p != d || p != s ) { // || p != rd || p != rs) {
      cout << RED << "MATRICES NOT EQUAL." << RESET << " i = " << i << endl;
      cout << RED << "d = " << d << endl;
      cout << RED << "s = " << s << endl;
      // cout << RED << "rd = " << rd << endl;
      // cout << RED << "rs = " << rs << endl;
      cout << RED << "p = " << p << endl;

      D.close();
      S.close();
      // RD.close();
      // RS.close();
      P.close();
      return 1;
    }
    ++i;

    if (S.eof() && D.eof() && P.eof()) // && RD.eof() && RS.eof())
      break;
    else if (S.eof() || D.eof() || P.eof()) { // || RD.eof() || RS.eof()) {
      cout << RED << "FILE SIZE NOT EQUAL." << RESET << endl;
      D.close();
      S.close();
      // RD.close();
      // RS.close();
      P.close();
      return 2;
    }
  }

  // cout << GREEN << "EQUAL" << RESET << endl;
  D.close();
  S.close();
  // RD.close();
  // RS.close();
  P.close();
  return 0;
}
