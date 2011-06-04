/*
This file shows some of the things you can do with Mocapy++'s MDArray.

The MDArray is a Multi-Dimensional array with many features. One feature is
that the shape can be modified at run time.
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mocapy.h"
using namespace mocapy;
using namespace std;

int main(void) {
  // Create a 2x2 MDArray with doubles initialized as zeros
  MDArray<double> myMDA01(vec(2,2));

  // Output the values of the MDArray
  cout << myMDA01 << endl;

  // Set index 0,0 to 1
  myMDA01.set(0,0,1);

  // Set index [1,1] to 5
  myMDA01.set(1,1,5);

  // Output the values of the MDArray
  cout << myMDA01 << endl;

  // Get value of index 0,0
  cout << myMDA01.get(0,0) << endl;

  // Create a 2x3x3 MDArray with doubles initialized as zeros
  MDArray<double> myMDA02(vec(2,3,3));

  // Output the values of the MDArray
  cout << myMDA02 << endl;

  // Set index [1,2,1] to 100. Note that you have to use the [] syntax when there are more than two dimensions.
  vector<uint> index = vec((uint)1, (uint)2, (uint)1);
  myMDA02[index] = 100;

  // Output the values of the MDArray
  cout << myMDA02 << endl;

  // Get the second subarray. The dimensions of the subarray are 3x3
  MDArray<double> sub =  myMDA02.get_view(1);
  cout << sub << endl;

  // Make another 3x3 array with random values
  MDArray<double> myMDA03(vec(3,3));
  RandomGen rg;
  myMDA03.randomize(&rg, 10);

  // Output the values of the MDArray
  cout << myMDA03 << endl;

  // Subtract the sub array with myMDA03
  MDArray<double> myMDA04 = myMDA03 - sub;

  // Output the values of the MDArray
  cout << myMDA04 << endl;

  // Arithmetics are faster if you can do them in place (overwrite the old values)
  myMDA03.add_inplace(sub);
  cout << myMDA03 << endl;

  // Compute determinant
  cout << myMDA03.det() << endl;
}
