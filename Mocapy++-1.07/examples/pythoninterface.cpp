#include "mocapy.h"

// Global DBN object
mocapy::DBN* dbn;

void makeDBN() {
  dbn = new mocapy::DBN();
  std::cout << "New DBN created" << std::endl;
}

void loadDBN(char* filename) {
  dbn->load(filename);
  std::cout << "DBN loaded" << std::endl;
}

std::vector<double> sample(int i) {
  std::pair<mocapy::Sequence, double> r = dbn->sample_sequence(i);
  return r.first.get_values_flat();
}


int main() {
}
