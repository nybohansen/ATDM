/*
 * hmm_mixed.cpp
 *
 *  Created on: 16. jun 2011
 *      Author: Kasper Nybo Hansen
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mocapy.h"
using namespace mocapy;
using namespace std;

int main(void) {

	DBN mdbn;

	mdbn.load("mixed_hmm2.dbn");

    std::vector<Node*> nodes = mdbn.getNodes0(); 
    
    Node *h1 = nodes[0]; 
    Node *o1 = nodes[1];     

	cout << "*** LOADED MODEL WITH PARAMETERS ***" << endl;
    cout << "h1: " << *h1 << endl;
    cout << "o1: " << *o1 << endl;

	return EXIT_SUCCESS;
}
