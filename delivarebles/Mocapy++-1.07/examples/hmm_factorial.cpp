/*
 * hmm_factorial.cpp
 *
 *  Copyright (C) 2008, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mocapy.h"
using namespace mocapy;

int main(void) {
   	mocapy_seed((uint)5556574);

	// HMM hidden and observed node sizes
	uint H_SIZE=5;
	uint O_SIZE=3;

	Node* th1 = NodeFactory::new_discrete_node(H_SIZE, "th1");
	Node* th2 = NodeFactory::new_discrete_node(H_SIZE, "th2");
	Node* th3 = NodeFactory::new_discrete_node(H_SIZE, "th3");

	Node* th01 = NodeFactory::new_discrete_node(H_SIZE, "th01");
	Node* th02 = NodeFactory::new_discrete_node(H_SIZE, "th02");
	Node* th03 = NodeFactory::new_discrete_node(H_SIZE, "th03");

	Node* to = NodeFactory::new_discrete_node(O_SIZE, "to");

	DBN tdbn;
	tdbn.set_slices(vec(th01, th02, th03, to), vec(th1, th2, th3, to));

	tdbn.add_intra("th1", "to");
	tdbn.add_intra("th2", "to");
	tdbn.add_intra("th3", "to");
	tdbn.add_inter("th1", "th01");
	tdbn.add_inter("th2", "th02");
	tdbn.add_inter("th3", "th03");
	tdbn.construct();

	return EXIT_SUCCESS;
}
