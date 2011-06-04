/*
 * utils.cpp
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
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <boost/random.hpp>
#include "utils.h"
#include <cstring>
#include <sys/stat.h> 

using namespace std;

namespace mocapy {

uint moc_seed=251177;

void mocapy_seed(uint s) {
  moc_seed=s;
}



vector<Node*> vec_concNode(vector<Node*> & v1, vector<Node*> & v2) {
	vector<Node*> v3;
	v3 = v1;
	for (uint i = 0; i < v2.size(); i++) {
		v3.push_back(v2[i]);
	}
	return v3;
}

/* reverse:  reverse string s in place */
void reverse(char s[])
{
    int c, i, j;

    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

/* itoa:  convert n to characters in s */
void itoa(int n, char s[]) {
	int i, sign;

	if ((sign = n) < 0) /* record sign */
		n = -n; /* make n positive */
	i = 0;
	do { /* generate digits in reverse order */
		s[i++] = n % 10 + '0'; /* get next digit */
	} while ((n /= 10) > 0); /* delete it */
	if (sign < 0)
		s[i++] = '-';
	s[i] = '\0';
	reverse(s);
}

typedef union {
	long L;
	float F;
} LF_t;

char *ftoa(float f, int *status) {
	long mantissa, int_part, frac_part;
	short exp2;
	LF_t x;
	char *p;
	static char outbuf[15];

	*status = 0;
	if (f == 0.0) {
		outbuf[0] = '0';
		outbuf[1] = '.';
		outbuf[2] = '0';
		outbuf[3] = 0;
		return outbuf;
	}
	x.F = f;

	exp2 = (unsigned char) (x.L >> 23) - 127;
	mantissa = (x.L & 0xFFFFFF) | 0x800000;
	frac_part = 0;
	int_part = 0;

	if (exp2 >= 31) {
		*status = 100;
		return 0;
	} else if (exp2 < -23) {
		*status = 101;
		return 0;
	} else if (exp2 >= 23)
		int_part = mantissa << (exp2 - 23);
	else if (exp2 >= 0) {
		int_part = mantissa >> (23 - exp2);
		frac_part = (mantissa << (exp2 + 1)) & 0xFFFFFF;
	} else
		/* if (exp2 < 0) */
		frac_part = (mantissa & 0xFFFFFF) >> -(exp2 + 1);

	p = outbuf;

	if (x.L < 0)
		*p++ = '-';

	if (int_part == 0)
		*p++ = '0';
	else {
		mocapy::itoa(int_part, p);
		while (*p)
			p++;
	}
	*p++ = '.';

	if (frac_part == 0)
		*p++ = '0';
	else {
		char m, max;

		max = sizeof(outbuf) - (p - outbuf) - 1;
		if (max > 7)
			max = 7;
		/* print BCD */
		for (m = 0; m < max; m++) {
			/* frac_part *= 10; */
			frac_part = (frac_part << 3) + (frac_part << 1);

			*p++ = (frac_part >> 24) + '0';
			frac_part &= 0xFFFFFF;
		}
		/* delete ending zeroes */
		for (--p; p[0] == '0' && p[-1] != '.'; --p)
			;
		++p;
	}
	*p = 0;

	return outbuf;
}

vector<char*> get_tokens(char* line, char sep) {
	uint MAX_LINE_SIZE(2048);

	vector<char*> tokens;
	char current[2048];
	int j(0);
	for (uint i=0; i<MAX_LINE_SIZE; i++) {
		char c = line[i];
//		cout << "c: " << c << " i: " << i << endl;
		if (c == sep || (c==0 && i>0) ) {
			if (j>0) {
				current[j++]=0;
				char* save = new char[j];
				strcpy(save, current);
				tokens.push_back(save);
				j=0;
			}
		}
		else
			current[j++] = c;

		if (c==0)
			break;


	}
	return tokens;
}


vector<MDArray<double> > data_loader(const char* filename, char sep, vector<uint> columns) {
	ifstream f;
	f.open(filename);
	if (!f) {
		cerr << "Could not open file " << filename << "\n" ;
		exit(0);
	}

	vector<MDArray<double> > seq_list;
	char line[2048];

	vector<vector<double> > sequence;
	while (f.getline(line, 2048)) {
		if (line[0] == '#')
			continue;

		vector<double> slice;

		vector<char*> tokens = get_tokens(line, sep);
/*
		for (uint i=0; i<tokens.size(); i++) {
			cout << "tok: " << i << " " << tokens[i] << endl;
		}
		*/

		if (tokens.empty()) {
			// new line means new sequence
			if (!sequence.empty()) {
				MDArray<double> newSeq = toMDArray(sequence);
				seq_list.push_back(newSeq);
				sequence.clear();
			}
			continue;

		}

		// Do selection and conversions
		vector<char*> selected;
		if (columns.empty()) {
			selected = tokens;
		}
		else {
			selected = take(tokens, columns);
		}

		for (uint i=0; i<selected.size(); i++) {
			// convert...
			slice.push_back(atof(selected[i]));
		}

		sequence.push_back(slice);

		// Clean up tokens
		for (uint i=0; i<tokens.size(); i++) {
			delete[] tokens[i];
		}
	}
    f.close();

	if (!sequence.empty()) {
		MDArray<double> newSeq = toMDArray(sequence);
		seq_list.push_back(newSeq);
	}

	return seq_list;
}



vector<MDArray<eMISMASK> > toMismask(vector<MDArray<double> > data) {
	vector<MDArray<eMISMASK> > mda(data.size());
	for (uint i = 0; i<data.size(); i++) {
		mda[i].copy_cast(data[i]);
	}
	return mda;
}



bool FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}




}


// This is for automake
void test_library() {
	cout << "Library OK" << endl;
}
