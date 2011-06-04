/*
 *  logfactorial.cpp
 *
 *  Copyright (C) 2008, Thomas Hamelryck, The Bioinformatics Centre,
 *  University of Copenhagen.
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



#include "LogFactorial.h"
#include <boost/math/special_functions/gamma.hpp>

namespace mocapy
{

  std::vector<double> *LogFactorial::cache=0;

LogFactorial::LogFactorial()
{
    if (!cache)
    {
      cache = new std::vector<double>(LOGFACTORIAL_MAXARG);

        // Fill the log factorial array
        for(uint i=0; i<LOGFACTORIAL_MAXARG; i++)
        {
            (*cache)[i]=boost::math::lgamma(i+1);
        }
    }
}

double LogFactorial::get(uint i)
{
    if(i<LOGFACTORIAL_MAXARG)
        return (*cache)[i];
    else
        return boost::math::lgamma(i+1);
}

}

