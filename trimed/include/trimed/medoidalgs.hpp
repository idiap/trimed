/* trimed, a C++ library and (optional) Python tool for obtaining the
*  medoid of a set
*  
*  Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
*  Written by James Newling <jnewling@idiap.ch>
*  
*  This file is part of trimed.
*  
*  trimed is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License version 3 as
*  published by the Free Software Foundation.
*  
*  trimed is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*  GNU General Public License for more details.
*  
*  You should have received a copy of the GNU General Public License
*  along with Foobar. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TRIMEDOID_H
#define TRIMEDOID_H

#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <trimed/arraydata.hpp>
#include <trimed/arrutilv2minmax.hpp>
#include <trimed/dugraph.hpp>
#include <trimed/uugraph.hpp>
/* where Data has public,
 * (1) ndata : number of datapoints
 * (2) set_distances(size_t i, TFloat * const distances) :
 * 0 <= i < ndata, after call distances[j] = \|data[i] - data[j]\| */
namespace shuffle {
std::unique_ptr<size_t[]> get_up_shuffled_range(size_t r0, size_t r1);
}

namespace medoidalgs {

template <class Data, typename TFloat>
/* Using a RAND type algorithm as per Eppstein */
void RAND(Data &data, TFloat &min_E, size_t &min_E_index, size_t maxcomputes);

template <class Data, typename TFloat>
/* Using the original triangle algorithm */
void trimed(Data &data, TFloat &min_E, size_t &min_E_i);

template <class Data, typename TFloat>
/* The TOPRANK1 algorithm of Okamoto with k = 1 */
void TOPRANK1(Data &data, TFloat &min_E, size_t &min_E_i);

template <class Data, typename TFloat>
/* The TOPRANK2 algorithm of Okamoto with k = 1 */
void TOPRANK2(Data &data, TFloat &min_E, size_t &min_E_i);

template <class Data, typename TFloat>
void medoid(Data &data, TFloat &min_E, size_t &min_E_i, std::string algorithm,
            bool capture_output, std::string &text, size_t seed,
            size_t maxcomputes);

template <typename TFloat>
void medoid_array(size_t ndata, size_t dimension, const TFloat *const data,
                  TFloat &min_E, size_t &min_E_i, std::string algorithm,
                  bool capture_output, std::string &text, size_t seed,
                  size_t maxcomputes) {
  arraydata::ArrayData<TFloat> arda(ndata, dimension, data);
  medoid(arda, min_E, min_E_i, algorithm, capture_output, text, seed,
         maxcomputes);
}

/* finds medoid of largest connected component of undirected unweighted graph */
void medoid_uugraph(size_t nedges_in, const size_t *const edges_in,
                    double &min_E, size_t &min_E_i, std::string algorithm,
                    bool capture_output, std::string &text, size_t seed,
                    size_t maxcomputes);

/* finds medoid of largest connected component of directed unweighted graph */
void medoid_dugraph(size_t nedges_in, const size_t *const edges_in,
                    double &min_E, size_t &min_E_i, std::string algorithm,
                    bool capture_output, std::string &text, size_t seed,
                    size_t maxcomputes);
}

#endif
