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

#ifndef ARRAYDATA_HPP
#define ARRAYDATA_HPP

#include <memory>
#include <trimed/arrutilv2l0.hpp>
#include <trimed/arrutilv2l1.hpp>
#include <trimed/arrutilv2l2.hpp>

namespace arraydata {
// special case for array of data.
template <typename TFloat> class ArrayData {
public:
  size_t ndata;
  size_t dimension;
  const TFloat *data;
  std::unique_ptr<TFloat[]> up_data_l22s;
  TFloat *data_l22s;

  ArrayData(size_t ndatai, size_t dimensioni, const TFloat *const datai)
      : ndata(ndatai), dimension(dimensioni), data(datai) {
    up_data_l22s.reset(new TFloat[ndata]);
    data_l22s = up_data_l22s.get();
    arrutilv2::set_rl22s(ndata, dimension, data, data_l22s);
  }

  void set_distances(size_t i, TFloat *const distances) {
    arrutilv2::set_rl2s(dimension, data + i * dimension, ndata, data,
                        data_l22s[i], data_l22s, distances);
  }
};
}

#endif
