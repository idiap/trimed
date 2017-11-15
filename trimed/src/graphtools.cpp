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

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <trimed/graphtools.hpp>
#include <vector>

namespace graphtools {


std::vector<size_t> getsensornet(size_t ndata, double phi_r) {

  size_t dimension = 2;
  std::vector<double> samples(ndata * dimension);
  std::vector<double> firstdimvals(ndata);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);
  for (size_t i = 0; i < ndata; ++i) {
    firstdimvals[i] = dis(gen);
    for (size_t j = 1; j < dimension; ++j) {
      samples[i * dimension + j] = dis(gen);
    }
  }
  std::sort(firstdimvals.begin(), firstdimvals.end());
  for (size_t i = 0; i < ndata; ++i) {
    samples[i * dimension] = firstdimvals[i];
  }

  double radius = phi_r / std::pow(static_cast<double>(ndata),
                                   1. / static_cast<double>(dimension));

  std::cout << "radius : " << radius << std::endl;
  std::vector<size_t> edges;

  double r2;
  double diff;
  size_t k;
  double diff0;
  for (size_t j = 0; j < ndata; ++j) {
    k = j + 1;
    diff0 = 0;
    while (k < ndata && diff0 * diff0 < radius * radius) {
      diff0 = samples[j * dimension + 0] - samples[k * dimension + 0];

      r2 = diff0 * diff0;
      for (size_t l = 1; l < dimension; ++l) {
        diff = (samples[j * dimension + l] - samples[k * dimension + l]);
        r2 += diff * diff;
      }
      if (r2 < radius * radius) {
        edges.push_back(j);
        edges.push_back(k);
      }
      ++k;
    }
  }

  return edges;
}

/* given
 * (1) nodes : ndata x 2  (longitude, latitude) and
 * (2) edges ndata x 2 x 2,
 * return edges in terms of
 * (1) node IDs : ndata x 2 of type size_t, and
 * (2) the lengths of the edges, using longitudinal correction (currently bad)
 */
std::pair<std::vector<size_t>, std::vector<double>>
convert_map(size_t nnodes, size_t nedges, const double *const nodes,
            const double *const edges) {

  std::vector<size_t> ID_edges(nedges * 2);
  std::vector<double> edge_lengths(nedges);

  double pi = atan(1) * 4;

  std::map<std::pair<double, double>, size_t> f_coord;
  for (size_t i = 0; i < nnodes; ++i) {
    f_coord[std::pair<double, double>(nodes[2 * i], nodes[2 * i + 1])] = i;
  }
  double delta_longitude;
  double delta_latitude;
  double mean_latitude;
  double delta_y2;
  double latcor;
  double delta_x2;
  double r2;

  for (size_t i = 0; i < nedges; ++i) {
    ID_edges[2 * i + 0] = f_coord[{edges[4 * i + 0], edges[4 * i + 1]}];
    ID_edges[2 * i + 1] = f_coord[{edges[4 * i + 2], edges[4 * i + 3]}];
    delta_longitude = edges[4 * i] - edges[4 * i + 2];
    delta_latitude = edges[4 * i + 1] - edges[4 * i + 3];
    mean_latitude = (edges[4 * i + 1] + edges[4 * i + 3]) / 2.;
    latcor = std::cos(mean_latitude / 90. * pi / 2.);
    delta_y2 = delta_longitude * latcor * delta_longitude * latcor;
    delta_x2 = delta_latitude * delta_latitude;
    r2 = delta_x2 + delta_y2;
    edge_lengths[i] = std::sqrt(r2);
  }

  return {ID_edges, edge_lengths};
}
}
