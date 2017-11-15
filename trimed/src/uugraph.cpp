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
#include <cstring>
#include <limits>
#include <numeric>
#include <trimed/uugraph.hpp>
#include <vector>

namespace graph {

std::unique_ptr<size_t[]> UUGraph::get_up_component_start_seach_by_priority() {
  std::unique_ptr<size_t[]> up_component_start_seach_by_priority(
      new size_t[nnodes]);
  auto component_start_seach_by_priority =
      up_component_start_seach_by_priority.get();
  std::iota(component_start_seach_by_priority,
            component_start_seach_by_priority + nnodes, 0);
  return up_component_start_seach_by_priority;
}

void UUGraph::initialise_allocate_edges() {
  up_n_edges.reset(new size_t[nnodes]);
  n_edges = up_n_edges.get();
  up_first_edge.reset(new size_t[nnodes + 1]);
  first_edge = up_first_edge.get();
  up_edges.reset(new size_t[2 * nedges]);
  edges = up_edges.get();
}

void UUGraph::set_edge_arrays(const size_t *const edges_in) {
  size_t internal_a;
  size_t internal_b;

  // (0) set up the number of edges each node contains.
  std::fill_n(n_edges, nnodes, 0);
  for (size_t e = 0; e < nedges; ++e) {
    internal_a = map_to_internalindex[edges_in[2 * e]];
    internal_b = map_to_internalindex[edges_in[2 * e + 1]];
    ++n_edges[internal_a];
    ++n_edges[internal_b];
  }

  // (1) set up the cumulative number of edges each node contains, starting with
  // 0.
  first_edge[0] = 0;
  for (size_t n = 0; n < nnodes; ++n) {
    first_edge[n + 1] = first_edge[n] + n_edges[n];
  }

  // (2) fill up the edges.
  std::unique_ptr<size_t[]> up_current_place(new size_t[nnodes]);
  size_t *current_place = up_current_place.get();
  std::memcpy(current_place, first_edge, sizeof(size_t) * nnodes);
  for (size_t e = 0; e < nedges; ++e) {
    internal_a = map_to_internalindex[edges_in[2 * e]];
    internal_b = map_to_internalindex[edges_in[2 * e + 1]];
    edges[current_place[internal_a]] = internal_b;
    ++current_place[internal_a];
    edges[current_place[internal_b]] = internal_a;
    ++current_place[internal_b];
  }
}

void UUGraph::test_input(size_t nedges_in, const size_t *const edges_in) {
  std::map<std::pair<size_t, size_t>, bool> seen;
  size_t a;
  size_t b;
  std::pair<size_t, size_t> ab;
  for (size_t i = 0; i < nedges_in; ++i) {
    a = edges_in[i * 2];
    b = edges_in[i * 2 + 1];
    if (a == b) {
      throw std::runtime_error("self-edges not acceptable");
    }
    if (a < b) {
      std::swap(a, b);
    }

    ab.first = a;
    ab.second = b;

    if (seen.count(ab) != 0) {
      throw std::runtime_error("no repeating of edges in input allowed");
    }
    seen[ab] = true;
  }
}

UUGraph::UUGraph(size_t nedges_in, const size_t *const edges_in)
    : BaseGraph(nedges_in, edges_in) {

  std::cout << "In UUGraph constructor" << std::endl;
#ifdef DEBUGMODE
  test_input(nedges_in, edges_in);
#endif
  initialiser(nedges_in, edges_in);
  std::cout << "UUGraph construction complete." << std::endl;

  // print();
}

void UUGraph::set_distances(size_t internalindex) {
  set_unweighted_distances(internalindex);
}

/* --------------------------------------------------------------------------------------------------------------------
 */

template class GraphComponent<UUGraph, double>;
}
