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
#include <trimed/dugraph.hpp>
#include <vector>

namespace graph {

std::unique_ptr<size_t[]> DUGraph::get_up_component_start_seach_by_priority() {
  return get_sink_priority();
}

// this order needs to be sink_priority for a directed graph, as per Python
// code.
std::unique_ptr<size_t[]> DUGraph::get_sink_priority() {

  std::unique_ptr<bool[]> up_has_had_visit(new bool[nnodes]);
  auto has_had_visit = up_has_had_visit.get();
  std::fill_n(has_had_visit, nnodes, 0);

  std::unique_ptr<size_t[]> up_post_visit_number(new size_t[nnodes]);
  auto post_visit_number = up_post_visit_number.get();
  size_t n_top;
  size_t child_i = 0;
  size_t child_n;
  bool unvisited_child_found;
  std::vector<size_t> stack;

  size_t n_post_visits = 0;
  for (size_t i = 0; i < nnodes; ++i) {
    if (has_had_visit[i] == false) {
      stack.push_back(i);
      has_had_visit[i] = true;
      while (stack.size() > 0) {
        n_top = stack.back();

        child_i = first_rev_edge[n_top];
        unvisited_child_found = false;
        while (unvisited_child_found == false &&
               child_i < first_rev_edge[n_top + 1]) {
          child_n = rev_edges[child_i];
          if (has_had_visit[child_n] == false) {
            unvisited_child_found = true;
            has_had_visit[child_n] = true;
            stack.push_back(child_n);
          } else {
            ++child_i;
          }
        }

        if (unvisited_child_found == false) {
          post_visit_number[n_top] = n_post_visits;
          ++n_post_visits;
          stack.pop_back();
        }
      }
    }
  }

  std::unique_ptr<size_t[]> up_sink_priority(new size_t[nnodes]);
  auto sink_priority = up_sink_priority.get();
  for (size_t i = 0; i < nnodes; ++i) {
    sink_priority[nnodes - 1 - post_visit_number[i]] = i;
  }

  return up_sink_priority;
}

void DUGraph::initialise_allocate_edges() {
  up_n_edges.reset(new size_t[nnodes]);
  n_edges = up_n_edges.get();
  up_first_edge.reset(new size_t[nnodes + 1]);
  first_edge = up_first_edge.get();
  up_edges.reset(new size_t[1 * nedges]);
  edges = up_edges.get();

  up_n_rev_edges.reset(new size_t[nnodes]);
  n_rev_edges = up_n_rev_edges.get();
  up_first_rev_edge.reset(new size_t[nnodes + 1]);
  first_rev_edge = up_first_rev_edge.get();
  up_rev_edges.reset(new size_t[1 * nedges]);
  rev_edges = up_rev_edges.get();
}

void DUGraph::set_edge_arrays(const size_t *const edges_in) {
  size_t internal_a;
  size_t internal_b;

  // (0) set up the number of edges and rev-edges for each node.
  std::fill_n(n_edges, nnodes, 0);
  std::fill_n(n_rev_edges, nnodes, 0);

  for (size_t e = 0; e < nedges; ++e) {
    internal_a = map_to_internalindex[edges_in[2 * e]];
    internal_b = map_to_internalindex[edges_in[2 * e + 1]];
    ++n_edges[internal_a];
    ++n_rev_edges[internal_b];
  }

  // (1) set up the cumulative number of edges and rev edges of each node,
  // starting with 0.
  first_edge[0] = 0;
  first_rev_edge[0] = 0;
  for (size_t n = 0; n < nnodes; ++n) {
    first_edge[n + 1] = first_edge[n] + n_edges[n];

    first_rev_edge[n + 1] = first_rev_edge[n] + n_rev_edges[n];
  }

  // (2) fill up the edges and rev edges.
  std::unique_ptr<size_t[]> up_current_place(new size_t[nnodes]);
  size_t *current_place = up_current_place.get();
  std::memcpy(current_place, first_edge, sizeof(size_t) * nnodes);

  std::unique_ptr<size_t[]> up_current_rev_place(new size_t[nnodes]);
  size_t *current_rev_place = up_current_rev_place.get();
  std::memcpy(current_rev_place, first_rev_edge, sizeof(size_t) * nnodes);

  for (size_t e = 0; e < nedges; ++e) {
    internal_a = map_to_internalindex[edges_in[2 * e]];
    internal_b = map_to_internalindex[edges_in[2 * e + 1]];
    edges[current_place[internal_a]] = internal_b;
    ++current_place[internal_a];
    rev_edges[current_rev_place[internal_b]] = internal_a;
    ++current_rev_place[internal_b];
  }
}

void DUGraph::test_input(size_t nedges_in, const size_t *const edges_in) {

  std::map<std::pair<size_t, size_t>, bool> seen;
  size_t a;
  size_t b;
  std::pair<size_t, size_t> ab;
  for (size_t i = 0; i < nedges_in; ++i) {
    a = edges_in[i * 2];
    b = edges_in[i * 2 + 1];
    ab.first = a;
    ab.second = b;

    if (seen.count(ab) != 0) {
      throw std::runtime_error("no repeating of edges in input allowed");
    }
    seen[ab] = true;
  }
}

DUGraph::DUGraph(size_t nedges_in, const size_t *const edges_in)
    : BaseGraph(nedges_in, edges_in) {

#ifdef DEBUGMODE
  test_input(nedges_in, edges_in);
#endif
  initialiser(nedges_in, edges_in);
  up_incoming_distances.reset(new size_t[nnodes]);
  incoming_distances = up_incoming_distances.get();

  // print();
  // printrev();
}

void DUGraph::printrev() {

  std::cout << "\nn_rev_edges and first_rev_edge" << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << n << " : " << n_rev_edges[n] << " : " << first_rev_edge[n]
              << std::endl;
  }

  std::cout << "\nrev_edges " << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << "( " << n << " ) : ";
    for (size_t i = first_rev_edge[n]; i < first_rev_edge[n + 1]; ++i) {
      std::cout << " " << rev_edges[i] << std::flush;
    }
    std::cout << std::endl;
  }
}

void DUGraph::set_distances(size_t internalindex) {

  // this->distances is set to outgoing distances.
  set_unweighted_distances(internalindex, distances, first_edge, edges);

  // set this->incoming_distances to incoming distances,
  set_unweighted_distances(internalindex, incoming_distances, first_rev_edge,
                           rev_edges); // TODO : problem solve : (1)

  // add incoming distances to distances to, to make distances the round-trip
  // distances,
  size_t component = components[internalindex];

  // std::cout << "\n" << "round trip from " << internalindex << " (component "
  // << component << ")" << std::endl;
  for (size_t i = component_starts[component];
       i < component_starts[component] + component_counts[component]; ++i) {
    distances[i] += incoming_distances[i];
    // std::cout << " [ " << i << " : " << distances[i] << " , " <<
    // incoming_distances[i] << " ] " << std::flush;
  }
}

/* -------------------------------------------*/

template class GraphComponent<DUGraph, double>;
}
