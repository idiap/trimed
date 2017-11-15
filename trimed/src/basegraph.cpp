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
#include <vector>
#include <trimed/basegraph.hpp>

namespace graph {

void set_components(size_t component, size_t node_i,
                    const size_t *const first_edge, const size_t *const edges,
                    size_t *const components, size_t *active_now,
                    size_t *active_next) {

  active_now[0] = node_i;
  size_t n_active_now = 1;
  components[node_i] = component;

  size_t n_active_next = 0;

  size_t start;
  size_t end;

  while (n_active_now > 0) {
    for (size_t i = 0; i < n_active_now; ++i) {
      start = active_now[i];
      for (size_t j = first_edge[start]; j < first_edge[start + 1]; ++j) {
        end = edges[j];

        if (components[end] == std::numeric_limits<size_t>::max()) {
          // if (components[end] != component){

          //#ifdef DEBUGMODE
          // if (components[end] != std::numeric_limits<size_t>::max()){
          // throw std::logic_error("Strange encroachment into another
          // component");
          //}
          //#endif

          components[end] = component;
          active_next[n_active_next] = end;
          ++n_active_next;
        }
      }
    }
    std::swap(active_now, active_next);
    n_active_now = n_active_next;
    n_active_next = 0;
  }
}

void set_all_components(size_t nnodes, const size_t *const first_edge,
                        const size_t *const edges, size_t *const components,
                        size_t &ncomponents,
                        const size_t *const component_start_seach_by_priority) {

  std::fill_n(components, nnodes, std::numeric_limits<size_t>::max());

  std::unique_ptr<size_t[]> up_active_now(new size_t[nnodes]);
  std::unique_ptr<size_t[]> up_active_next(new size_t[nnodes]);

  size_t *active_now = up_active_now.get();
  size_t *active_next = up_active_next.get();

  size_t component = 0;

  size_t node_i;

  for (size_t i = 0; i < nnodes; ++i) {
    node_i = component_start_seach_by_priority[i];
    if (components[node_i] == std::numeric_limits<size_t>::max()) {
      set_components(component, node_i, first_edge, edges, components,
                     active_now, active_next);
      ++component;
    }
  }
  ncomponents = component;
}

void BaseGraph::set_components() {
  /* (1) set components of all nodes */
  std::unique_ptr<size_t[]> up_component_start_seach_by_priority =
      get_up_component_start_seach_by_priority();
  set_all_components(nnodes, first_edge, edges, components, ncomponents,
                     up_component_start_seach_by_priority.get());

  /* (2) some initialisation */
  up_component_counts.reset(new size_t[ncomponents]);
  component_counts = up_component_counts.get();
  up_component_starts.reset(new size_t[ncomponents]);
  component_starts = up_component_starts.get();

  /* (3) set counts of components */
  std::fill_n(component_counts, ncomponents, 0);
  std::fill_n(component_starts, ncomponents,
              std::numeric_limits<size_t>::max());
  for (size_t n = 0; n < nnodes; ++n) {
    ++component_counts[components[n]];
    if (component_starts[components[n]] == std::numeric_limits<size_t>::max()) {
      component_starts[components[n]] = n;
    }
  }
}

void BaseGraph::print() {

  std::cout << "nnodes : " << nnodes << std::endl;
  std::cout << "nedges : " << nedges << std::endl;
  std::cout << "ncomponents : " << ncomponents << std::endl;

  std::cout << "\nmap_to_inputindex : " << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << n << " : " << map_to_inputindex[n] << std::endl;
  }

  std::cout << "\nmap_to_internalindex : " << std::endl;
  for (auto &x : map_to_internalindex) {
    std::cout << x.first << " : " << x.second << std::endl;
  }

  std::cout << "\nn_edges and first_edge" << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << n << " : " << n_edges[n] << " : " << first_edge[n]
              << std::endl;
  }

  std::cout << "\nedges " << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << "( " << n << " ) : ";
    for (size_t i = first_edge[n]; i < first_edge[n + 1]; ++i) {
      std::cout << " " << edges[i] << std::flush;
    }
    std::cout << std::endl;
  }

  std::cout << "\ncomponents " << std::endl;
  for (size_t n = 0; n < nnodes; ++n) {
    std::cout << n << " : " << components[n] << std::endl;
  }

  std::cout << "\ncomponent_counts & component_starts" << std::endl;
  for (size_t c = 0; c < ncomponents; ++c) {
    std::cout << c << " : " << component_counts[c] << " : "
              << component_starts[c] << std::endl;
  }
}

void BaseGraph::link_by_input_ID_order() {
  size_t internalindex = 0;
  for (auto &x : map_to_internalindex) {
    x.second = internalindex;
    map_to_inputindex[internalindex] = x.first;
    ++internalindex;
  }
}

void BaseGraph::link_so_contiguous_components() {

  std::map<size_t, size_t> inputindex_to_component;
  for (size_t n = 0; n < nnodes; ++n) {
    inputindex_to_component[map_to_inputindex[n]] = components[n];
  }

  std::unique_ptr<size_t[]> up_cum_component_counts(new size_t[ncomponents]);
  up_cum_component_counts[0] = 0;
  for (size_t c = 1; c < ncomponents; ++c) {
    up_cum_component_counts[c] =
        up_cum_component_counts[c - 1] + component_counts[c - 1];
  }

  std::vector<size_t> current_position(ncomponents);
  std::memcpy(current_position.data(), up_cum_component_counts.get(),
              sizeof(size_t) * ncomponents);

  for (size_t n = 0; n < nnodes; ++n) {
    map_to_internalindex[map_to_inputindex[n]] =
        current_position[components[n]];
    ++current_position[components[n]];
  }

  for (auto &x : map_to_internalindex) {
    map_to_inputindex[x.second] = x.first;
  }
}

// std::vector<size_t> new_map_to_inputindex (nnodes);
//
// for (size_t n = 0; n < nnodes; ++n){
// new_map_to_inputindex[n] = current_position[components[n]];
//++current_position[components[n]];
//}

// size_t acounter = 0;
// for (auto & x : map_to_internalindex){
// x.second = new_map_to_inputindex[acounter];
// map_to_inputindex[new_map_to_inputindex[acounter]] = x.first;
//++acounter;
//}

void BaseGraph::initialiser(size_t nedges_in, const size_t *const edges_in) {

  nedges = nedges_in;

  /* establish the input indices */
  for (size_t ei = 0; ei < 2 * nedges; ++ei) {
    map_to_internalindex[edges_in[ei]] = std::numeric_limits<size_t>::max();
  }
  /* number of distinct nodes */
  nnodes = map_to_internalindex.size();

  /* some initialisation  */
  up_map_to_inputindex.reset(new size_t[nnodes]);
  map_to_inputindex = up_map_to_inputindex.get();
  up_components.reset(new size_t[nnodes]);
  components = up_components.get();
  up_distances.reset(new size_t[nnodes]);
  distances = up_distances.get();

  /* initial linking of input and internal indices */
  link_by_input_ID_order();

  /* initial setting of edges and components (virtuals) */
  initialise_allocate_edges();
  set_edge_arrays(edges_in);

  /* set the components (non-virtual) */
  set_components();

  /* now that the components have been established, we can reset internal
   * indices so that they are contiguous by component. This will allow for
   * improved memory hits in traversal */
  link_so_contiguous_components();

  /* reset edges and components with the contiguous interal indexing */
  set_edge_arrays(edges_in);
  set_components();

  // std::cout << "\ncomponent_counts & component_starts" << std::endl;
  // for (size_t c = 0; c < 20; ++c){
  // std::cout << c << " : " << component_counts[c] << " : " <<
  // component_starts[c] << std::endl;
  //}
}

void BaseGraph::set_unweighted_distances(size_t internalindex) {
  set_unweighted_distances(internalindex, distances, first_edge, edges);
}

// generic (can be used by forward and round-trip)
void BaseGraph::set_unweighted_distances(size_t internalindex,
                                         size_t *const x_distances,
                                         const size_t *const x_first_edge,
                                         const size_t *const x_edges) {

  size_t component = components[internalindex];

  std::unique_ptr<size_t[]> up_minimal_X(
      new size_t[component_counts[component]]);
  std::unique_ptr<size_t[]> up_minimal_Y(
      new size_t[component_counts[component]]);
  size_t *minimals = up_minimal_X.get();
  size_t *minimals_p1 = up_minimal_Y.get();

  std::fill_n(x_distances + component_starts[component],
              component_counts[component], std::numeric_limits<size_t>::max());
  x_distances[internalindex] = 0;

  minimals[0] = internalindex;
  size_t n_minimal = 1;
  size_t n_minimal_p1 = 0;
  size_t front_distance = 1;

  // To help make code clearer we use these, although not nec. Hopefully
  // compiler will make them vanish (check). If not come and efficienticate.
  size_t source;
  size_t target;

  while (n_minimal > 0) {
    for (size_t i = 0; i < n_minimal; ++i) {
      source = minimals[i];
      for (size_t j = x_first_edge[source]; j < x_first_edge[source + 1]; ++j) {
        target = x_edges[j];
        if (x_distances[target] > front_distance) {
          x_distances[target] = front_distance;
          minimals_p1[n_minimal_p1] = target;
          ++n_minimal_p1;
        }
      }
    }
    ++front_distance;
    std::swap(minimals, minimals_p1);
    n_minimal = n_minimal_p1;
    n_minimal_p1 = 0;
  }
}

BaseGraph::BaseGraph(size_t nedges_in, const size_t *const edges_in) {
  (void)nedges_in;
  (void)edges_in;
  // TODO : some type of base test
}
}
