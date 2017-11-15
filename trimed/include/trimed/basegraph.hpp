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

#ifndef BASEGRAPH_H
#define BASEGRAPH_H

#include <iostream>
#include <map>
#include <memory>

#define DEBUGMODE

namespace graph {

void set_components(size_t component, size_t node_i,
                    const size_t *const first_edge, const size_t *const edges,
                    size_t *const components, size_t *active_now,
                    size_t *active_next);

void set_all_components(size_t nnodes, const size_t *const first_edge,
                        const size_t *const edges, size_t *const components,
                        size_t &ncomponents);

/* base class for {undirected, directed} x {unweighted, weighted} graphs */
class BaseGraph {
public:
  BaseGraph(size_t nedges_in, const size_t *const edges_in);

  virtual ~BaseGraph() {}

  // maps inputindex to internalindex (nodes).
  std::map<size_t, size_t> map_to_internalindex;

  // maps internalindex to inputindex.
  size_t *map_to_inputindex;

  // total number of nodes.
  size_t nnodes;

  // total number of edges.
  size_t nedges;

  // total number of components.
  size_t ncomponents;

  // the number of edges emanating from each node.
  size_t *n_edges;

  // map to first edge in edges emanating from internalindex, cumulative version
  // of above.
  size_t *first_edge;

  // the edges, contiguous by node.
  size_t *edges;

  // the component to which nodes belong.
  size_t *components;

  // the size of each component.
  size_t *component_counts;
  size_t *component_starts;

  size_t *distances;
  virtual void set_distances(size_t internalindex) = 0;

protected:
  // undirected versions need 2*edges_in while directed need 2*edges_in
  virtual void initialise_allocate_edges() = 0;
  virtual void set_edge_arrays(const size_t *const edges_in) = 0;
  virtual void test_input(size_t nedges_in, const size_t *const edges_in) = 0;
  virtual std::unique_ptr<size_t[]>
  get_up_component_start_seach_by_priority() = 0;

  void set_components();
  void link_by_input_ID_order();
  void link_so_contiguous_components();
  void initialiser(size_t nedges_in, const size_t *const edges_in);
  void print();
  void set_unweighted_distances(size_t internalindex);
  void set_unweighted_distances(size_t internalindex, size_t *const x_distances,
                                const size_t *const x_first_edge,
                                const size_t *const x_edges);

  std::unique_ptr<size_t[]> up_map_to_inputindex;
  std::unique_ptr<size_t[]> up_n_edges;
  std::unique_ptr<size_t[]> up_first_edge;
  std::unique_ptr<size_t[]> up_edges;
  std::unique_ptr<size_t[]> up_component_counts;
  std::unique_ptr<size_t[]> up_component_starts;
  std::unique_ptr<size_t[]> up_components;
  std::unique_ptr<size_t[]> up_distances;
};

template <typename TGraph, typename TFloat> class GraphComponent {
public:
  size_t ndata;

  GraphComponent(TGraph *agraphi, size_t componenti) {
    agraph = agraphi;
    component = componenti;
    ndata = agraph->component_counts[component];
  }

  void set_distances(size_t i, TFloat *const distances) {
    agraph->set_distances(i + agraph->component_starts[component]);
    for (size_t j = 0; j < ndata; ++j) {
      distances[j] = static_cast<TFloat>(
          agraph->distances[j + agraph->component_starts[component]]);
    }
  }

private:
  TGraph *agraph;
  size_t component;
};
}

#endif
