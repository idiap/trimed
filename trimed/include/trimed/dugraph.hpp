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

#ifndef DUGRAPH_HPP
#define DUGRAPH_HPP

#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <trimed/basegraph.hpp>

namespace graph {

/* directed-unweighted graph, ie a graph which is directed with unweighted edges
 */
class DUGraph : public BaseGraph {

public:
  virtual void initialise_allocate_edges() override final;
  virtual void set_edge_arrays(const size_t *const edges_in) override final;
  virtual void test_input(size_t nedges_in,
                          const size_t *const edges_in) final override;
  DUGraph(size_t nedges_in, const size_t *const edges_in);
  virtual void set_distances(size_t internalindex) final override;
  std::unique_ptr<size_t[]>
  get_up_component_start_seach_by_priority() final override;

  // the number of reverse edges emanating from each node.
  size_t *n_rev_edges;
  // map to first reverse edge in edges emanating from internalindex, cumulative
  // version of above.
  size_t *first_rev_edge;
  // the reverse  edges, contiguous by node.
  size_t *rev_edges;

  size_t *incoming_distances;

private:
  std::unique_ptr<size_t[]> up_n_rev_edges;
  std::unique_ptr<size_t[]> up_first_rev_edge;
  std::unique_ptr<size_t[]> up_rev_edges;
  std::unique_ptr<size_t[]> up_incoming_distances;
  std::unique_ptr<size_t[]> get_sink_priority();
  void printrev();
};
}

#endif
