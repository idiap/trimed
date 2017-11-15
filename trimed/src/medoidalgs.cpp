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
#include <iostream>
#include <memory>
#include <trimed/arraydata.hpp>
#include <trimed/medoidalgs.hpp>
#include <trimed/uugraph.hpp>

namespace shuffle {
/* returns range [r0, r1) in random order */
std::unique_ptr<size_t[]> get_up_shuffled_range(size_t r0, size_t r1) {
  std::unique_ptr<size_t[]> up_shuffled_range(new size_t[r1 - r0]);
  auto shuffled_range = up_shuffled_range.get();
  std::iota(shuffled_range, shuffled_range + r1 - r0, r0);
  // shuffle:
  for (size_t i = 1; i < r1 - r0; ++i) {
    std::swap(shuffled_range[rand() % i], shuffled_range[i]);
  }
  return std::move(up_shuffled_range);
}
}

namespace medoidalgs {

template <class Data, typename TFloat>
void RAND(Data &data, TFloat &min_E, size_t &min_E_index, size_t maxcomputes) {

  std::cout << "In algorithm RAND for finding (approximate) medoid"
            << std::endl;
  min_E = std::numeric_limits<TFloat>::max(); // the mimumum energy found so far

  std::unique_ptr<TFloat[]> up_sums(new TFloat[data.ndata]);
  auto sums = up_sums.get();
  std::fill_n(sums, data.ndata, 0);

  std::unique_ptr<TFloat[]> up_distances(new TFloat[data.ndata]);
  auto distances = up_distances.get();

  std::unique_ptr<bool[]> up_computed(new bool[data.ndata]);
  auto computed = up_computed.get();
  std::fill_n(computed, data.ndata, false);

  // shuffle order.
  auto up_visiting_order = shuffle::get_up_shuffled_range(0, data.ndata);

  size_t min_sum_index;   // the index minimising the sum
  TFloat E_min_sum_index; // the energy at of min_sum_i

  size_t ncomputes = 0;

  size_t i;   // the random sample to compute
  TFloat E_i; // the energy at the random sample computed.

  std::cout << "Round\tncomputes\tmin_E\tmin_E_index" << std::endl;
  size_t round = 0;
  while (ncomputes < maxcomputes) {
    i = up_visiting_order[round];
    data.set_distances(i, distances);
    computed[i] = true;
    ++ncomputes;
    arrutilv2::addto(data.ndata, distances, sums);
    // candidate 1 : i
    arrutilv2::set_sum(data.ndata, distances, E_i);
    if (E_i < min_E) {
      min_E = E_i;
      min_E_index = i;
      std::cout << round << "\t" << ncomputes << "\t" << min_E << "\t"
                << min_E_index << std::endl;
    }

    // candidate 2 : min_sum_index (if not already checked)
    arrutilv2::set_argmin(data.ndata, sums, min_sum_index);
    if (computed[min_sum_index] == false) {
      data.set_distances(min_sum_index, distances);
      computed[min_sum_index] = true;
      ++ncomputes;
      arrutilv2::set_sum(data.ndata, distances, E_min_sum_index);
      if (E_min_sum_index < min_E) {
        min_E = E_min_sum_index;
        min_E_index = min_sum_index;
        std::cout << round << "\t" << ncomputes << "\t" << min_E << "\t"
                  << min_E_index << std::endl;
      }
    }
    ++round;
  }
  std::cout << round << "\t" << ncomputes << "\t" << min_E << "\t"
            << min_E_index << std::endl;
  std::cout << "Complete." << std::endl;
}

template <class Data, typename TFloat>
void trimed(Data &data, TFloat &min_E, size_t &min_E_i) {

  std::cout << "In algorithm trimed for finding medoid" << std::endl;

  // initialise min E found so far to be largest possible value
  min_E = std::numeric_limits<TFloat>::max();
  TFloat E_i;

  // initialise the lower bounds on energies for all samples to zero
  std::unique_ptr<TFloat[]> up_lowers(new TFloat[data.ndata]);
  auto lowers = up_lowers.get();
  std::fill_n(lowers, data.ndata, 0);

  // track which samples have had their energies computed
  std::unique_ptr<bool[]> up_computed(new bool[data.ndata]);
  auto computed = up_computed.get();
  std::fill_n(computed, data.ndata, false);

  std::unique_ptr<TFloat[]> up_distances(new TFloat[data.ndata]);
  auto distances = up_distances.get();

  // shuffle order
  auto up_visiting_order = shuffle::get_up_shuffled_range(0, data.ndata);

  size_t i;
  size_t ncomputes = 0;
  size_t nsaves = 0;
  std::cout << "round\tmin_E\tncalcs\tnsaves\tmin_E_i" << std::endl;
  for (size_t io = 0; io < data.ndata; ++io) {
    i = up_visiting_order[io];

    // if the lower bound on sample i is less than the lowest energy found so
    // far,
    if (lowers[i] < min_E) {

      // compute all N distances.
      data.set_distances(i, distances);
      ++ncomputes;
      // compute the sum of all distances (ie energy of sample i).
      arrutilv2::set_sum(data.ndata, distances, E_i);

      computed[i] = true;
      // if the energy of sample i is lower than the lowest found so far, update
      // lowest.
      if (E_i < min_E) {
        min_E = E_i;
        min_E_i = i;
        std::cout << io << "\t" << min_E << "\t" << ncomputes << "\t" << nsaves
                  << "\t" << min_E_i << std::endl;
      }

      // update the lower bounds for all samples using the N distances computed
      // to sample i.
      for (size_t ip = i + 1; ip < data.ndata; ++ip) {
        lowers[ip] = std::max<TFloat>(
            lowers[ip], std::abs(data.ndata * distances[ip] - E_i));
      }

    }

    // the lower bound test saved N distance calculations, record this.
    else {
      ++nsaves;
    }
  }
  std::cout << data.ndata << "\t" << min_E << "\t" << ncomputes << "\t"
            << nsaves << "\t" << min_E_i << std::endl;
}

// The TOPRANK algorithm of Okamoto, alpha' as per correspondence, log base e.
// For k = 1.
template <class Data, typename TFloat>
void TOPRANK1(Data &data, TFloat &min_E, size_t &min_E_i) {

  std::cout << "In algorithm TOPRANK1 for finding (w.h.p) medoid" << std::endl;

  std::unique_ptr<TFloat[]> up_sums(new TFloat[data.ndata]);
  auto sums = up_sums.get();
  std::fill_n(sums, data.ndata, 0);

  std::unique_ptr<TFloat[]> up_distances(new TFloat[data.ndata]);
  auto distances = up_distances.get();

  std::unique_ptr<bool[]> up_computed(new bool[data.ndata]);
  auto computed = up_computed.get();
  std::fill_n(computed, data.ndata, false);

  // shuffle order.
  auto up_visiting_order = shuffle::get_up_shuffled_range(0, data.ndata);

  // std::vector<size_t> up_visiting_order(data.ndata);
  // std::iota(up_visiting_order.begin(), up_visiting_order.end(), 0);

  TFloat ndata_float = static_cast<TFloat>(data.ndata);
  TFloat ndata_2over3_logterm =
      std::pow(ndata_float, 2. / 3.) * std::cbrt(std::log(ndata_float));
  //(0) take l to be min(N, N^(2/3) * N^(1/3)), could be anything for w.hpp>.p
  // result.
  size_t l = std::min(data.ndata, static_cast<size_t>(ndata_2over3_logterm));
  TFloat l_float = static_cast<TFloat>(l);
  std::cout << "ndata : " << data.ndata << std::endl;
  std::cout
      << "computed l (as per Lemma 2, page 191 : n^(2/3) (log(n))^(1/3)) : "
      << l << std::endl;

  TFloat max_distance;
  TFloat delta_hat = std::numeric_limits<TFloat>::max();

  size_t ndcalcs = 0;

  //(1) get sums from l random points, and set delta_hat.
  std::cout << "computing sums for l randomly selected points " << l
            << std::endl;
  for (size_t i = 0; i < l; ++i) {
    data.set_distances(up_visiting_order[i], distances);
    ndcalcs += 1;
    arrutilv2::addto(data.ndata, distances, sums);

    max_distance = 0.0;
    for (size_t j = 0; j < data.ndata; ++j) {
      max_distance = std::max<TFloat>(max_distance, distances[j]);
    }

    delta_hat = std::min<TFloat>(delta_hat, 2. * max_distance);
  }

  //(2) set a_hat_1 to minimum estimated energy (sum[j] / l),
  TFloat a_hat_1 = std::numeric_limits<TFloat>::max();
  for (size_t j = 0; j < data.ndata; ++j) {
    a_hat_1 = std::min(a_hat_1, sums[j] / l_float);
  }

  //(3) set threshold to a_hat_1 + 2*f(l)*delta_hat, where f(l) is alpha'
  // sqrt(log n / l )
  // Correspondence with the authors suggests taking alpha' to be 4 sqrt(2.5) /
  // 3, which is 2 sqrt(10) / 3.
  // However, taking 1 seems to work be fine.
  TFloat alpha_prime = 1.0; //(2.*std::sqrt(10.))/3.;
  TFloat f_of_l = alpha_prime * std::sqrt(std::log(ndata_float) / l_float);
  TFloat threshold = a_hat_1 + 2 * delta_hat * f_of_l;

  std::cout << "\na_hat_1 : " << a_hat_1 << "\t delta_hat : " << delta_hat
            << "\t f_of_l : " << f_of_l << "\t threshold : " << threshold
            << std::endl;
  //(4) get the minimum below threshold.
  TFloat E_j;
  min_E = std::numeric_limits<TFloat>::max();
  for (size_t j = 0; j < data.ndata; ++j) {
    if (sums[j] / l_float < threshold) {
      data.set_distances(j, distances);
      ++ndcalcs;
      arrutilv2::set_sum(data.ndata, distances, E_j);
      if (E_j < min_E) {
        min_E = E_j;
        min_E_i = j;
      }
    }
  }
  std::cout << "l\tncomputes\tmin_E\tmin_E_i" << std::endl;
  std::cout << l << "\t" << ndcalcs << "\t" << min_E << "\t" << min_E_i
            << std::endl;
}

// The TOPRANK2 algorithm of Okamoto, log base e. For k = 1.
template <class Data, typename TFloat>
void TOPRANK2(Data &data, TFloat &min_E, size_t &min_E_i) {

  std::cout << "In algorithm TOPRANK2 for finding approximate medoid"
            << std::endl;

  std::unique_ptr<TFloat[]> up_sums(new TFloat[data.ndata]);
  auto sums = up_sums.get();
  std::fill_n(sums, data.ndata, 0);

  std::unique_ptr<TFloat[]> up_distances(new TFloat[data.ndata]);
  auto distances = up_distances.get();

  std::unique_ptr<bool[]> up_computed(new bool[data.ndata]);
  auto computed = up_computed.get();
  std::fill_n(computed, data.ndata, false);

  // shuffle order.
  auto up_visiting_order = shuffle::get_up_shuffled_range(0, data.ndata);

  size_t ndcalcs = 0;

  TFloat ndata_float = static_cast<TFloat>(data.ndata);

  // number of samples used to obtain Eppstein sum. Grows by q, so l1 = l0 + q
  // after initial round.
  size_t l0 = 0;
  size_t l1 =
      static_cast<size_t>(std::sqrt(ndata_float)); // choosing this is an art...
  std::cout << "l1 set to " << l1 << " which is sqrt(n) where n is "
            << ndata_float << "  (" << data.ndata << ")" << std::endl;
  // as recommended in Okamoto et al.
  size_t q = static_cast<size_t>(std::log(ndata_float));

  // given l0, how many samples are below threshold and would need to be
  // computed?
  size_t p0 = 0;
  size_t p1 = 0;

  // current round.
  size_t round = 0;

  TFloat alpha_prime = 0;
  TFloat f_of_l = 0;
  TFloat threshold = 0;

  TFloat max_distance;
  TFloat delta_hat = std::numeric_limits<TFloat>::max();

  std::cout << "l [p0 - p1] = " << std::flush;

  // to determine a good l1, we roll until increasing l1 is more costly than the
  // subsquent sub-threshold computation would be.
  while (l1 < data.ndata && (round < 3 || (int(q) < int(p0) - int(p1)))) {
    if (round >= 2) {
      std::cout << l1 << "[" << int(p0) - int(p1) << "] ... " << std::flush;
    }
    ++round;
    for (size_t i = l0; i < l1; ++i) {
      ndcalcs += 1;
      data.set_distances(up_visiting_order[i], distances);
      max_distance = 0;
      for (size_t j = 0; j < data.ndata; ++j) {
        sums[j] += distances[j];
        max_distance = std::max(max_distance, distances[j]);
      }
      delta_hat = std::min(delta_hat, 2 * max_distance);
    }

    // set a_hat_1 to minimum estimated energy,
    TFloat a_hat_1 = std::numeric_limits<TFloat>::max();
    for (size_t j = 0; j < data.ndata; ++j) {
      a_hat_1 = std::min(a_hat_1, sums[j] / static_cast<TFloat>(l1));
    }

    // set threshold to a_hat_1 + 2*f(l)*delta_hat, as per TOPRANK1.
    alpha_prime = 1; //(2.*std::sqrt(10.))/3.;
    f_of_l = alpha_prime * std::sqrt(std::log(static_cast<TFloat>(data.ndata)) /
                                     static_cast<TFloat>(l1));
    threshold = a_hat_1 + 2 * delta_hat * f_of_l;

    std::cout << "a_hat_1 :  " << a_hat_1 << "\tdelta_hat : " << delta_hat
              << "\tf_of_l :  " << f_of_l << "\tthreshold : " << threshold
              << std::endl;
    // count the number below threshold.
    p0 = p1;
    p1 = 0;
    for (size_t j = 0; j < data.ndata; ++j) {
      if (sums[j] / static_cast<TFloat>(l1) < threshold) {
        ++p1;
      }
    }
    l0 = l1;
    l1 += q;
  }

  std::cout << l1 << "[" << int(p0) - int(p1) << "]" << std::endl;

  // compute true distances at all sub-threshold values.
  min_E = std::numeric_limits<TFloat>::max();
  TFloat E_j;
  for (size_t j = 0; j < data.ndata; ++j) {
    if (sums[j] / static_cast<TFloat>(l1) < threshold) {
      data.set_distances(j, distances);
      ++ndcalcs;
      E_j = 0;
      for (size_t jp = 0; jp < data.ndata; ++jp) {
        E_j += distances[jp];
      }
      if (E_j < min_E) {
        min_E = E_j;
        min_E_i = j;
      }
    }
  }
  std::cout << "l_final\tncomps\tmin_E\tmin_E_i" << std::endl;
  std::cout << l1 << "\t" << ndcalcs << "\t" << min_E << "\t" << min_E_i << "\t"
            << std::endl;
}

template <class Data, typename TFloat>
void medoid(Data &data, TFloat &min_E, size_t &min_E_i, std::string algorithm,
            bool capture_output, std::string &text, size_t seed,
            size_t maxcomputes) {

  srand(seed);

  text = "";
  std::stringstream buffer;
  auto cout_buff = std::cout.rdbuf();
  if (capture_output == true) {
    std::cout.rdbuf(buffer.rdbuf());
  }

  if (algorithm.compare("trimed") == 0) {
    trimed(data, min_E, min_E_i);
  }

  else if (algorithm.compare("RAND") == 0) {
    RAND(data, min_E, min_E_i, maxcomputes);
  }

  else if (algorithm.compare("TOPRANK1") == 0) {
    TOPRANK1(data, min_E, min_E_i);
  }

  else if (algorithm.compare("TOPRANK2") == 0) {
    TOPRANK2(data, min_E, min_E_i);
  }

  else {
    throw std::runtime_error("Unrecognised algorithm, " + algorithm + ".");
  }

  if (capture_output == true) {
    text = buffer.str();
    std::cout.rdbuf(cout_buff);
  }

  else {
    text = "capture_output was false, output seen already";
  }
}

/* get medoid of largest connected component on an unsigned graph */
template <typename TGraph>
void medoid_graph(size_t nedges_in, const size_t *const edges_in, double &min_E,
                  size_t &min_E_i, std::string algorithm, bool capture_output,
                  std::string &text, size_t seed, size_t maxcomputes) {

  TGraph agraph(nedges_in, edges_in);
  /* get max component */
  size_t max_component = 0;
  size_t max_component_size = 0;
  for (size_t c = 0; c < agraph.ncomponents; ++c) {
    if (agraph.component_counts[c] > max_component_size) {
      max_component_size = agraph.component_counts[c];
      max_component = c;
    }
  }

  graph::GraphComponent<TGraph, double> acomponent(&agraph, max_component);

  medoid(acomponent, min_E, min_E_i, algorithm, capture_output, text, seed,
         maxcomputes);
  /* map index to original index */
  min_E_i =
      agraph
          .map_to_inputindex[min_E_i + agraph.component_starts[max_component]];
}

void medoid_uugraph(size_t nedges_in, const size_t *const edges_in,
                    double &min_E, size_t &min_E_i, std::string algorithm,
                    bool capture_output, std::string &text, size_t seed,
                    size_t maxcomputes) {
  medoid_graph<graph::UUGraph>(nedges_in, edges_in, min_E, min_E_i, algorithm,
                               capture_output, text, seed, maxcomputes);
}

void medoid_dugraph(size_t nedges_in, const size_t *const edges_in,
                    double &min_E, size_t &min_E_i, std::string algorithm,
                    bool capture_output, std::string &text, size_t seed,
                    size_t maxcomputes) {
  medoid_graph<graph::DUGraph>(nedges_in, edges_in, min_E, min_E_i, algorithm,
                               capture_output, text, seed, maxcomputes);
}

/* get medoid of array of double data */
template void medoid_array(size_t ndata, size_t dimension,
                           const double *const data, double &min_E,
                           size_t &min_E_i, std::string algorithm,
                           bool capture_output, std::string &text, size_t seed,
                           size_t maxcomputes);
}
