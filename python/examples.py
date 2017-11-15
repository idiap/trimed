#######################################################################
# trimed, a C++ library and (optional) Python tool for obtaining the
# medoid of a set
# 
# Copyright (c) 2017 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# 
# This file is part of trimed.
# 
# trimed is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
# 
# trimed is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Foobar. If not, see <http://www.gnu.org/licenses/>.
#######################################################################

"""
Some examples using pytrimed.
"""

import sys
import random
import numpy as np
import numpy.random as npr

#where is pytrimed.so ? Make sure this is correct.
sys.path.append("../build/python")
import pytrimed

from IPython.core.debugger import Tracer 

def get_n_computed(X, alg):
  if alg == "trimed":
    return int(X.split("\n")[-2].split()[2])
  else:
    return int(X.split("\n")[-2].split()[1])


def basic_experiment():
  ndata = 10**4
  good = 0
  nruns = 1
  alg = "trimed"
  d = 2
  data = npr.rand(ndata,d).astype(np.float64)
  print "d = ", d, "  alg = ", alg
  X = pytrimed.pmedoid(data, datatype = "points", algorithm = alg, capture_output = False, maxcomputes = None)
  print "\n"


def experiment2():  
  """
  prototype of curve in Figure 2 (left). Extend algs, dims amd ndatas lists/arrays to get full figure
  """

  import matplotlib.pyplot as pl
  pl.ion()
  
  c_outp = True
  n_calcs = {}
  
  algs = ["trimed"]#, "TOPRANK1"]:
  dims = [2,3,5]
  ndatas = np.array([int(x) for x in 10**np.linspace(2, 4.6, 7)])
  
  for alg in algs:
    n_calcs[alg] = {}
    for d in dims:
      n_calcs[alg][d] = []
      for ndata in ndatas:
        print "\n\nd = ", d, "  alg = ", alg, "  ndata = ", ndata
        data = npr.rand(ndata,d).astype(np.float64)
        X = pytrimed.pmedoid(data, datatype = "points", algorithm = alg, capture_output = c_outp, maxcomputes = None)
        n_calcs[alg][d].append(get_n_computed(X["text"], alg))


  
  pl.clf()
  for alg in algs:
    for d in dims:
      pl.plot(ndatas, n_calcs[alg][d], label = "alg=%s d=%s"%(alg, d), linestyle=":", marker = "o")
  pl.legend(loc = "upper left")
  pl.show()
  
  
  
def experiment_U_Sensor_net():
  """
  As per Table 1 of paper, column trimed (change algorithm argument for different algorithm
  """
  ndata = 3.6*10**5
  edges = pytrimed.pgetsensornet(ndata, 1.25)
  edges = np.array(edges, dtype = np.uint64)
  edges = edges.reshape(-1,2)
  X = pytrimed.pmedoid(edges, datatype = "uu_edges", algorithm = "trimed", capture_output = False, maxcomputes = None)

def experiment_D_Sensor_net():
  """
  As per Table 1 of paper, column trimed (change algorithm argument for different algorithm
  """
  ndata = 3.6*10**5
  edges = pytrimed.pgetsensornet(ndata, 1.45)
  edges = np.array(edges, dtype = np.uint64)
  edges = edges.reshape(-1,2)  
  rev_edges = edges.copy()
  rev_edges[:,0] = edges[:,1]
  rev_edges[:,1] = edges[:,0]
  nedges = edges.shape[0]
  original = (npr.rand(nedges) > 0.5).reshape(-1,1)
  edges = np.array(edges*original + rev_edges*(1 - original), dtype = np.uint64)

  X = pytrimed.pmedoid(edges, datatype = "du_edges", algorithm = "trimed", capture_output = False, maxcomputes = None)


		
