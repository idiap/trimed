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

import numpy as np
import numpy.random as npr
import random
import multiprocessing as mpr
from libcpp.string cimport string
from libcpp.vector cimport vector 
from libcpp.memory cimport unique_ptr # can you get this to work with unique_ptr to array ?
from libcpp.pair cimport pair
from libcpp cimport bool
cimport cython
cimport cython.floating

cdef extern from "trimed/medoidalgs.hpp" namespace "medoidalgs":		
  void medoid_array(size_t ndata, size_t dimension, const double * const data, double & min_E, size_t & min_E_i, string algorithm, bool caputre_output, string & text, size_t seed, size_t maxcomputes) except+;
  void medoid_uugraph(size_t nedges_in, const size_t * const edges_in, double & min_E, size_t & min_E_i, string algorithm, bool caputre_output, string & text, size_t seed, size_t maxcomputes) except+;
  void medoid_dugraph(size_t nedges_in, const size_t * const edges_in, double & min_E, size_t & min_E_i, string algorithm, bool capture_output, string & text, size_t seed, size_t maxcomputes) except+;
  

def dangerwrap(f):
  """
  I assume f is a function which returns an object and takes no parameters
  """
  event  = mpr.Event()
  q = mpr.Queue()
  
  def signalling_f():
    try:
      q.put(f())
    
    except Exception as e:
      print "Caught exception in dangerwrap:"
      print e
      q.put(e)
      event.set()
      return None
      
    event.set()
  
  f_process = mpr.Process(target = signalling_f)
  f_process.start()
  try:
    event.wait()
  
  except KeyboardInterrupt:
    f_process.terminate()
    f_process.join()
    raise KeyboardInterrupt("Caught KeyboardInterrupt in dangerwrap")
  
  return q.get()


def base_medoid(nrows, ncols, datatype, cython.double [:] pointdata, size_t [:] edgedata, algorithm, capture_output, seed, maxcomputes):
   
  cdef double min_E = 0
  cdef size_t min_E_i = 0
  cdef string text = ""
  if datatype == "uu_edges":
    medoid_uugraph(nrows, &edgedata[0], min_E, min_E_i, algorithm, capture_output, text, seed, maxcomputes)	
  
  if datatype == "du_edges":
    medoid_dugraph(nrows, &edgedata[0], min_E, min_E_i, algorithm, capture_output, text, seed, maxcomputes)
  
  elif datatype == "points":
    medoid_array(nrows, ncols, &pointdata[0], min_E, min_E_i, algorithm, capture_output, text, seed, maxcomputes)
  
  return {'text':text, 'min_E': min_E, 'min_E_i': min_E_i}
  
	
def pmedoid(data, datatype = None, algorithm = "trimed", capture_output = False, seed = None, maxcomputes = None):
  """
  data
  data to find medoid of, format depends on datatype (see next entry)
  
  datatype
  "points" / "uu_edges" / "du_edges"
  
  algorithm: 
  "trimed" / "RAND" / "TOPRANK1" / "TOPRANK2"
  
  capture_output:
  capture the output, returned as string in dict 
  
  seed:
  set the random seed
  
  maxcomputes:
  used in RAND algorithm, max number of points to compute. 
  
  """

  algorithms = ["trimed", "RAND", "TOPRANK1", "TOPRANK2"]
  if algorithm not in algorithms:
    errm =  "algorithm should be one of [ "
    for x in algorithms :
      errm += x 
      errm += " "
    errm += "] and not "
    errm += algorithm 
    raise RuntimeError(errm)

  if algorithm == "RAND" and not maxcomputes:
    raise RuntimeError("with RAND, maxcomputes should be provided")
    
  if seed == None:
    seed = npr.randint(10**9)
    
  if datatype == None:
    raise RuntimeError("datatype is required explicitly")
  
  pointdata = None
  edgedata = None
  nrows, ncols = data.shape	
  if datatype == "points":
    if data.dtype != np.float64:
      raise RuntimeError("Currently points datatype only support double (float64)")
    
    pointdata = data
    edgedata = np.array([7], dtype = np.uint64)
    
  elif datatype == "uu_edges" or datatype == "du_edges" :
    if data.dtype != np.uint64:
      raise RuntimeError("Currently du_edges / uu_edges datatype only support size_t (uint64)")
    if ncols != 2:
      raise RuntimeError("the array passed in needs to have 2 columns if datatype is `du_edges' or `uu_edges' ")
    
    pointdata = np.array([7.], dtype = np.float64)
    edgedata = data
  
  else:
    raise RuntimeError("Unrecognised datatype, should be edges or points")
    
  if maxcomputes == None:
    maxcomputes = np.uint64(2**64 - 1)
  
  
  return dangerwrap(lambda : base_medoid(nrows, ncols, datatype, pointdata.ravel(), edgedata.ravel(), algorithm, capture_output, seed, maxcomputes))
	


###########################

cdef extern from "trimed/graphtools.hpp" namespace "graphtools":		
  vector[size_t] getsensornet(size_t ndata, double phi_r) except +;
  pair[vector[size_t], vector[double]] convert_map(size_t nnodes, size_t nedges, const double * const nodes, const double * const edges) except +;
  
cdef extern from "<utility>" namespace "std":
  vector[size_t] move(vector[size_t]) # Cython has no function templates
    
def base_convert_map(nnodes, nedges, cython.double [:] nodes, cython.double [:] edges):
  cdef pair[vector[size_t], vector[double]] retval
  retval = convert_map(nnodes, nedges, &nodes[0], &edges[0])
  return {'ID_edges': np.array(retval.first).reshape(-1, 2), 'edge_lengths': np.array(retval.second)}
  
def pconvert_map(nodes, edges):
  nnodes = nodes.shape[0]
  nedges = edges.shape[0]
  return dangerwrap(lambda : base_convert_map(nnodes, nedges, nodes.ravel(), edges.ravel()))
  
def base_getsensornet(ndata, phi_r):
  cdef vector[size_t] up = move(getsensornet(ndata, phi_r)) #unique_ptr up = 	
  print "vector made"
  return np.array(up)
  
  
def pgetsensornet(ndata, phi_r):
  """
  generate sensor network data (see examples.py for use).
  """
  
  return dangerwrap(lambda : base_getsensornet(ndata, phi_r))
  
  
