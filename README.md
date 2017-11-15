# TRIMED

A C++ library and (optional) Python tool for obtaining the medoid of a set, based on
the paper https://arxiv.org/abs/1605.06950.


## PREREQUISITES

* CMake 
* for the Python library: Cython and Python


## CONFIGURE WITH CMAKE


Create a build directory:
```
mkdir build; cd build;
```

If you do NOT want the Python library, 

```
cmake -DBUILD_PYTHON_LIB=NO ..
```

If you do want the Python library, 


```
cmake ..
```

## BUILD

The library can be built, from the `build` directory 

```
make -j5

```

## EXAMPLES

In file `python/examples.py`
