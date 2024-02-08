# PROJ - Robust Shortest Path

This repository presents a collaborative work by Arthur Divanovic and Axel Navarro.

It contains all the necessary files for interpreting, using, and testing various resolution methods for the robust shortest path (RSP) problem (with uncertinaty in the objective and the constraints). The instances used for performance assessment are the ones provided in the 9th DIMACS challenge.

## Table of Contents

1. [Introduction](#1-introduction)
2. [Installation](#2-installation)
3. [Structure and Documentation](#3-structure-and-documentation)
4. [Use](#4-use)

## 1. Introduction

The goal of this repository is to gather all the functions required to solve the RSP problem thanks to two three big categories of resolution methods: some metaheuristics, dual reformulation methods and cutting planes methods.

This repository can be divided into two main folders:

- **src**: functions used for the solving of the RSP problem.
- **results**: results of the tests, grouped by instance.

## 2. Installation

This repository can be cloned directly from this webpage.

## 3. Structure and Documentation

### 3.1 src folder

This folder contains all the useful files to launch and evaluate the proposed resolution methods.

The file `main.jl` gathers all the necessary imports used throughout the project.

The rest of the folder is divided into two sub-folders:

#### 3.1.a Methods

- `antcolony.jl`: Adaptation of the Ant-Colony Heuristic for the RSP problem.
- `branchcut.jl`: Implementation of a resolution by branch & cut.
- `cutting.jl`: Implementation of a resolution by a cutting plane method.
- `djikstra.jl`: Implementation of an heuristic based on Djikstra's shortest path algorithm.
- `static.jl`: Implementation of the resolution of the static version of the RSP problem.

#### 3.1.b Utils

- `eval.jl`: Gathers the functions necessary for the rapid evaluation of a solution (objective and contraint violation).
- `graph.jl`: Defines the Graph structure and a parser to transform DIMACS data into Graph objects.
- `results.jl`: Gathers the functions to display and save the experimental results of all the methods considered.

### 3.2 Results

The Results folder contains `.txt` files. The names of the files correspond to the instances in the Data folder. In each result file, the best RSP solution found at the en dof the execution of a method is stored, along with the parameters employed for the method.

## 4. Use

Here is an example. 

1. Initialize a Graph object

```julia
g = parse_file("data/800_USA-road-d.COL.gr")
```

2.Initialize the parameters
```julia
save = true
time_limit = 8 * 60
```

3.Apply the method
```julia
obj_value, path, resolution_time = branch_and_cut_resolution(g, save, time_limit)
```

4. Display the results (path, objective value, respect of the robust constraint...)
```julia
display_results(g, path, resolution_time)
```
