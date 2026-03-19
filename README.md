## Multi-objective optimum design of geometrically non-linear shallow dome structures using differential evolution-based algorithms

This repository presents a comparative study of multi-objective metaheuristic algorithms (MMIPDE, SHAMODE, and SHAMODE-WO) applied to structural optimization problems involving spatial truss structures.

The implemented methods are evaluated on benchmark truss problems with increasing complexity (24, 30, and 60 bars), allowing for a systematic analysis of algorithmic performance in multi-objective structural design.

---

## Algorithms

### MMIPDE
A population-based multiobjective metaheuristic that uses iterative parameter distribution estimation to guide the search process (Wansasueb et al., 2020). It incorporates Pareto dominance and maintains an external archive of non-dominated solutions.

### SHAMODE
A multiobjective evolutionary algorithm based on differential evolution strategies with adaptive parameter control (Panagant et al., 2019).

### SHAMODE-WO
A variant of SHAMODE with modified operators or reduced mechanisms, used to assess the contribution of specific algorithmic components (Panagant et al., 2019).

---

## Benchmark Problems

The algorithms are evaluated on three structural optimization problems:

- 24-bar truss
- 30-bar truss
- 60-bar truss

These problems represent increasing levels of structural complexity and dimensionality.

---

## Repository Structure

The repository is organized as follows:

MMIPDE:
- mmipde_24bars_truss/
- mmipde_30bars_truss/
- mmipde_60bars_truss/

SHAMODE / SHAMODE-WO:
- shamodes_24bars_truss/
- shamodes_30bars_truss/
- shamodes_60bars_truss/

Each folder includes:
- the implementation of the algorithm
- the structural problem definition
- constraint handling routines
- evaluation functions
- scripts for execution

---

## How to Run

1. Open MATLAB

You can run the code either using the MATLAB IDE or the command window.

- Using the MATLAB IDE:
  Open the desired folder and run the main script (e.g., `MMIPDE.m` or `main_shamodes.m`) by clicking "Run".

- Using the command window:

  cd 'path_to_repository/folder_name'

  % Run MMIPDE: MMIPDE

  % OR run SHAMODE: main_shamodes

Replace `folder_name` with one of the following:

- mmipde_24bars_truss
- mmipde_30bars_truss
- mmipde_60bars_truss
- shamodes_24bars_truss
- shamodes_30bars_truss
- shamodes_60bars_truss

---

## Experimental Setup

Typical parameters are defined within each algorithm's main script and may vary slightly between implementations:

- Population size (nsol)
- Number of iterations (nloop)
- Archive size (narchive)
- Number of independent runs (nrun)

Additionally, the number of objective functions (`NumFO`) and constraints (`QuantRestr`) must be properly defined in the `dados_do_problema.m` file according to the specific structural problem being considered.

The implementations also incorporate components of the NUMA-TF framework, as proposed by Rangel (2019).

```bibtex
@mastersthesis{rafaelrangel,
author = "Rangel, R. L.",
title = "Educational Tool for Structural Analysis of Plane
Frame Models with Geometric Nonlinearity",
year = "2019",
school = "Programa de Pós-graduação em Engenharia Civil, PUC-Rio", 
note ="In Portuguese"
}
```
---

## Output

The algorithms generate:

- Pareto-optimal solution sets
- Optimization history
- Result files in .mat format

---

## Platform Notes

Some implementations rely on compiled MEX files (e.g., structural solvers):

- Windows: precompiled .mexw64 files included
- Linux/macOS: recompilation required

Example:

mex Trelica_3d.c

---

## Reproducibility

To reproduce the results:

- Use consistent parameter settings
- Execute multiple independent runs
- Store and compare Pareto fronts
- (Optional) fix random seeds

---

## Citation

If you use this repository, please cite:

```bibtex
@article{mmipde2020,
  title={Multiobjective meta-heuristic with iterative parameter distribution estimation for aeroelastic design of an aircraft wing},
  author={Wansasueb, K. and Pholdee, N. and Panagant, N. and Bureerat, S.},
  journal={Engineering with Computers},
  pages={1--19},
  year={2020}
}
```

```bibtex
@article{panagant2019novel,
  title={A novel self-adaptive hybrid multi-objective meta-heuristic for reliability design of trusses with simultaneous topology, shape and sizing optimisation design variables},
  author={Panagant, Natee and Bureerat, Sujin and Tai, Kang},
  journal={Structural and Multidisciplinary Optimization},
  volume={60},
  number={5},
  pages={1937--1955},
  year={2019},
  publisher={Springer}
}
```
---

## Authors and Affiliations

- João Marcos de Paula Vieira, PPGMC/UFJF
- José Pedro Gonçalves Carvalho, PEC/COPPE/UFRJ
- Dênis Emanuel da Costa Vargas, DM/CEFET-MG
- Érica da Costa Reis Carvalho, DCOMP/UFSJ
- Patrícia Habib Hallak, MAC/UFJF 
- Afonso Celso de Castro Lemonge, MAC/UFJF
