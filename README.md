# Interactive Evolutionary Multiobjective Optimization of Primer Design with Uncertain Objectives

This repository contains the implementation for the paper [Interactive Evolutionary Multiobjective Optimization of Primer Design with Uncertain Objectives](https://dl.acm.org/doi/abs/10.1145/3638529.3654167) (**GECCO 2024**) for finding optimal primer designs interactively, with uncertainty in the objective, and for estimating melting temperatures. 

### Data and Results
* Results of the experiments can be found in the `primer_opt_results` folder.
* The dataset used for the tests is provided in the `primer_data` folder.

### Clone Repository
`gh repo clone amrzr/MOPrimer_Interactive`

### Requirements:
* The project relies on the desdeo-emo and desdeo-problems packages for the optimization algorithm used and defining the optimization problem.
* Python 3.7 or up

### Installation Process:
* Create a conda environment by
`conda env create -f environment.yml`

* Activate the conda environment
` conda activate moprimer_interactive`

### Optimizing Primer Design
* Run the `primer_optimize_intr_refpnt.py` file to find optimal primer designs interactively with reference point type preference. The plots of the solutions are saved as .html files and can be viewed with a web browser.
*  You can read other DNA sequences by uncommenting the fastaparser lines and providing the ID. Or you can simply define the DNA sequence.
