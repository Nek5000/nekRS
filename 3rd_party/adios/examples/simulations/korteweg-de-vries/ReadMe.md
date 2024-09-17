### ADIOS2 korteweg-de-vries example

Solves the initial value problem for the Korteweg de-Vries equation via the Zabusky and Krustal scheme.

#### Running the example

```
./adios2_simulations_kortewegDeVries N t_max δ
   N is number of spacial gridpoints (∆t chosen from ∆x via the stability condition)
   t_max is max simulation time
   δ is  an interaction parameter;
   
e.g., ./adios2_simulations_kortewegDeVries 512 10 0.022
```

#### View data generated:

```
$ conda create -n adios2_simulations_kortewegDeVries
$ conda activate adios2_simulations_kortewegDeVries
$ conda install -c williamfgc adios2 -c conda-forge
$ conda install pip
$ pip install diagram
```
