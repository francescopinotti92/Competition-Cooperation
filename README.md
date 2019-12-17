# Competition-Cooperation
Here we provide the code required to reproduce results reported in the paper 

[Interplay between competitive and cooperative interactions in a three-player pathogen system](https://arxiv.org/abs/1912.07289)

F. Pinotti, F. Ghanbarnejad, P. Hövel, C. Poletto, arXiv:1912.07289 2019

**When you use this code, please cite the above reference**

Further details on terms of use: see LICENSE

# How to

1) Mean-Field simulations code/ contains python code required to integrate the mean-field equations for a single population and for two coupled populations for any given choice of initial conditions and model parameters. Also, it contains the code necessary to perform linear stability analysis for the single-population case

2) Stochastic simulations code/ contains the c++ code required to perform stochastic simulations on either Érdos-Renyj or modular networks. 

  2-a) Compile with ```make``` (c++11 compiler required)
  
  2-b) Run with ```./launch_simulation args```, where the ```args``` are:
    
    - Number of nodes
    
    - Maximum simulation time
    
    - Transmissibility B1 
    
    - Transmissibility B2
    
    - Transmissibility A
    
    - Cooperative factor C1
    
    - Cooperative factor C2
    
    - Recovery probability
    
    - Average network degree
    
    - Number of communities 
    
    - Probability of within-community stub
    
    - Index of community where B1 is initially injected 
    
    - Index of community where B2 is initially injected     
    
    - Index of community where A is initially injected     
