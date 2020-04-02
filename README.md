# tauRAMD
Computation of the drug-target relative residence times from RAMD simulations:
this script helps to analyse results of Random Acceleration Molecular Dynamics (RAMD) simulatiuons to estimate relative residence time of a protein-ligand complex. 

    Author: Daria Kokh
    v.1.0
    Copyright (c) 2019
    Released under the GNU Public Licence, v2 or any higher version

    Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
    Schloss-Wolfsbrunnenweg 35
    69118 Heidelberg, Germany
    
    PREREQUISITE
    
    RAMD simulations were performed using either Gromacs or NAMD software (in the later case )
    
    USAGE
          python tauRAMD-v1.py  input_file[s]
          input files must contain a set of lines extracted from the gromacs output. Each line contains the number of steps executed before dissociation 
          and has the following format: “XX/YYYY.out:==== RAMD ==== GROMACS will be stopped after 874650 steps.”
          or if RAMD simulations were done using NAMD software: "EXIT: 462950  > LIGAND EXIT EVENT DETECTED: STOP SIMULATION"
    OUTPUT
         residence time with the standard deviation computed for each input_file and an image with histogram representation of the bootstrapping output
    ''')
    
This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020  Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).
