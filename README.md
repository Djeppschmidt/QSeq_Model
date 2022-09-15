# QSeq_Model
Analysis and models for QSeq paper

The file folder "Scripts" has two script files. QSeq_betadiv_analysis.R has the scripts used to compare QSeq methods on previously published datasets. benchmarking.R contains the scripts for replicating the simulations in the paper.

The file folder "Data" contains phyloseq community datasets derived from each of the three studies, each ending in "comdat.RDS". Files containing "initsksp" are the initial Reference matrices for the benchmarking workflow. files with "model" are the results of specific parameterized runs of the simulation. These are each included so that the figures from the paper can be reproduced exactly; however the script file also shows how to freshly generate random data, and rerun the simulation from the Reference data generation step.
