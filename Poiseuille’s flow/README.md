# Model Inputs

Given the Fortran-based architecture of this model, all inputs are expected to be provided in the form of text files. Within the input directory, a Matlab code is available to facilitate the conversion of model inputs from Matlab matrices to text files, ensuring compatibility with the model.


## yp.txt

This file detail the y location of colloids within the system:
The x location of colloids can be specified in the main code


# Model Outputs

Three .dat files are the outputs of the fortran code. 
Outputs directory contains a Matlab code to extract the outputs of the simulation as matrices.

The output matrices will be:
1. attached.mat: has information about attached colloids
2. pass.mat: has information about the colloids that exit the domain without being attached.

both of these matrices have 8 columns for different quantities. Here is a description for the important quantities for attached.mat matrix:


column 1: attachment x position

column 2: attachment y position

column 3: total residence time


