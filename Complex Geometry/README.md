# Model Inputs

Given the Fortran-based architecture of this model, all inputs are expected to be provided in the form of text files. Within the input directory, a Matlab code is available to facilitate the conversion of model inputs from Matlab matrices to text files, ensuring compatibility with the model.


## Ufx.txt and Ufy.txt

These files represent the velocity field in both the x and y directions. It is essential to adhere to the correct format, where Ufx(x, y) is the appropriate structure (NOT Ufx(y, x)), a stipulation that also applies to Ufy.


## D.txt 

This input file serves to describe the geometry, encompassing information about the locations and sizes of grains.
Note: This model can handle up to 3 different grain sizes within the geometry.

## xp.txt and yp.txt

These files detail the locations of colloids within the system:

xp.txt contains the x-coordinate colloid locations.
yp.txt contains the y-coordinate colloid locations.


# Model Outputs

Three .dat files are the outputs of the fortran code. 
Outputs directory contains a Matlab code to extract the outputs of the simulation as matrices.

The output matrices will be:
1. attached.mat: has information about attached colloids
2. pass.mat: has information about the colloids that exit the domain without being attached.

both of these matrices have 23 columns for different quantities. Here is a description for some of these quantities for attached.mat matrix:


column 1: attachment x position

column 2: attachment y position

column 3: total residence time

column 4: number of garins this colloid interacted with

column 5: near surface residence time (NSRT) on the grain of attachment (s)

column 6: NSRT on 1st grain this colloid interacted with

column 7: NSRT on 2nd grain this colloid interacted with

column 8: NSRT on 3rd grain this colloid interacted with

column 9: NSRT on 4th grain this colloid interacted with

column 10: accumulativge NSRT on all garins this colloid interacted with

column 19: The round number that this colloid attached at (since this model can be applied to periodic geometry, this value is important to compute the colloid traveled distance)

