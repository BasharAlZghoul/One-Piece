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




