1. # One Piece Module for Particle Tracking

## Model Description
The presented module represents a sophisticated particle tracking model meticulously crafted for the micro-scale analysis of colloidal particles within complex porous media. This model excels in capturing colloids attached to medium grain surfaces under varying conditions, including both favorable and unfavorable scenarios.

Moreover, the module possesses the capability to compute several critical quantities, including but not limited to:

1. Attachment Location (Retention Profile)
2. Bulk Residence Time
3. Near Surface Residence Time (the duration a colloid spends within the near-surface zone, precisely 200 nm from the grain surface, to which it attaches)
4. The Number of Grains Interacted with Prior to Attachment
5. Breakthrough Curve
6. Attachment Angle (relative to the front flow stagnation zone of the corresponding grain)

This module, designed for particle tracking, serves as an invaluable tool for in-depth analyses of colloidal behavior within complex porous environments, offering a comprehensive suite of features for advanced research and understanding.



## Model Inputs

### dkdkd

## How to Compile and Run this Model

This model is written in Parallel using MPI in Fortran to reduce computational time as much as possible. It can be compiled using any suitable compiler. 



