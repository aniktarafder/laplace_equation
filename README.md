# laplace_serial
The code was originally written by **John Urbanic, PSC**. I modified how the arrays are defined and added one subroutine to plot **Temperature** field in **Tecplot** format.

This code solves laplace heat equation in a 2D domain using jacobian iterations. 
* Left and bottom boundary --> 0 Temperature.
* Top boundary   --> Temperature varies from 0 to 100 linearly from left to right.
* Right boundary --> Temperature varies from 0 to 100 linearly from bottom to top.

Just clone the **repo** --> **compile** --> **run**

The code will generate 1 **.dat** file of temperature fields.

* **To clone**: git clone https://github.com/aniktarafder/laplace_serial.git
* **To compile** : gfortran -o laplace_serial.exe laplace_serial.f90
* **To run**     : ./laplace_serial.f90

**Course Link:** https://www.psc.edu/resources/training/xsede-hpc-workshop-may-4-5-2021-mpi/
