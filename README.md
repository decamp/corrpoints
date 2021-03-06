### CorrPoints:
Fits transform functions to sets of 2D CORRespondence POINTS. This library currently supports the following function families:

- linear
- polynomial (implementations for 2, 3 and N degrees)
- projection 

**Why?**

Say you need to warp an image to correct for camera distortion. You may not have the functions that define the lens optics,
but in a photograph taken by the camera, you might identify several points and where they should be located in the corrected
photo. Given a small set of such points, this library will generate a smooth function that can be used to map all the other
points.


### Build:
$ ant


### Runtime:
After build, add all jars in **lib** and **target** directories to your project.


### Dependencies:

**jlapack**

JLapack is a pure java version of the BLAS and LAPACK math libraries. 

Distributed under a BSD license, which is located in the "LICENSE.txt" file.
May be located online by looking for "f2j" and "lapack". 

This version is machine generated and, unfortunately, the process converts Fortran directly into Java bytecode without creating any source code. The version of jlapack I have included uses the non-strict math functions, so if numerical accuracy
is required, you may opt to locate a version that uses strict math. It also normally comes as four separate jars, which are located in 
the "resources_build" folder, but I've combined them into a single jar for my own convenience.

---
Author: Philip DeCamp
