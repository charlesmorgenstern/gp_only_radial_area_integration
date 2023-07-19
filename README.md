# gp_only_radial_area_integration
Test gp-only integration to find area of unit circle 

plotpaths(n,n2)
Plot domain, mock charge density, and gradient paths.
Use n gradient paths with n2 points each.

getpaths(n,n2,intradius)
Generate n gradient paths with n2 points starting on interior circle of radius intradius.

int_area(n,n2,intradius)
Approximate the area of a unit circle with gp-only integration. 
Use first order right hand rule to approximate the inetegral of the curvature.
Use second order trapezoidal rule to approximate the integral of l1(s) down each path.
Use n gradient paths with n2 points starting on interior circle of radius intradius.

errortable()
Create a table of relative error and order of convergence for the method.
Use 10 gradient paths. Very n the number of points on each path.
Fix radius of interior circle where paths begin at .05.


