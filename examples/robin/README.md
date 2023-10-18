# Unsteady conduction with Robin boundary conditions.

This is a 1-element unsteady conduction example that
illustrates the use of the convection heat-transfer
boundary condition, which is based on the Newton
law of cooling.

For this boundary condition, the thermal flux at
the wall is proportional to the difference between
the wall temperature, Tw, and the ambient temperature,
Tinf (T_infinity), with proportionality coefficient
given by the heat transfer coefficient, hc.  That is,
at the wall,

      -div k grad T = hc(T-Tinf)

The thermal boundary condition is specified in Nek5000
as "c  " when hc and Tinf are given in userbc(), or 
"C  " when hc and Tinf are given as constant parameters
in the .rea file.  For example, if Element 1, Side 2, 
were to have Tinf=5 and hc=9, the .rea file would have

  ***** THERMAL BOUNDARY CONDITIONS *****
 :    :  :  
 :    :  :  
 C    1  2  5.000000      9.000000     0.0000000E+00 0.0000000E+00 0.0000000E+00
 :    :  :  
 :    :  :  

This condition could also be set via userbc() with the following
lines in the .rea and .usr files.

robin.rea:

  ***** THERMAL BOUNDARY CONDITIONS *****
 :    :  :  
 :    :  :  
 c    1  2  5.000000      9.000000     0.0000000E+00 0.0000000E+00 0.0000000E+00
 :    :  :  

robin.usr:

      subroutine userbc (ix,iy,iz,iside,ieg)
      :
      :
      tinf=5.0
      hc=9.0
      :


The latter approach is used in the robin example given
here.

