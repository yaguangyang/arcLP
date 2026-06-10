

   AUTHORS:

       Yaguang Yang
       yaguang.yang@verizon.net


   REFERENCE:

    -  Y. Yang, “Two computationally efficient polynomial-iteration 
       infeasible interior-point algorithms for linear programming”, 
       Numerical Algorithms, Vol. 79 (2018), pp. 957–992.
    -  Y. Yang, “ArcLp: A MATLAB implementation of an O(\sqrt{n}L)
       arc-search infeasible interior-point algorithm for linear 
       programming", Journal of Open Research Software, Vol. xx (2026),
       pp. xxx-xxx.


   SOFTWARE REVISION DATE:

       V1.0, February 2026

   SOFTWARE LANGUAGE: 

       MATLAB



=====================================================================
PACKAGE
=====================================================================

The directory contains the following files

README              :  this file
arcLP.m           :  MATLAB program using arc-search to solve linear 
                    :  programming problems  
mehrotra.m          :  MATLAB program using Mehrotra's algorithm to solve 
                    :  linear programming problems
extract.m           :  MATLAB program used to extract (A,b,c) from Netlib
                    :  test problems
infeasible_method.m :  MATLAB program calls either arcLp.m or 
                    :  mehrotra.m and returns iter, obj, and infe
lp_*.mat            :  Netlib test problems which contain detailed 
                    :  information about the problems, including (A,b,c)

infeasibleexample.m is a MATLAB code which generates Figure 2 of the paper.

=====================================================================
HOW TO INSTALL
=====================================================================

Download and unpack the zip archive. A folder containing the package 
files will be created.

================
RUN THE PROGRAMS
================

1. Open MATLAB and change to the package folder
2. Load a standard LP problem listed in the paper, for example,
	>> load lp_adlittle
3. Extract matrix A, and vectors b and c by calling the MATLAB function 
   extract.m
	>> extract
4. Run arcLP.m as follows, assuming degeneracy is not a problem
	>> [x,obj,kk,infe]=arclp(A,b,c);
5. Run mehrotra.m as follows, assuming degeneracy is not a problem
	>> [x,obj,kk,infe]=mehrotra(A,b,c);

If degeneracy is a problem, Steps 4 and 5 are replaced by the following 
4'. Run arcLP.m as follows 
	>> [x,obj,kk,infe,m1,n1]=arclp(A,b,c,1);
5'. Run mehrotra.m as follows 
	>> [x,obj,kk,infe]=mehrotra(A,b,c,1);

====================
A degeneracy example
====================
>> load lp_degen3.mat
>> extract
>> [x,obj,kk,infe]=arcLP(A,b,c);
>> [x,obj,kk,infe]=arcLP(A,b,c,1);
The arcLP(A,b,c) without option d=1 cannot find the solution 
after 32 iterations because the problem is degenerate. But if
arcLP(A,b,c,1) is called, it finds the solution in 24 
iterations with the infeasibility measure infe = 6.9339e-09.
  
