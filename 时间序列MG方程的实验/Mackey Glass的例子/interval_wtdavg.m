%% interval_wtdavg.m 

%% Function used to implement the iterative procedure descibed in 
%% Theorem 9-1 to compute the maximum and minimum of a weighted average 
%% where both the "z_i" and "w_i" are interval sets.

%% Written by Nilesh N. Karnik - August 9,1998
%% For use with MATLAB 5.1 or higher.

%% Outputs : "l_out" and "r_out" (scalars) are, respectively, the
%% left and the right end-points of the wighted average "Y", which 
%% itself is an interval type-1 set.

%% Inputs : "c" and "s" are both M-vectors, containing, respectively, 
%% the centers and spreads of the "Z_i"s, and "h" and "delta" are  
%% M-vectors containing the centers and spreads of the "W_i"s. For each,
%% "W_i", "h >= delta", so that the support of each "W_i" is positive. 
%% This program is valid for minimum or product t-norms.

%% NOTE : Since type-reduction calculations can be represented as 
%% extended weighted average calculations, this program can also be used 
%% for type-reduction calculations for interval type-2 FLS's (see Chapter 9).

%% Uses "adapt.m".


function[l_out,r_out] = interval_wtdavg(c,s,h,delta)

lower = h - delta ;
upper = h + delta ;

l_out = adapt(c-s,lower,upper,-1) ; 
r_out = adapt(c+s,lower,upper,1) ; 

return ;

