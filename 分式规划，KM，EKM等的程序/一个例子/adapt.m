%% adapt.m  

%% Function used to implement the iterative procedure described in
%% Theorem 9-1 to compute the maximum and minimum of a
%% weighted average, where the "y_i" s are crisp numbers and
%% the "w_i"s are interval sets that take values from some
%% real interval "[w_lower,w_upper]". 
 
%% Written by Nilesh N. Karnik - August 9,1998
%% For use with MATLAB 5.1 or higher.

%% Outputs : "outextreme" is the extreme value of the output. "count" 
%% is the number of iterations required to reach the optimum, and 
%% "theta" is the combination of "w_l"s that achieves the optimum.
 
%% Inputs : "ypoint", "lower" and "upper" are all M-dimensional vectors. 
%% "ypoint" contains the "y_l"s. "lower" and "upper" 
%% contain, respectively, the "w_lower" and "w_upper" values for each
%% weight "w_l". If "maxflag > 0" (scalar), "S" is maximized, else 
%% it is minimized.
 
%% Uses "../OPERATIONS/trimvec.m".

%% NOTE : The "addpath" and "rmpath" commands in the beginning and end 
%% of the function must be modified, depending upon the directory in 
%% which the function "trimvec.m" is stored.


function[outextreme,count,theta] = adapt(ypoint,lower,upper,maxflag)

%addpath ../OPERATIONS ;   %% for the function "trimvec.m"

[ysort,lower_sort,upper_sort] = trimvec(ypoint,lower,upper,1) ;
ly=length(ysort) ;

hl = (lower_sort+upper_sort)/2 ;
S = sum(ysort.*hl)/sum(hl) ;   % starting point

eps = 1e-5 ;   %%% small quantity to avoid floating point equality problems

count = 0 ;
theta = hl ;
S_new = S + 10*eps ;

if ((abs(S-ysort(1)) < eps) | (abs(S-ysort(ly)) < eps)),
   outextreme = S ;

else

   while (abs(S-S_new) > eps) ,
      count = count + 1;

      if count > 1,
         S = S_new ;
      end   %%% if count

      in1 = find(ysort > (S-eps)) ;
      min1 = min(in1) ;

      if min1 > 2,
         in2 = [1 : min1-1] ;
      else in2 = 1;
      end   %%% if
   
      if maxflag > 0,
         theta(in1) = upper_sort(in1) ;
         theta(in2) = lower_sort(in2) ;
      else
         theta(in1) = lower_sort(in1) ;
         theta(in2) = upper_sort(in2) ;

                     % To avoid division by zero if all lower_sort=0
         if abs(S - ysort(min1)) < eps,
                theta(min1) = upper_sort(min1) ; 
         end  %%% if abs(S_new ...

      end    %%% if maxflag
   
   
      S_new = sum(ysort.*theta)/sum(theta) ;

   end    %%% while

   outextreme = S_new ;

end % if ((abs(S-ysort(1)) < eps) .......

%rmpath ../OPERATIONS ;  %% Remove added path

return ;
