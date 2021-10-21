function[outextreme,count,theta] = adapt(ypoint,lower,upper,maxflag)
[ysort,lower_sort,upper_sort] = trimvec(ypoint,lower,upper,1) ;
ly=length(ysort) ;

hl = (lower_sort+upper_sort)/2 ;
S = sum(ysort.*hl)/sum(hl) ;   % starting point

eps = 1e-5 ;  

count = 0 ;
theta = hl ;
S_new = S + 10*eps ;

if ((abs(S-ysort(1)) < eps) || (abs(S-ysort(ly)) < eps))
   outextreme = S ;

else

   while (abs(S-S_new) > eps) 
      count = count + 1;

      if count > 1
         S = S_new ;
      end   %%% if count

      in1 = find(ysort > (S-eps)) ;
      min1 = min(in1) ;

      if min1 > 2
         in2 = 1 : min1-1;
      else
          in2 = 1;
      end   
   
      if maxflag > 0
         theta(in1) = upper_sort(in1) ;
         theta(in2) = lower_sort(in2) ;
      else
         theta(in1) = lower_sort(in1) ;
         theta(in2) = upper_sort(in2) ;

                     % To avoid division by zero if all lower_sort=0
         if abs(S - ysort(min1)) < eps
                theta(min1) = upper_sort(min1) ; 
         end  %%% if abs(S_new ...

      end    %%% if maxflag
   
   
      S_new = sum(ysort.*theta)/sum(theta) ;

   end    %%% while

   outextreme = S_new ;

end % if ((abs(S-ysort(1)) < eps) .......

return ;
