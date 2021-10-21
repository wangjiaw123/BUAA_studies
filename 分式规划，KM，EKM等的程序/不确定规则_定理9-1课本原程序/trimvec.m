%% trimvec.m 
 
%% This function trims a vector, i.e., sorts it and removes repeated 
%% elements. The vector "x1" is sorted and repeated elements in it
%% are removed. Vectors "y1" and "y2" are just rearranged like "x1". 
%% The elements of "y1" and "y2" corresponding to the repeated 
%% elements of "x1" are either summed or only the maximum of them 
%% is kept, depending upon the value of the flag "max_or_sum".
 
%% Written by Nilesh N. Karnik - August 9,1998
%% For use with MATLAB 5.1 or higher.

%% Outputs : The trimmed vectors "xtr", "ytr1" and "ytr2" (see the 
%% last paragraph).

%% Inputs : "x1", "y1" and "y2" are vectors of the same size. "x1" 
%% is sorted and trimmed. "y1" and "y2" are just rearranged in the 
%% same way as "x1". For both "y1" and "y2" : if "max_or_sum < 0" 
%% (scalar), only the maximum of the elements corresponding to the
%% repeated elements in "x1" is kept, else the elements corresponding 
%% to the repeated elements in "x1" are summed. The default value
%% for "max_or_sum" is -1.

%% NOTE 1 : This function can also be used to just trim one vector, as 
%% "xtr = trimvec(x1)", or when there is only one y-vector, as 
%% "[xtr,ytr1] = trimvec(x1,y1)". In the latter case, the default value 
%% for "max_or_sum" is -1 (i.e., by default, the maximum of the 
%% elements in "y1" corresponding to the repeated elements in "x1" is 
%% kept).
 
%% NOTE 2 : This function is used by the following M-files :
%% extend_binary.m
%% extend_nary1.m
%% extend_nary2.m
%% ../INTERVAL/adapt.m

function[xtr,ytr1,ytr2] = trimvec(x1,y1,y2,max_or_sum) 

if nargin == 1,
   y1 = ones(size(x1)) ;
   y2 = ones(size(x1)) ;
   max_or_sum = -1 ;
elseif nargin == 2,
   y2 = ones(size(x1)) ;
   max_or_sum = -1 ;
elseif nargin == 3,
   max_or_sum = -1 ;
end    % if nargin

[x,in1] = sort(x1) ;
y = y1(in1) ;
y22 = y2(in1) ;
lx = length(x) ;

xd = diff(x) ;

in2 = find(xd == 0) ;
lin2 = length(in2) ;



if ~isempty(in2), 

    if in2(1) > 1,
       xtr = x(1:in2(1)-1) ;
       ytr1 = y(1:in2(1)-1) ;
       ytr2 = y22(1:in2(1)-1) ;
    else 
       xtr = [] ;
       ytr1 = [] ;
       ytr2 = [] ;
    end   %%% if in2(1)
    
    
    for i1 = 1 : lin2,
        in3 = find(x == x(in2(i1))) ;
    
        if i1 > 1,
           if x(in2(i1)) ~= x(in2(i1-1)),
                   xtr = [xtr x(in2(i1))] ;
                 
                   if max_or_sum < 0,
                          ytr1 = [ytr1 max(y(in3))] ;
                          ytr2 = [ytr2 max(y22(in3))] ;
                   else
                          ytr1 = [ytr1 sum(y(in3))] ;
                          ytr2 = [ytr2 sum(y22(in3))] ;
                   end  %%% if max_or_sum
           end % if x(..)
        else
           xtr = [xtr x(in2(i1))] ;

           if max_or_sum < 0,
                  ytr1 = [ytr1 max(y(in3))] ;
                  ytr2 = [ytr2 max(y22(in3))] ;
           else
                  ytr1 = [ytr1 sum(y(in3))] ;
                  ytr2 = [ytr2 sum(y22(in3))] ;
           end  %%% if max_or_sum
        end   %%% if i1
        

        if i1 < lin2,
           if in2(i1+1) - in2(i1) > 2,
              xtr = [xtr x(in2(i1)+2 : in2(i1+1)-1)] ;
              ytr1 = [ytr1 y(in2(i1)+2 : in2(i1+1)-1)] ;
              ytr2 = [ytr2 y22(in2(i1)+2 : in2(i1+1)-1)] ;
           end  %%
        end   %% if i1
    end  %% for i1
    
    if in2(lin2) < lx-1,
           xtr = [xtr x(in2(lin2)+2 : lx)] ;
           ytr1 = [ytr1 y(in2(lin2)+2 : lx)] ;
           ytr2 = [ytr2 y22(in2(lin2)+2 : lx)] ;
    end %% if

else
    
    xtr = x ;
    ytr1 = y ;
    ytr2 = y22 ;

end  %%% if ~isempty(in2)

return ;
