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
end   

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
