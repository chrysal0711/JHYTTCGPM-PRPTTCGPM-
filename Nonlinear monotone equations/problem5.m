% Matlab Model by Jianghua Yin (Nov.,2015,Yulin)
% Copyright (C) 2015 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem5(n,x,mode)
% adapted from W. La Cruz. A spectral algorithm for large-scale systems of nonlinear 
%monotone equations. Numerical Algorithms, 76(4): 1109-1130, 2017.

if mode==1
%   Fx=ones(n,1);
%   for i=1:n
%     Fx(i)=x(i)-sin(abs(x(i)-1));
%   end
    out=2*x-sin(abs(x));
elseif mode==2
%    n=length(x);
    out=max(x,0);
end


    
        
    
    


