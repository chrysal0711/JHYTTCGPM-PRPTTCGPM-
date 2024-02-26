% Matlab Model by Jianghua Yin (Nov.,2015,Yulin)
% Copyright (C) 2015 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem6(n,x,mode)
% adapted from La Cruz, W., Raydan, M.: Nonmonotone spectral methods for large-scale nonlinear systems. Optim.
%Methods Softw. 18, 583¨C599 (2003)

if mode==1
%   Fx=ones(n,1);
%   for i=1:n
%     Fx(i)=x(i)-sin(abs(x(i)-1));
%   end
    out=log(abs(x)+1)-x/n;
elseif mode==2
%    n=length(x);
    out=max(x,0);
end


    
        
    
    


