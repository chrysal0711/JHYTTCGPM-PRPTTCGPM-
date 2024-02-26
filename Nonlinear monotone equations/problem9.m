% Matlab Model by Jiachen Jin (Dec.,2021,Nanning)
% Copyright (C) 2021 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.                     


function out= problem9(n,x,mode)
% Mor¨¦, J. J., Garbow, B. S., & Hillstrom, K. E. (1981). Testing unconstrained optimization software. 
%ACM Transactions on Mathematical Software, 7(1), 17-41.

if mode==1
  Fx=ones(n,1);
  Fx(1)=2*x(1)+(x(1)+(1/(n+1)))^3/(2*(n+1)^2)-x(2);
  for i=2:n-1
    Fx(i)=2*x(i)+(x(i)+(i/(n+1)))^3/(2*(n+1)^2)-x(i-1)+x(i+1);
  end
  Fx(n)=2*x(n)+(x(n)+(n/(n+1)))^3/(2*(n+1)^2)-x(n-1);
  out=Fx;
elseif  mode==2
%     p=zeros(n,1);
%     for i=1:n
%       if x(i)>=0;
%          p(i)=x(i);
%       else
%         p(i)=0;
%       end
%     end
    out=max(x,0);
%     if  min(x)>-1   %x(1)+x(2)+...+x(n)<=n && x(i)>=0,i=1,2...,n
%         out=x;
%     else
%         out=quadprog(speye(n),-x,ones(1,n),n,[],[],-ones(n,1));
%         %out=quadprog(eye(n),-x,[ones(1,n);zeros(n-1,n)],[n;zeros(n-1,1)],[],[],-1*ones(n,1),[]);
%     end
end

    
        
    
    


