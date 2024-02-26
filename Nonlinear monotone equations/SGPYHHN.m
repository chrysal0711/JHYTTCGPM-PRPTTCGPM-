% Matlab Model by Jianghua Yin (Mar.,2019,Nanning)
% Copyright (C) 2019 Jian Group
% All Rights Reserved

%% Derivative-Free Projection Method

% NO=1;
% method='PRP+';
 function [Tcpu,NF,Itr,NG] = SGPYHHN(NO,method) 
format long
Step=0:4;  % 生成一个5维行向量（0,1,2,3,4）
nextstep=Step(1);
% disp(date);
finish=0;
Itr=0;  %% 迭代
NF=1;   %% 目标函数计算次数
%NP=1;  %% 投影计算次数
tic     %% 初始化时钟
while finish==0
    loop=1;
    while loop==1
        switch nextstep
            case Step(1)
                % k=1;
                %% Step 0 初始化 参数定义
                %% 设置初始点
                [nprob,n,x0]=init(NO);
                epsilon=1e-6;
                epsilon1=1e-7;
                %% 线搜索参数
                sigma=0.01; % sigma=0.01
                gamma=0.8;    % gamma=1, initial guess for the steplength
                rho=0.4;    % rho=0.9
                %% 初始方向
                Fk=feval(nprob,n,x0,1);  % feval把n,x0赋值到一个定义好的函数nprob的情形1中
                dk=(-1)*Fk;
                nextstep=Step(2);
%               tic;
            case Step(2)
                %% Step 1 终止准则和线搜索
                if norm(Fk,2)<=epsilon || Itr>2000  %跳出循环
                    %                     disp('恭喜啦，迭代达到精度要求！')
                    %                     disp('xk是KKT点！！')
                    nextstep=Step(5);  % 跳出
                    break;
                else  % 线搜索
                 %% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    t =gamma;  % initial steplength
                    L1=1;
                    while L1==1
                        z_new=x0+t*dk;
                        Fz_new=feval(nprob,n,z_new,1);
                        NF=NF+1;
                        eta_k=max(0.001,min(0.8,norm(Fz_new,2)));  % 0.001,8
                      %% check the Armijo-type line search condition
                        if (-Fz_new'*dk < sigma*t*eta_k*norm(dk,2)^2 && t>10^(-10))  % the Armijo-type line search condition violated
                            L1=1;
                            t=t*rho;
                        else
                            L1=0;
                        end
                    end       % 终止while L1==1
                 %% End Armijo-type line search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     t
                    zk=x0+t*dk;  % zk=z_new;
                    Fzk= feval(nprob,n,zk,1);  %  Fzk=Fz_new;
%                     if norm(Fzk,2)<=epsilon
%                        nextstep=Step(5); % 跳出
%                        break;
%                     else
%                        nextstep=Step(3);
%                     end
                end  % 终止 if norm(Fk,2)<=epsilon || Itr>2000  
                nextstep=Step(3);
            case Step(3)
                xik=Fzk'*(x0-zk)/norm(Fzk)^2;
                zk1=x0-1.9*xik*Fzk;     % 此处通常的投影法系数是1
                x1=feval(nprob,n,zk1,2);     % 投影得到新的迭代点x_{k+1}
                %NP=NP+1;
                Fk0=Fk;
                Fk=feval(nprob,n,x1,1);   % 产生新的函数值F_{k+1}
                NF=NF+1;
                nextstep=Step(4);
            case Step(4)
                %% Step 2 更新
                Itr=Itr+1;
%               Fk0=feval(nprob,n,x0,1);
                sk=x1-x0;
                yk=Fk-Fk0;
                switch method
                    case 'NDYSPCGP'
                        mu=0.5;                                                                                  %0.1  和1没有区别  就继续压小  0.05稍微慢一点
                        theta=(abs(Fk0'*Fk))/(norm(Fk)*norm(Fk0));
                        if  -dk'*Fk0>dk'*yk
                            seg=-dk'*Fk0;
                        else
                            seg=dk'*yk;
                        end
                        if  Fk'*dk<mu*seg                         %限制条件变得更严格 
                            betak=(theta*(Fk'*Fk-(norm(Fk))/(norm(Fk0))*abs(Fk0'*Fk))+(1-theta)*(Fk'*Fk))/(seg);    
                            dk=(-1)*Fk+betak*dk;
                        else
                        betak=0;
                        dk=(-1)*Fk+betak*dk;  
                        end
                    case 'WYL'
                        betak=Fk'*(Fk-norm(Fk)/norm(Fk0)*Fk0)/norm(Fk0)^2; % WYL
                        dk=-1*Fk+betak*dk; 
                    case 'PRP+'
                        betak=max(0,Fk'*(Fk-Fk0)/norm(Fk0)^2);
                        dk=-1*Fk+betak*dk; 
                    case 'NN'
                        %Sun M, Liu J. New hybrid conjugate gradient projection method for the convex constrained equations[J].
                        %Calcolo, 2016, 53(3): 399-411.
                        mu=1.4;
                        betak=(norm(Fk)^2-max(0,norm(Fk)/norm(Fk0)*Fk'*Fk0))/max(mu*norm(dk)*norm(Fk),dk'*yk);
                        dk=-1*Fk+betak*dk; 
                    case 'MNN'
                        mu=1.4;
                        betak=(norm(Fk)^2-max(0,norm(Fk)/norm(Fk0)*Fk'*Fk0))/(max(mu*norm(dk)*norm(Fk),-dk'*Fk0));%max(norm(Fk0)^2,
                        dk=-1*Fk+betak*dk; 
                   case 'MN'
                        mu=1.4;
                        betak=(norm(Fk)^2-max(Fk'*Fk0,max(0,norm(Fk)/norm(Fk0)*Fk'*Fk0)))/(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2));%dk'*yk));
                        dk=-1*Fk+betak*dk; 
                   case 'YJH'
                        mu=1.4;
                        betak=(norm(Fk)^2-max(0,norm(Fk)/norm(Fk0)*Fk'*Fk0))/max(mu*norm(dk)*norm(Fk),dk'*yk);
                        dk=-(1+betak*Fk'*dk/norm(Fk)^2)*Fk+betak*dk; 
                   case 'JYJ' %投JCAM使用的共轭参数
                        mu=1.4;
                        betak=min(abs(Fk'*(Fk-Fk0))/norm(Fk0)^2,(norm(Fk)^2-norm(Fk)/norm(Fk0)*Fk'*Fk0)/(2*mu*norm(dk)*norm(Fk)));
                        dk=-1*Fk+betak*dk; 
                   case 'IFR'
                        betak=abs(Fk'*dk)/(-Fk0'*dk)*norm(Fk)^2/norm(Fk0)^2;
                        dk=-1*Fk+betak*dk; 
                   case 'YLS'
                        r=0.01;
                        yk1=yk+r*sk;
                        thetak=norm(sk)^2/(sk'*yk1);
                        dk=-1*thetak*Fk;
                   case'MNEWP' %参考自 Li M. A three term Polak-Ribière-Polyak conjugate gradient method close to the 
                       %memoryless BFGS quasi-Newton method, Journal of Industrial & Management Optimization, 2018.
                        tk=min(0.3,max(0,1-yk'*sk/norm(yk)^2));
                        mu=0.2;   % 数值效果还行 mu=0.2
                        thetak=tk*Fk'*dk/max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk),dk'*yk));  % 原始的分母是 norm(Fk0)^2  
                        betak=Fk'*yk/max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk),dk'*yk))-norm(yk)^2*Fk'*dk/(max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk),dk'*yk)))^2; %max(mu*norm(dk)*norm(yk),dk'*yk)
                        dk=-Fk+betak*dk+thetak*yk;
                   case'MNEW'
                        yk1=yk+0.01*dk;
                        tk=min(0.3,max(0,1-yk1'*sk/norm(yk1)^2));
                        mu=0.2;   % 数值效果还行 mu=0.2
                        fenmu=max(mu*norm(dk)*norm(yk1),dk'*yk1);
                        thetak=tk*Fk'*dk/fenmu;  % 原始的分母是 norm(Fk0)^2  
                        betak=Fk'*yk1/fenmu-norm(yk1)^2*Fk'*dk/fenmu^2; %(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2))
                        dk=-Fk+betak*dk+thetak*yk1;
                   case'HTTCG'
                        tk=min(0.3,max(0,1-yk'*sk/norm(yk)^2)); % yk1=yk+0.01*dk的数值效果没有直接使用yk效果好
                        mu=0.2;   
                        fenmu=max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk),dk'*yk));
                        thetak=tk*Fk'*dk/fenmu;  % 原始的分母是 norm(Fk0)^2  
                        betak=Fk'*yk/fenmu-norm(yk)^2*Fk'*dk/fenmu^2; %(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2))
                        dk=-Fk+betak*dk+thetak*yk;
                   case'HTTCGN'  %效果不错
%                         yk0=Fzk-Fk0;
                        lambdak=1+max(0,-yk'*sk/(norm(Fk)*norm(sk)^2));
                        yk1=yk+lambdak*norm(Fk)*sk;
                        tk=min(0.3,max(0,1-yk'*sk/norm(yk)^2)); % yk1=yk+0.01*dk的数值效果没有直接使用yk效果好
                        mu=0.2;   
                        fenmu=max(norm(Fk0)^2,max(mu*norm(dk)*norm(yk1),dk'*yk1));
                        thetak=tk*Fk'*dk/fenmu;  % 原始的分母是 norm(Fk0)^2  
                        betak=Fk'*yk1/fenmu-norm(yk1)^2*Fk'*dk/fenmu^2; %(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2))
                        dk=-Fk+betak*dk+thetak*yk1;
                   case'MMNEW'  % 效果没MNEWP和MNEW好
                        yk1=yk+0.01*sk;
                        tk=min(0.3,max(0,1-yk1'*sk/norm(yk1)^2));
                        mu=0.5;   % 数值效果还行 mu=0.2
                        fenmu=max(mu*norm(sk)*norm(yk),sk'*yk);
                        thetak=tk*Fk'*sk/fenmu;  % 原始的分母是 norm(Fk0)^2  
                        betak=Fk'*yk1/fenmu-norm(yk1)^2*Fk'*sk/fenmu^2; %(max(mu*norm(dk)*norm(Fk),norm(Fk0)^2))
                        dk=-Fk+betak*sk+thetak*yk1;
                   case 'PDY' %Liu J.K.,Feng Y.M.,(2018)A derivative-free iterative method for nonlinear monotone
                              %equations with convex constraints,Numerical Algorithms
                        tk=1+max(0,-dk'*yk/norm(dk)^2);
                        wk=yk+tk*dk;
                        betak=norm(Fk)^2/(dk'*wk);
                        c=1;
                        thetak=c-Fk'*dk/(dk'*wk);
                        dk=-thetak*Fk+betak*dk;
                   case 'TPRP' %2015 A three-terms Polak–Ribière–Polyak conjugate gradient algorithm for large-scale nonlinear equations
                        mu=0.01;
                        dk=-Fk+(Fk'*yk*dk-Fk'*dk*yk)/max(mu*norm(dk)*norm(yk),norm(Fk0)^2);
                end
                x0=x1; 
                if norm(dk)<epsilon1
                    nextstep=Step(5);
                    break;
                end
%                 t=t*(Fk11'*dk0)/(Fk'*dk);
%                 if (Fk'*dk>epsilon1)
%                     x1=x0+t*(-Fk);
%                 end
                nextstep=Step(2);
            case Step(5)
                loop=0;
%                 finish=1;
%                 break;
        end
    end
    if Itr>=2001 %|| norm(gk,2)>epsilon || isnan(Itr)==1  || abs(t)<epsilon  norm(gk,2)>epsilon 
      Tcpu=NaN;
      NF=NaN;
      Itr=NaN;
      NG=NaN;
    else
      Tcpu=toc; 
      NG=norm(Fk);
    end
    finish=1;
%      Tcpu=toc;
end
% [Tcpu Itr NF norm(Fk) norm(dk)]
%Tcpu=cputime-t0;
