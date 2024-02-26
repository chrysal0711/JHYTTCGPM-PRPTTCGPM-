% Matlab Model by Jianghua Yin (Nov.,2015,Yulin)
% Copyright (C) 2015 Jian Group
% All Rights Reserved
%%% 谱梯度投影法
% NO=1;
% method='PRP+';
 function [Tcpu,NF,Itr,NG] = SGP(NO,method) 
format long
Step=0:4;
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
                gamma=1;  % gamma=1, initial guess for the steplength
                rho=0.5;    % rho=0.9
                %% 初始方向
                Fk=feval(nprob,n,x0,1); % feval把n,x0赋值到一个定义好的函数nprob的情形1中
                dk=(-1)*Fk;
                nextstep=Step(2);
%               tic;
            case Step(2)
                %% Step 1 终止准则和线搜索
                if norm(Fk,2)<=epsilon || Itr>2000  %跳出循环
                    %                     disp('恭喜啦，迭代达到精度要求！')
                    %                     disp('xk是KKT点！！')
                    nextstep=Step(5); % 跳出
                    break;
                else % 线搜索
                 %% %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    t =gamma;  % initial steplength
                    L1=1;
                    while L1==1
                        z_new=x0+t*dk;
                        Fz_new=feval(nprob,n,z_new,1);
                        NF=NF+1;
                        % check the Armijo-type line search condition
                        if (-Fz_new'*dk < sigma*t*norm(dk,2)^2 && t>10^(-10)) % the Armijo-type line search condition violated
                            L1=1;
                            t=t*rho;
                        else
                            L1=0;
                        end
                    end       % 终止while L1==1
                 %% %%% End Armijo-type line search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     t
                    zk=x0+t*dk;
                    Fzk= feval(nprob,n,zk,1);
                    nextstep=Step(3);
%                     if norm(Fzk,2)<=epsilon
%                        nextstep=Step(5); % 跳出
%                        break;
%                     else
%                        nextstep=Step(3);
%                     end
                end  % 终止 if norm(Fk,2)<=epsilon || Itr>2000      
            case Step(3)
                xik=Fzk'*(x0-zk)/norm(Fzk)^2;
                zk1=x0-xik*Fzk;
                x1=feval(nprob,n,zk1,2);     % 投影得到新的迭代点x_{k+1}
                %NP=NP+1;
                Fk0=Fk;
                Fk=feval(nprob,n,x1,1);
                NF=NF+1;
                nextstep=Step(4);
            case Step(4)
                %% Step 2 更新
                Itr=Itr+1;
                yk=Fk-Fk0;
                switch method
                    case 'JYJ' %投JCAM使用的共轭参数
                        mu=1.4;
                        betak=min(abs(Fk'*(Fk-Fk0))/norm(Fk0)^2,(norm(Fk)^2-norm(Fk)/norm(Fk0)*Fk'*Fk0)/(mu*norm(dk)*norm(Fk)));
                        dk=-1*Fk+betak*dk; 
                    case 'YLS'  %Yu Z, Lin J, Sun J, et al. Spectral gradient projection method for monotone nonlinear equations with convex constraints[J]. Applied numerical mathematics, 2009, 59(10): 2416-2423.
                        r=0.01;
                        sk=x1-x0;
                        yk1=yk+r*sk;
                        betak=norm(sk)^2/(sk'*yk1);
                        dk=-1*betak*Fk;
                    case 'IFR'
                        betak=abs(Fk'*dk)/(-Fk0'*dk)*norm(Fk)^2/norm(Fk0)^2;
                        dk=-1*Fk+betak*dk; 
                end
%                 dk0=dk;
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
