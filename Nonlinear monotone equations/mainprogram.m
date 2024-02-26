% Matlab Model by Jianghua Yin (Dec.,2019,Yulin)
% Copyright (C) 2019 Jian Group
% All Rights Reserved
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that the copyright notice and this
% permission notice appear in all supporting documentation.  

% This is a demo of JHYTTCGPM(which respond to the method PRPTTCGPM1-4) method for solving the following constrained
% nonlinear monotone equations
% F(x)=0, x\in C
% where C is a nonempty closed convex set.

% using the algorithm modified three-term conjugate gradient projection method, described in the following paper
%
%  A conjugate gradient projection method with restart procedure for solving constraint equations and image restorations
% -----------------------------------------------------------------------
% Copyright (2019): Jianghua Yin
% ----------------------------------------------------------------------
%
% The first version of this code by Jianghua Yin, Dec., 16, 2019
clear all 
clc;
format long
format compact
%% 建立latex表格文件
fid_tex=fopen('mytext.txt','w'); 
%% 初始化参数
%np=168;    %问题集的个数
np=189;    %问题集的个数
%168+21 = 189
ns=6;     %算法的个数
T=zeros(np,ns);F=zeros(np,ns);N=zeros(np,ns); %G=zeros(np,ns);
%% 运行
for i=1:np
  % 1:np
[name,n,a]=init(i);
i
% name=getname(nprob);
[T1,NFF1,NI1,G1] = line1(i,'ATTCGP');     % acceleration
%[T2,NFF2,NI2,G2] = line2(i,'PDY');        % Liu J, Feng Y. (2018) PDY, Numerical Algorithms
[T2,NFF2,NI2,G2] = SGPYJJ(i,'HTTCG');        % Liu J, Feng Y. (2018) PDY, Numerical Algorithms
[T3,NFF3,NI3,G3] = SGPYHHM1(i,'PRPTTCGPM1');      % Sun M., Liu J., (2016) HCGP, Calcolo
[T4,NFF4,NI4,G4] = SGPYHHM2(i,'PRPTTCGPM2');     % Gao P.T., He C.J., Liu Y., (2019) ATTCGP, Applied Mathematics and Computation
[T5,NFF5,NI5,G5] = SGPYHHM3(i,'PRPTTCGPM3');    % Sun M., Tian M.Y., (2019)
[T6,NFF6,NI6,G6] = SGPYHHM4(i,'PRPTTCGPM4');

T(i,:) = [T1,T2,T3,T4,T5,T6];
F(i,:) = [NFF1,NFF2,NFF3,NFF4,NFF5,NFF6];
N(i,:) = [NI1,NI2,NI3,NI4,NI5,NI6];
% G(i,:) = [G1,G2,G3,G4]; %,G5];
b=n/10000;
%% 输入到 latex 
fprintf(fid_tex,'%s $x_{%d}(%d)$ &%d/%d/%.3f/%.2e&%d/%d/%.3f/%.2e&%d/%d/%.3f/%.2e&%d/%d/%.3f/%.2e&%d/%d/%.3f/%.2e&%d/%d/%.3f/%.2e \\\n',...
                name,b,b,NI1,NFF1,T1,G1,NI2,NFF2,T2,G2,NI3,NFF3,T3,G3,NI4,NFF4,T4,G4,NI5,NFF5,T5,G5,NI6,NFF6,T6,G6);%,NI7,NFF7,NGG7,T7);
end
%% 关闭文件
fclose(fid_tex);

%% 画图
clf;   %clf删除当前图形窗口中、
       %%句柄未被隐藏(即它们的HandleVisibility属性为on)的图形对象。
figure(1);
%subplot(2,2,1);
perf(T,'logplot');
%title('时间性能');
%set(gca,'ylim',[0.3,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('ATTCGP' ,'HTTCG','JHYTTCGPM1','JHYTTCGPM2','JHYTTCGPM3','JHYTTCGPM4');
%subplot(2,2,2);
figure(2);
perf(F,'logplot');
%title('目标函数计算性能');
% set(gca,'ylim',[0.1,1]);
xlabel('\tau','Interpreter','tex');                 
ylabel('\rho(\tau)','Interpreter','tex');            
legend('ATTCGP' ,'HTTCG','JHYTTCGPM1','JHYTTCGPM2','JHYTTCGPM3','JHYTTCGPM4');
%subplot(2,2,3);
figure(3);
perf(N,'logplot');
%title('迭代次数性能');
%set(gca,'ylim',[0.5,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('ATTCGP' ,'HTTCG','JHYTTCGPM1','JHYTTCGPM2','JHYTTCGPM3','JHYTTCGPM4'); 