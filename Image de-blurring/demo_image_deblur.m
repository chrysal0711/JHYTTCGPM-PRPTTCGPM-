% This demo shows the reconstruction of a sparse image
% of randomly placed plus an minus ones
% 
clear all
close all
clf
randn('seed',1); % 1 for the experiments in the paper
rand('seed',1); % 1 for the experiments in the paper
addpath Images
% the test images set
A={ 'lena.bmp' 'boat.bmp' 'bridge.bmp' 'flower.bmp' 'rose.bmp'...
   'flower.bmp' 'fourviere.bmp'};
% f = double(imread('cameraman.png'));
% f = double(imread('lena.png'));
% f = double(imread('barbara.png'));
% f = double(imread('man.bmp'));
f=double(imread(A{2})); 
[m,n] = size(f);
scrsz = get(0,'ScreenSize');
figure(1)
%set(1,'Position',[10 scrsz(4)*0.05 scrsz(3)/4 0.85*scrsz(4)])
% subplot(1,2,1)
imagesc(f)
colormap(gray(255))
axis off
axis equal
title('Original image','FontName','Times','FontSize',22)
% title(sprintf('Original: %4d  %4d',m,n),'FontSize',22)
% create observation operator; in this case 
% it will be a blur function composed with an
% inverse weavelet transform
disp('Creating observation operator...');

middle = n/2 + 1;

% uncomment the following lines for Experiment 1 (see paper)
% sigma = 0.56;
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= 1; 
%    end
% end


% uncomment the following lines for Experiment 2 (see paper)
% sigma = sqrt(2);
sigma = sqrt(0.01);
h = zeros(size(f)); % h is the same size as f and all zeros.
for i=-4:4
   for j=-4:4
      h(i+middle,j+middle)= (1/(1+i*i+j*j));
   end
end

% uncomment the following lines for Experiment 3 (see paper)
% sigma = sqrt(8);
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= (1/(1+i*i+j*j));
%    end
% end

% % center and normalize the blur
h = fftshift(h);  % ��Ƶ��ͼ��λ����Ƶ����Ƶ��ͼ����
% fftshift is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
h = h/sum(h(:));  % sum(h(:)): sum all elements of h.

% definde the function handles that compute 
% the blur and the conjugate blur.
R = @(x) real(ifft2(fft2(h).*fft2(x))); % fft2(X) returns the two-dimensional Fourier transform of matrix X.
% ifft2(F) returns the two-dimensional inverse Fourier transform of matrix F
RT = @(x) real(ifft2(conj(fft2(h)).*fft2(x)));


% define the function handles that compute 
% the products by W (inverse DWT) and W' (DWT)
wav = daubcqf(2);
W = @(x) midwt(x,wav,3);
WT = @(x) mdwt(x,wav,3);

%Finally define the function handles that compute 
% the products by A = RW  and A' =W'*R' 
A = @(x) R(W(x));
AT = @(x) WT(RT(x));

fid=fopen('mytext.txt','w');
% generate noisy blurred observations
y = R(f) + sigma*randn(size(f));
figure(2)
% subplot(1,2,2)
imagesc(y)
colormap(gray(255))
axis off
axis equal
title('Blurred image','FontName','Times','FontSize',22)
% title(['Blurred: ',num2str(blur_label),'dB'],'FontName','Times','FontSize',22)

% regularization parameter
tau = .35;

% set tolA
tolA = 1.e-5;

% Run CGD until the relative change in objective function is no
% larger than tolA
disp('Starting CGDESCENT_CS_image ...')
[x_CGDCS,theta_debias,obj_CGDCS,t_CGDCS,debias_CGDCS,mses_CGDCS]= ...
	CGDESCENT_CS_image(y,A,tau,...
	'Debias',0,...
	'AT',AT,... 
    'True_x',WT(f),...
	'Initialization',AT(y),... %0,...%
	'StopCriterion',1,...
	'ToleranceA',tolA,...
     'Verbose',0);
 SNR_CGDCS = 20*log10(norm(f,'fro')/norm(f-W(x_CGDCS),'fro'));
 Itr_CGDCS = length(t_CGDCS);
 fprintf(1,'CGDCS: Iter=%4d, obj=%10.3e, times=%6.2f, MSE=%10.4e, SNR=%6.2f\n', Itr_CGDCS,...
      obj_CGDCS(end),t_CGDCS(end),mses_CGDCS(end)/m/n,SNR_CGDCS);
%  CGD_label=[roundn(t_CGDCS(end)) roundn(20*log10(norm(f,'fro')/norm(f-W(x_CGDCS),'fro')))];
  
% Run HTTCGP until the relative change in objective function is no
% larger than tolA
disp('Starting ATTCGP_CS_image ...')
[x_ATTCGPCS,theta_debias,obj_ATTCGPCS,t_ATTCGPCS,debias_ATTCGPCS,mses_ATTCGPCS]= ...
	ATTCGP_CS(y,A,tau,...
	'Debias',0,...
	'AT',AT,... 
    'True_x',WT(f),...
	'Initialization',AT(y),... %0,...%
	'StopCriterion',1,...
	'ToleranceA',tolA,...
     'Verbose',0);
 SNR_ATTCGPCS = 20*log10(norm(f,'fro')/norm(f-W(x_ATTCGPCS),'fro'));
 Itr_ATTCGPCS = length(t_ATTCGPCS);
 fprintf(1,'ATTCGPCS: Iter=%4d, obj=%10.3e, times=%6.2f, MSE=%10.4e, SNR=%6.2f\n', Itr_ATTCGPCS,...
      obj_ATTCGPCS(end),t_ATTCGPCS(end),mses_ATTCGPCS(end)/m/n,SNR_ATTCGPCS);
%   ATTCGP_label=[roundn(t_ATTCGPCS(end)) roundn(20*log10(norm(f,'fro')/norm(f-W(x_ATTCGPCS),'fro')))];
  
% Run ATTCGP until the relative change in objective function is no
% larger than tolA
disp('Starting PRPTTCGPM_CS_image ...')
[x_PRPTTCGPMCS,theta_debias,obj_PRPTTCGPMCS,t_PRPTTCGPMCS,debias_PRPTTCGPMCS,mses_PRPTTCGPMCS]= ...
	PRPTTCGPM_CS(y,A,tau,...
	'Debias',0,...
	'AT',AT,... 
    'True_x',WT(f),...
	'Initialization',AT(y),... %0,...%
	'StopCriterion',1,...
	'ToleranceA',tolA,...
     'Verbose',0);
 SNR_PRPTTCGPMCS=20*log10(norm(f,'fro')/norm(f-W(x_PRPTTCGPMCS),'fro'));
 Itr_PRPTTCGPMCS = length(t_PRPTTCGPMCS);
 fprintf(1,'JHYTTCGPMCS: Iter=%4d, obj=%10.3e, times=%6.2f, MSE=%10.4e, SNR=%6.2f\n', Itr_PRPTTCGPMCS,...
      obj_PRPTTCGPMCS(end),t_PRPTTCGPMCS(end),mses_PRPTTCGPMCS(end)/m/n,SNR_PRPTTCGPMCS);
%  ATTCGP_label=[roundn(t_ATTCGPCS(end)) roundn(20*log10(norm(f,'fro')/norm(f-W(x_ATTCGPCS),'fro')))];


% disp('Starting NDYSPCGP_CS_image ...')
% [x_NDYSPCGPCS,theta_debias,obj_NDYSPCGPCS,t_NDYSPCGPCS,debias_NDYSPCGPCS,mses_NDYSPCGPCS]= ...
% 	NDYSPCGP_CS(y,A,tau,...
% 	'Debias',0,...
% 	'AT',AT,... 
%     'True_x',WT(f),...
% 	'Initialization',AT(y),... %0,...%
% 	'StopCriterion',1,...
% 	'ToleranceA',tolA,...
%      'Verbose',0);
%  SNR_NDYSPCGPCS=20*log10(norm(f,'fro')/norm(f-W(x_NDYSPCGPCS),'fro'));
%  Itr_NDYSPCGPCS = length(t_NDYSPCGPCS);
%  fprintf(1,'NDYSPCGPCS: Iter=%4d, obj=%10.3e, times=%6.2f, MSE=%10.4e, SNR=%6.2f\n', Itr_NDYSPCGPCS,...
%       obj_NDYSPCGPCS(end),t_NDYSPCGPCS(end),mses_NDYSPCGPCS(end)/m/n,SNR_NDYSPCGPCS);

psnr_CGDCS = mpsnr(f,W(x_CGDCS));
psnr_ATTCGPCS = mpsnr(f,W(x_ATTCGPCS));
psnr_PRPTTCGPMCS = mpsnr(f,W(x_PRPTTCGPMCS));
%psnr_NDYSPCGPCS = mpsnr(f,W(x_NDYSPCGPCS));
mssim_CGDCS = ssim(f, W(x_CGDCS));       % [mssim, ssim_map] = ssim(f, g);
mssim_ATTCGPCS = ssim(f, W(x_ATTCGPCS));
mssim_PRPTTCGPMCS = ssim(f, W(x_PRPTTCGPMCS));
%mssim_NDYSPCGPCS = ssim(f, W(x_NDYSPCGPCS));

fprintf(fid,'%d/%.2f/%.2f/%.2f & %d/%.2f/%.2f/%.2f & %d/%.2f/%.2f/%.2f\\\\\n', ... 
        Itr_CGDCS,t_CGDCS(end),psnr_CGDCS,mssim_CGDCS,Itr_ATTCGPCS,t_ATTCGPCS(end), ...
        psnr_ATTCGPCS,mssim_ATTCGPCS,Itr_PRPTTCGPMCS,t_PRPTTCGPMCS(end),...
        psnr_PRPTTCGPMCS,mssim_PRPTTCGPMCS);
fclose(fid);
  
% %% ================= Plotting results ==========
figure(3)
% subplot(1,5,3)
% if prod(size(theta_debias))~=0 % prod(size(theta_debias)) is the product of the elements of the vector (size(theta_debias)
%    imagesc(W(theta_debias))
% else
   imagesc(W(x_CGDCS))
% end
colormap(gray)
axis off
axis equal
title('CGD','FontName','Times','FontSize',22)

figure(4)
% subplot(1,5,4)
imagesc(W(x_ATTCGPCS))
colormap(gray)
axis off
axis equal
title('ATTCGP','FontName','Times','FontSize',22)

% 
figure(5)
% subplot(1,5,5)
imagesc(W(x_PRPTTCGPMCS))
colormap(gray)
axis off
axis equal
title('JHYTTCGPM4','FontName','Times','FontSize',22)

% figure(6)
% imagesc(W(x_NDYSPCGPCS))
% colormap(gray)
% axis off
% axis equal
% title('NDYSPCGP','FontName','Times','FontSize',22)
% 
% figure(11)
% scrsz = get(0,'ScreenSize');
% lft = 0.55*scrsz(3)-10;
% btm = 0.525*scrsz(4);
% wdt = 0.45*scrsz(3);
% hgt = 0.375*scrsz(4);
% 
% set(2,'Position',[lft btm wdt hgt])
% plot(obj_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(obj_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(obj_SGCS,'g:','LineWidth',1.8);
% plot(obj_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','SGCS','IST');
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Objective function','FontName','Times','FontSize',16)
% xlabel('Iterations','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times');
% set(leg,'FontSize',16);
% set(gca,'FontName','Times');
% set(gca,'FontSize',16);
%     
% 
% figure(12)
% scrsz = get(0,'ScreenSize');
% lft = 0.55*scrsz(3)-10;
% btm = 0.025*scrsz(4);
% wdt = 0.45*scrsz(3);
% hgt = 0.375*scrsz(4);
%  
% set(3,'Position',[lft btm wdt hgt])
% plot(times_QP_BB_mono,obj_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(times_QP_BB_notmono,obj_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(t_SGCS,obj_SGCS,'g:','LineWidth',1.8);
% plot(times_IST,obj_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','SGCS','IST');
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Objective function','FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times');
% set(leg,'FontSize',16);
% set(gca,'FontName','Times');
% set(gca,'FontSize',16);
%     
% 
% figure(13)
% 
% plot(times_QP_BB_mono,mses_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(times_QP_BB_notmono,mses_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(t_SGCS,mses_SGCS,'g:','LineWidth',1.8);
% plot(times_IST,mses_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','SGCS','IST');
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Deconvolution MSE','FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times');
% set(leg,'FontSize',16);
% set(gca,'FontName','Times');
% set(gca,'FontSize',16);

% fprintf(1,'\nSGCS monotone; cpu: %6.2f secs (%d iterations)\n',...
%         t_SGCS,length(t_SGCS))
% fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
%           obj_SGCS(end),(1/(n*m))*norm(x_SGCS-f)^2)      



