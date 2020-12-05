% function demo_l2_TV
% TV based image restoration using TwIST
% cameraman, blur uniform 9*9, BSNR = 40 dB
%
% This demo illustrates the computation of  
%
%     xe = arg min 0.5 ||Ax-y||^2 + tau TV(x)
%             x
% 
% with the TwIST algorithm. 
%function demo_l2_TV
% TV based image restoration using TwIST
% cameraman, blur uniform 9*9, BSNR = 40 dB
%
% This demo illustrates the computation of  
%
%     xe = arg min 0.5 ||Ax-y||^2 + tau TV(x)
%             x
% 
% with the TwIST algorithm. 
%
% The proximal operator, i.e., the solution of
%
%     xe = arg min 0.5 ||x-y||^2 + tau TV(x)
%             x
% is  given the Chambolle algorithm
% 
% A. Chambolle, “An algorithm for total variation minimization and
% applications,” Journal of Mathematical Imaging and Vision, vol. 20,
% pp. 89-97, 2004.
%
%
% For further details about the TwIST algorithm, see the paper:
%
% J. Bioucas-Dias and M. Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image 
% Restoration",  IEEE Transactions on Image processing, 2007.
% 
% and
% 
% J. Bioucas-Dias and M. Figueiredo, "A Monotonic Two-Step 
% Algorithm for Compressive Sensing and Other Ill-Posed 
% Inverse Problems", submitted, 2007.

%
% Paper: J. Bioucas- Dias and M.  Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image Restoration", 
% Submitted to IEEE Transactions on Image processing,  2007.
%
% (available at   http://www.lx.it.pt/~bioucas/publications.html)
%
%
% 
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior Técnico, Outubro, 2007

% clear all
% close all



%%%%%%%%%%%%%%%%%%%%% Original Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = double(imread('cameraman.tif'));
%x = double(imread('phantom.tif'));
%x=256*phantom(256);
x=ImgGenerada.*Circuloinit;

N=length(x);
colormap(gray);

% remove mean (not necessary)
% mu=mean(x(:));
% x=x-mu;

% define the blur operator
% N must be even
middle=N/2+1;
% blurr matrix
B=zeros(N);
% Uniform blur
lx=4; %blur x-size
ly=4; %blurr y-size
B((middle-ly):(middle+ly),(middle-lx):(middle+lx))=1;
B=PSFmod;
%circularly center
B=fftshift(B);
%normalize
B=B/sum(sum(B));

% convolve
y = ImgGENFilt;
% y = x_twist222;
% y=xoptTwIST;
% y=C02tFILT(:,:,23);
% y=xorg;

% figure(1); colormap gray; 
% imshow(x,[],'InitialMagnification',1024); axis off;
% title('Original image')
% figure(2); colormap gray; 
% imshow(y,[],'InitialMagnification',1024); axis off;
% title('Noisy and blurred image')
% smoothing parameter (empirical setting)

K=fft2(B);
KC=conj(K);

% handle functions for TwIST
%  convolution operators
A = @(x) real(ifft2(K.*fft2(x)));
AT = @(x) real(ifft2(KC.*fft2(x)));

% denoising function;
tv_iters = 5;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
% TV regularizer;
Phi = @(x) TVnorm(x);
%Phi = @(x) sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));
varx = var(y(:)); 	% approximate var of x
x0 = real(ifft2(KC./(abs(KC).^2+10*0.0002^2/varx).*fft2(y)));
for taus=[0.000001:0.00005:0.0005]
tau = taus;%2e-2*0.2^2/0.56^2;
% extreme eigenvalues (TwIST parameter)
lam1=1e-10;    
% TwIST is not very sensitive to this parameter
% rule of thumb: lam1=1e-4 for severyly ill-conditioned% problems
%              : lam1=1e-1 for mildly  ill-conditioned% problems
%              : lam1=1    when A = Unitary matrix
% ------------  TV experiments ---------------------
% start with the wiener filter


tolA = 1e-8;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST(y,A,tau,...
         'AT', AT, ...
         'lambda',lam1,...
         'True_x', x,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'Initialization',x0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 0);

% figure(3); colormap gray; 
% imshow(x_twist,[],'InitialMagnification',1024); axis off;
% drawnow;
% title('TwIST restored image')

figure(4);
imshow([x y; normalizar(x0) normalizar(x_twist)],[],'InitialMagnification',1024);
title("x y x0 TwIST restored image"+tau)
end

% -- IST (lam1=1) ---------------------------
% stop criterium:  the objective function 
% falls below obj_twist(end)
%
% IST takes too long and thus we run only 200 iterations
[x_ist,dummy,obj_ist,...
  times_ist,dummy,mse_ist]= ...
         TwIST(y,A,tau,...
         'AT', AT, ...
         'lambda',1,...
         'True_x', x,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'MaxiterA',200,...
         'Initialization',x0,...
         'StopCriterion',3,...
       	 'ToleranceA',obj_twist(end),...
         'Verbose', 0);
figure(5);
imshow([x y; normalizar(x0) normalizar(x_ist)],[],'InitialMagnification',1024);
title("x y TwIST IST restored image"+tau)


% figure(5)
% subplot(2,1,1)
% semilogy(times_twist,obj_twist,'r',...
%          times_ist,obj_ist,'b','LineWidth',2)
% legend('TwIST','IST')
% st=sprintf('tau = %2.2e, sigma = %2.2e',tau,sigma),...
% title(st)
% ylabel('Obj. function')
% xlabel('CPU time (sec)')
% 
% grid
% subplot(2,1,2)
% plot(times_twist(2:end),mse_twist(2:end),'r',...
%          times_ist(2:end),mse_ist(2:end),'b','LineWidth',2)
% legend('TwIST','IST')
% ylabel('MSE')
% xlabel('CPU time (sec)')


fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));
fprintf(1,'IST     CPU time - %f\n', times_ist(end));
%% 
figure
imshow([x,y;normalizar(x_twist), normalizar(x_ist)],[],'InitialMagnification',1024);
