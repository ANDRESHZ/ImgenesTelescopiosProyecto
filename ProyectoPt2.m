clc; clf; clear all; close all;
%% Incializacion 
addpath(genpath('TwIST_v2'));
addpath(genpath('l1magic'));
addpath(genpath('funcionesRegularizar'));
fftnc = @(x) fftshift(fftn(fftshift(x)));
ifftnc = @(x) ifftshift(ifftn(ifftshift(x)));
normalizar=@(imgX)(imgX-min(imgX(:)))/(max(imgX(:))-min(imgX(:)));
vec=@(data) data(:);
n=79;
Circuloinit = zeros(n, n); 
[xp, yp] = meshgrid(1:n, 1:n); 
Circuloinit((xp - n/2).^2 + (yp - n/2).^2 <= (38.5).^2) = 1;
n=79;
% data="TempCorrC1017";
%data="TempCorrC1025";
%data="TempCorrC1032";
data="TempCorrC1044";
load (data)
N=length(C02t)-1;
C02tN=zeros(79,79,N);
C13tN=C02tN;
C20tN=C02tN;
C31tN=C02tN;
for j1 = [1 : N] %Normalizamos todos los frames
    minimo=min(min(C02t{j1}));
    C02tN(:,:,j1)=(C02t{j1}-minimo)/(max(max(C02t{j1}))-minimo);
    minimo=min(min(C13t{j1}));
    C13tN(:,:,j1)=(C13t{j1}-minimo)/(max(max(C13t{j1}))-minimo);
    minimo=min(min(C20t{j1}));
    C20tN(:,:,j1)=(C20t{j1}-minimo)/(max(max(C20t{j1}))-minimo);
    minimo=min(min(C31t{j1}));
    C31tN(:,:,j1)=(C31t{j1}-minimo)/(max(max(C31t{j1}))-minimo);
end
%% Filtrar imagen
xorg=C02tN(:,:,25);%aqui se ajustan la imagenes a observar, si desa una secuancia use un For
x=xorg;
xfilt=imnlmfilt(xorg,'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
x=xfilt;
figure
imshow([xorg x],[],'InitialMagnification',1024)
C02tFILT=zeros(79,79,N);
C13tFILT=zeros(79,79,N);
C20tFILT=zeros(79,79,N);
C31tFILT=zeros(79,79,N);
for im=1:N %% filtramos identicamente a todas la imagenes
    C02tFILT(:,:,im)=imnlmfilt(C02tN(:,:,im),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
    C13tFILT(:,:,im)=imnlmfilt(C13tN(:,:,im),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
    C20tFILT(:,:,im)=imnlmfilt(C20tN(:,:,im),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
    C31tFILT(:,:,im)=imnlmfilt(C31tN(:,:,im),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
end
%% Calculos de PSF basado en la imagen ON-line o Movil.
% imwrite(x,"Evidencias\x1.bmp")
[i, j] = find(ismember(x, max(x(:))));
xC=x(i-2:i+2,j-2:j+2);
minxC=min(xC(:));
PorBus=(1-minxC)/5; %porcentaje de buesqueda 33% podria variar respecto a la iluminacion generla de la imagen.
dist=zeros(1,4);
for iv=1:1:min((78-i),i-1)
    flag=0;
    if(x(i-iv,j)>=minxC*(1-PorBus)) %mirar arriba
     dist(1)=iv;
     flag=1;
    end
    if(x(i+iv,j)>=minxC*(1-PorBus)) %mirar abajo
     dist(2)=iv;
     flag=1;
    end
    if(x(i,j-iv)>=minxC*(1-PorBus)) %izquierda
     dist(3)=iv;
     flag=1;
    end
     if(x(i,j-iv)>=minxC*(1-PorBus)) %derecha
     dist(4)=iv;
     flag=1;
    end
    if (flag==0) %si no existe en ninguna direccion finalizamos
        break;
    end
end

distProm=int8(mean(dist));
PSF=zeros(n,n);
RecorPSF=x(i-distProm:i+distProm,j-distProm:j+distProm);
RecorPSF=(RecorPSF-min(RecorPSF(:)))/(max(RecorPSF(:))-min(RecorPSF(:)));
%RecorPSF=(RecorPSF-min(RecorPSF(:)));
PSF(40-distProm:40+distProm,40-distProm:40+distProm)=RecorPSF;
imshow([PSF],[],'InitialMagnification',1024)
%% Recortar PSF para generar border suaves
CiculoCorte = zeros(n, n); 
[xp, yp] = meshgrid(1:n, 1:n); 
CiculoCorte((xp - n/2).^2 + (yp - n/2).^2 <= (distProm).^2) = 1;
sigma = double((distProm*0.7));
gaussCirc = fspecial('gaussian', 79, 2.3); 
gaussCirc=gaussCirc/max(max(gaussCirc));% normalizar
PSFmod=(gaussCirc.*(CiculoCorte-1)*-1)+(CiculoCorte.*PSF);
figure
imshow([xorg, x,CiculoCorte.*PSF, PSFmod, gaussCirc],[],'InitialMagnification',1024)
PSFFinal= PSF.*CiculoCorte;
%% filtro inverso
ImgGENFilt=xfilt;
ImgGENR=C13tFILT(:,:,20);
clc;
yy = abs(ifftnc(fftnc(ImgGENFilt)./fftnc(PSFmod+0.01)));
yy=normalizar(yy);
figure(9)
imshow([ImgGENFilt, yy],[],'InitialMagnification',1024)
for kkk=0.0001:0.005:0.5
    yyREg=normalizar(deconvreg(ImgGENFilt,PSFmod,kkk));
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyREg,0,Circuloinit,0);
    figure(10)
    imshow([ ImgGENFilt ImgGENR;yy yyREg],[],'InitialMagnification',1024)
    title("K="+kkk+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    pause(0.05)
end

%% Wiener

yyWie2=normalizar(deconvwnr(ImgGENR,PSFmod));
figure
imshow([ImgGENFilt, yyWie2],[],'InitialMagnification',1024)

signal_var = var(ImgGENFilt(:));
NSR = 0.0007 / signal_var;
for kkk=0.001:0.005:0.6
    NSR =kkk;
    yyWie=deconvwnr(ImgGENR,PSFmod,NSR,0.712);
    yyWie=normalizar(yyWie);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyWie,1,Circuloinit,0);
    figure(20)
    imshow([ ImgGENFilt ImgGENR;yyWie2 yyWie],[],'InitialMagnification',1024)
    title("NSR_{input}="+kkk+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    pause(0.05)
end
%% Richardson Lucy
pause(0.5)
for kkk=1:1:30
    yyRL=deconvlucy(ImgGENFilt,PSFmod,kkk,0.0005,Circuloinit,0.002);
    yyRLmod=normalizar(yyRL.*Circuloinit);
    yyRL=deconvlucy(ImgGENFilt,PSFFinal,kkk,0.0005,Circuloinit,0.00002);
    yyRLfinal=normalizar(yyRL.*Circuloinit);
    yyRL=deconvlucy(ImgGENR,PSFmod,kkk,0.0005,Circuloinit,0.00002);
    yyRLgauss=normalizar(yyRL.*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyRLmod,1,Circuloinit,0);
    figure(21)
    imshow([ImgGENFilt, yyRLmod;yyRLfinal,yyRLgauss],[],'InitialMagnification',1024)
    title("#_{iter}="+kkk+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    pause(0.05)
end
%% SVD
g = vec(fftnc(ImgGENFilt));
time0 = clock;
% gruido = vec(fftnc(Iruidoc));
A = diag(vec(fftnc(PSFmod)));
[U,s,V] = csvd(A);
%% Picard Graph
g = vec(fftnc(ImgGENFilt));
sigma = s;
figure()
semilogy(sigma(1:end),'b.')
hold on 
valse=abs(U'*g);
semilogy(valse(1:end),'go')
y = abs(U'*g)./sigma;
semilogy(y(1:end),'yx')
[~, idx] = min(y(1:6241));
plot(idx,y(idx),"rd",'MarkerSize',10,'LineWidth',2)
hold off
grid("on")
xlabel("k")
legend("$\sigma_k$","$|U^{\mathrm{T}}g|$","$\frac{|U^{\mathrm{T}}g|}{\sigma_k}$", "Interpreter","latex")
title("valor $k_{opt}=$"+idx,"Interpreter","latex")
%% DP Principio de Discrepancia 
g = vec(fftnc(ImgGENFilt));
g2 = vec(fftnc(ImgGENR));
for delta=[10:3:44]
time0 = clock;
[x_DP,lambda] = discrep(U,s,V,g,delta);
[x_DP2,lambda2] = discrep(U,s,V,g2,delta);
x_ka=ifftnc(reshape(x_DP,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    x_ka2=normalizar(abs(ifftnc(reshape(x_DP2,79,79))).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(35)
    imshow([ImgGENFilt ImgGENR; x_ka x_ka2],[],"InitialMagnification",400)
    title("\Delta_{DPF}="+lambda+" \Delta_{DPR}="+lambda2+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
    pause(0.005)
end

%% K TSVD
g = vec(fftnc(ImgGENFilt));
g2 = vec(fftnc(ImgGENR));
for K=[1200:100:4000]
    time0 = clock;
    [x_k,rho,eta] = tsvd(U,s,V,g,K);
    [x_k2,rho,eta] = tsvd(U,s,V,g2,K);
    x_ka=ifftnc(reshape(x_k,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    x_ka2=ifftnc(reshape(x_k2,79,79));
    x_ka2=normalizar(abs(x_ka2).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(36)
    imshow([ImgGENFilt ImgGENR; x_ka x_ka2],[],"InitialMagnification",400)
     title("K_{TSVD}="+K+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
     pause(0.005)
end

%% L1 Curve TSVD
g = vec(fftnc(ImgGENFilt));
g2 = vec(fftnc(ImgGENR));
figure
time0 = clock;
[reg_corner,rho,eta,reg_param]=l_curve(U,s,g,'tsvd');
figure
[reg_corner2,rho2,eta2,reg_param2]=l_curve(U,s,g2,'tsvd');
    [x_kLcurv,rho,eta] = tsvd(U,s,V,g,reg_corner);
    [x_kLcurv2,rho,eta] = tsvd(U,s,V,g2,reg_corner2);
    x_ka=ifftnc(reshape(x_kLcurv,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    x_ka2=normalizar(abs(ifftnc(reshape(x_kLcurv2,79,79))).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    [SNR_ind2,BRISQUE2,NIQE2,PIQE2,NRMSE2,SNR2,PSNR2,SSIM2]=ErroresIMG(x_ka2,1,Circuloinit,0);
    figure(37)
    imshow([ImgGENFilt ImgGENR; x_ka x_ka2],[],"InitialMagnification",400)
    title({["TSVD K_{LCurveF}="+reg_corner+"=>|SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE];...
        ["TSVD K_{LCurveR}="+reg_corner2+"=>|SNR_{ind}=="+SNR_ind2+" |BRISQUE="+BRISQUE2+" |PSNR="+NIQE2+" |SSIM="+PIQE2]});
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));

%% GCV TSVD Generalized cross-validation 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_min,G,reg_param] = gcv(U,s,g,'tsvd');
[x_kGCV,rho,eta] = tsvd(U,s,V,g,reg_min);
    x_ka=ifftnc(reshape(x_kGCV,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(38)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("K_{GCV}="+reg_min+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% NCP TSVD Normalized cumulative periodogram
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_minNCP,G,reg_param] = ncp2(U,s,g,'tsvd');
[x_kNCP,rho,eta] = tsvd(U,s,V,g,reg_minNCP);
    x_ka=ifftnc(reshape(x_kNCP,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(329)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
     title("K_{NCP}="+reg_minNCP+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% QOC Quasi-optimality criterion 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_minQOC,Q,reg_param] = quasiopt(U,s,g,'tsvd');
[x_kQOC,rho,eta] = tsvd(U,s,V,g,reg_minQOC);
    x_ka=ifftnc(reshape(x_kQOC,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(352)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("K_{QOC}="+reg_minQOC+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));

%% TIKONOV
g = vec(fftnc(ImgGENFilt));
g2 = vec(fftnc(ImgGENR));
for lambda=(0.11:0.01:0.8)
    time0 = clock;
[x_lambda,rho,eta] = tikhonov(U,s,V,g,lambda);
[x_lambda2,rho2,eta2] = tikhonov(U,s,V,g2,lambda);
x_tik=ifftnc(reshape(x_lambda,79,79));
x_tik=normalizar(abs(x_tik).*Circuloinit);
x_tik2=normalizar(abs(ifftnc(reshape(x_lambda2,79,79))).*Circuloinit);
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_tik,1,Circuloinit,0);
figure(402)
imshow([ ImgGENFilt ImgGENR; x_tik x_tik2],[],"InitialMagnification",400)
title("\lambda_{TIKHONOV}="+lambda+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
     pause(0.005)
end
%% L1 Curve Tikhonov
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[lambdaTK,rho,eta,reg_param]=l_curve(U,s,g,'Tikh');
[x_lambdaLCURV,rho,eta] = tikhonov(U,s,V,g,lambdaTK);
    x_ka=ifftnc(reshape(x_lambdaLCURV,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(422)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{Lcurve}="+lambdaTK+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% GCV Tikhonov Generalized cross-validation 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintk,G,reg_param] = gcv(U,s,g,'tikh');
[x_tkGCV,rho,eta] = tikhonov(U,s,V,g,reg_mintk);
    x_ka=ifftnc(reshape(x_tkGCV,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(400)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{GCV}="+reg_mintk+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% NCP Tikhonov Normalized cumulative periodogram
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintkNCP,G,reg_param] = ncp2(U,s,g,'tikh');
[x_tkNCP,rho,eta] = tikhonov(U,s,V,g,reg_mintkNCP);
    x_ka=ifftnc(reshape(x_tkNCP,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(403)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{NCP}="+reg_mintkNCP+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% QOC Quasi-optimality criterion 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintkQOC,Q,reg_param] = quasiopt(U,s,g,'tikh');
[x_tkQOC,rho,eta] = tikhonov(U,s,V,g,reg_mintkQOC);
    x_ka=ifftnc(reshape(x_tkQOC,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,1,Circuloinit,0);
    figure(404)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{QOC}="+reg_mintkQOC+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));


%% L1 MAgic
x=ImgGENFilt;
largescale = 1;
n = 78;
II = x(1:0+n,1:0+n);
N = n*n;
I = II/norm(II(:));
%I = I - mean(I(:));
x = reshape(I,N,1);
% num obs
K = 6000;
% permutation P and observation set OMEGA
P = randperm(N)';
q = randperm(N/2-1)+1;
OMEGA = q(1:K/2)';
% measurement matrix
if (largescale)
  A = @(z) A_f(z, OMEGA, P);
  At = @(z) At_f(z, N, OMEGA, P);
  % obsevations
  b = A(x);
  % initial point
  x0 = At(b);
else
  FT = 1/sqrt(N)*fft(eye(N));
  A = sqrt(2)*[real(FT(OMEGA,:)); imag(FT(OMEGA,:))];
  A = [1/sqrt(N)*ones(1,N); A];
  At = [];
  % observations
  b = A*x;
  % initial point
  x0 = A'*b;
end
figure
imshow(II/max(max(II)))

IMGh=II/max(max(II));
HH=0;
for j=[2,3,6,8:1:11]%iteraciones=(multiplos de 4)-1 
    HH=HH+1
    epsilon =(4+(j*2))*3^(-j)
       
    tvI = sum(sum(sqrt([diff(I,1,2) zeros(n,1)].^2 + [diff(I,1,1); zeros(1,n)].^2 )));
    disp(sprintf('Original TV = %.3f', tvI));
    time0 = clock;
    xp =  tvqc_logbarrier(x0, A, At, b, epsilon, 1e-4,3, 1e-8, 100);
    Ip = reshape(xp, n, n);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
    tam=size(IMGh);
    if tam(2)==n*4
        if HH==4
            IMG=IMGh;
        else
            IMG=[IMG;IMGh];
        end
        
        IMGh=Ip/max(max(Ip));
    else
        IMGh=[IMGh,Ip/max(max(Ip))];
    end   
end
IMG=[IMG;IMGh];
figure
imshow(imresize(IMG,3,'box'));
title("L1 Magic TV \epsilon=["+(4+(j*2))*3^(-2)+" al "+ (4+(j*2))*3^(-11)+"]")

%% Twist Manual
x=ImgGENFilt;
alpha = 0.5;
beta = 0.25;
iterTWIST=300;
iteraTV=5;
lambdas=[0.86:-0.02:0.74];
lambdas = (4+(j*2))*2.^(-[3+4/7:1/7:(8-17/7)]);
for im=15%N-10
x=ImgGENFilt;
y=x;
IMGh=[y];
HH=0;
for j=lambdas
    HH=HH+1;
%     lambda = (4+(j*2))*2^(-j)
    lambda=j;
%     time0 = clock;
    x_twist = TWIST_manual(x,y,alpha,beta,iterTWIST,lambda,iteraTV);
%     disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
    %mostrar imagenes
    tam=size(IMGh);
    if tam(2)==length(y)*4
        if HH==4
            IMG=IMGh;
        else
            IMG=[IMG;IMGh];
        end
        
        IMGh=x_twist/max(max(x_twist));
    else
        IMGh=[IMGh,x_twist/max(max(x_twist))];
    end   
end
IMG=[IMG;IMGh];
figure(503);
imshow(imresize(IMG,3,'box'));
title("TwIST L pesado  \lambdas del "+lambdas(1)+" al "+ lambdas(end))
end


%% TwIST CON TV
x1=ImgGENR;
N=length(x1);
B=PSFmod;
B=fftshift(B);%circular y centrar
B=B/sum(sum(B));%normalize
y = ImgGENFilt;
% y = x_twist222;
% y=xoptTwIST;
% y=C02tFILT(:,:,23);
% y=xorg;
K=fft2(B);
KC=conj(K);
% operadores de convolución
%  convolution operators
A = @(x) real(ifft2(K.*fft2(x)));
AT = @(x) real(ifft2(KC.*fft2(x)));
tv_iters = 5;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
% TV regularizer;
Phi = @(x) TVnorm(x);
%Phi = @(x) sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));
varx = var(y(:)); 	% approximate var of x
x0 = real(ifft2(KC./(abs(KC).^2+10*0.0002^2/varx).*fft2(y)));
x02 = real(ifft2(KC./(abs(KC).^2+10*0.0002^2/varx).*fft2(x1)));
x00=normalizar(x0);
%%
for taus=[0.000001:0.000025:0.0005]
    tau = taus;
    lam1=1e-10;    
    tolA = 1e-8;
    [x_twistTV]=TwIST(y,A,tau,'AT', AT,'lambda',lam1,'True_x', x1,'Psi',...
             Psi,'Phi',Phi,'Monotone',1,'Initialization',x0,'StopCriterion',1,'ToleranceA',tolA,'Verbose', 0);
     [x_twistTV2]=TwIST(x1,A,tau,'AT', AT,'lambda',lam1,'True_x', y,'Psi',...
         Psi,'Phi',Phi,'Monotone',1,'Initialization',x02,'StopCriterion',1,'ToleranceA',tolA,'Verbose', 0);
    x_twistTV=normalizar(x_twistTV);
    x_twistTV2=normalizar(x_twistTV2);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_twistTV,x1,Circuloinit,0);
    figure(505);
    imshow([x1 y; x_twistTV2  x_twistTV],[],'InitialMagnification',1024);
    title("\tau="+tau+" |SNR_{ind}="+SNR_ind+" |BRISQUE="+BRISQUE+" |NIQE="+NIQE+" |PIQE="+PIQE);
    pause(0.005)
end
