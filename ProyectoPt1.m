% clc; clf; clear all; close all;
%% Incializacion 
addpath(genpath('TwIST_v2'));
addpath(genpath('l1magic'));
addpath(genpath('funcionesRegularizar'));
fftnc = @(x) fftshift(fftn(fftshift(x)));
ifftnc = @(x) ifftshift(ifftn(ifftshift(x)));
vec=@(data) data(:);
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
for im=1:N %% filtramos identicamente a todas la imagenes
    C02tFILT(:,:,im)=imnlmfilt(C02tN(:,:,im),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
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
%% MIRANDO INFO
n=79;
N=30;
dataold=data;
RUIDOS=zeros(79,79,2*N);
RUIDOS2=zeros(79,79,2*N);
RUIDOS3=zeros(79,79,2*N);
RUIDOS4=zeros(79,79,2*N);
Mascara=zeros(79,79);
Mascara2=zeros(79,79);
Mascara3=zeros(79,79);
Mascara4=zeros(79,79);
ruidoAVG=zeros(79,79);
ruidoAVG2=zeros(79,79);
ruidoAVG3=zeros(79,79);
ruidoAVG4=zeros(79,79);
for d=1:1:2
    if d==1
       %data="TempCorrC1017";
       data="TempCorrC1025";
       %data="TempCorrC1032";
       %data="TempCorrC1044";
    else d==2
        data="TempCorrC1032";
    end
    load (data);
    for im=1:N
        MINIMAL=min(min(C02t{im}));
        MAXIMUN=max(max(C02t{im}));
        C02tNN=(C02t{im}-MINIMAL)/(MAXIMUN-MINIMAL);
        Mascara=((C02tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*1.0)+Mascara;
        RUIDOS(:,:,im+(d*N))=(C02tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*C02tNN;
        MINIMAL=min(min(C13t{im}));
        MAXIMUN=max(max(C13t{im}));
        C13tNN=(C13t{im}-MINIMAL)/(MAXIMUN-MINIMAL);
        Mascara2=((C13tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*1.0)+Mascara2;
        RUIDOS2(:,:,im+(d*N))=(C13tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*C13tNN;
        MINIMAL=min(min(C20t{im}));
        MAXIMUN=max(max(C20t{im}));
        C20tNN=(C20t{im}-MINIMAL)/(MAXIMUN-MINIMAL);
        Mascara3=((C20tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*1.0)+Mascara3;
        RUIDOS3(:,:,im+(d*N))=(C20tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*C20tNN;
        MINIMAL=min(min(C31t{im}));
        MAXIMUN=max(max(C31t{im}));
        C31tNN=(C31t{im}-MINIMAL)/(MAXIMUN-MINIMAL);
        Mascara4=((C31tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*1.0)+Mascara4;
        RUIDOS4(:,:,im+(d*N))=(C31tNN<=(-MINIMAL/(MAXIMUN-MINIMAL))).*C31tNN;
    end
end
for im=1:size(RUIDOS4,3)
    ruidoAVG=ruidoAVG+RUIDOS(:,:,im);
    ruidoAVG2=ruidoAVG2+RUIDOS2(:,:,im);
    ruidoAVG3=ruidoAVG3+RUIDOS3(:,:,im);
    ruidoAVG4=ruidoAVG4+RUIDOS4(:,:,im);
end
ruidoAVG=ruidoAVG./(Mascara+((Mascara<=0).*1.0));
ruidoAVG2=ruidoAVG2./(Mascara2+((Mascara2<=0).*1.0));
ruidoAVG3=ruidoAVG3./(Mascara3+((Mascara3<=0).*1.0));
ruidoAVG4=ruidoAVG4./(Mascara4+((Mascara4<=0).*1.0));
figure(33)
imshow([ruidoAVG ruidoAVG2 ruidoAVG3 ruidoAVG4],[],'InitialMagnification',1024)
data=dataold;
load (data)
%% 
n=79;
N=30;
normalizar=@(imgX)(imgX-min(imgX(:)))/(max(imgX(:))-min(imgX(:)));
Circuloinit = zeros(n, n); 
[xp, yp] = meshgrid(1:n, 1:n); 
Circuloinit((xp - n/2).^2 + (yp - n/2).^2 <= (38.5).^2) = 1;
gaussCirc = fspecial('gaussian', 79, 2); 
gaussCirc=normalizar(gaussCirc);
ImgGenerada=zeros(n,n);
ImgGenerada(25:26,20:21)=ones(2,2).*0.5;
ImgGenerada(28:28,22:22)=1;
ImgGenerada(60:61,17:18)=ones(2,2).*0.9;
ImgGenerada(69:69,44:44)=1;
ImgGenerada(71:71,42:42)=1;
ImgGenerada(49:49,33:33)=0.8;
ImgGenerada(45:45,10:10)=1.3;
ImgGenerada(16:17,23:24)=ones(2,2).*0.75;
ImgGenerada(36:36,38:38)=1.5;
ImgGenerada(34:35,40:41)=ones(2,2).*0.75;
ImgGenerada(40,40)=1.6;
ImgGenerada(66:70,53:57)=[0,0,1,0,0;0,0,1,0,0;1,1,1,1,1;0,0,1,0,0;0,0,1,0,0].*0.83;
ImgGenerada(31,53:56)=[1;1;1;1].*0.63;
ImgGenerada(47:51,56)=[1 1 1 1 1].*0.93;
ImgGenerada((xp - 55).^2 + (yp - 18).^2 <= (4).^2)=0.27;
ImgGenerada=normalizar(ImgGenerada);

imgGAUS=normalizar(abs(ifftnc(fftnc(gaussCirc).*fftnc(ImgGenerada))));

ImgGeneradablu=abs(ifftnc(fftnc(PSFmod).*fftnc(ImgGenerada)));
ImgGeneradablu=normalizar(ImgGeneradablu).*Circuloinit;
ruidoFONDO=ruidoAVG3+((ruidoAVG3<=0).*ruidoAVG3(1,1));
ImgGEN=(ImgGeneradablu+((ImgGeneradablu<=0).*ruidoAVG3(1,1)/6)+ruidoFONDO);
ImgGEN=normalizar(ImgGEN);
ImgGENR=normalizar(imnoise(ImgGEN,'gaussian',0.0005,0.00002).*Circuloinit);
ImgGENR=ImgGENR+(ImgGENR==0).*ImgGEN(1,1);
ImgGENFilt=imnlmfilt(ImgGENR,'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
ImgGAUSFilt=imnlmfilt(normalizar(imnoise(imgGAUS,'gaussian',0.0005,0.00002).*Circuloinit),'ComparisonWindowSize',3,'SearchWindowSize',21,"DegreeOfSmoothing",0.02);
figure
imshow([ImgGenerada imgGAUS ImgGAUSFilt ImgGeneradablu; ImgGEN ImgGENR ImgGENFilt normalizar(C20t{7})],[],'InitialMagnification',1024)

%% filtro inverso
clc;
yy = abs(ifftnc(fftnc(ImgGENFilt)./fftnc(PSFmod+0.01)));
yy=normalizar(yy);
figure(9)
imshow([ImgGENFilt, yy],[],'InitialMagnification',1024)

NRMSEopt=456789;
for kkk=0.0001:0.005:0.4
    yyREg=normalizar(deconvreg(ImgGENFilt,PSFmod,kkk).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyREg,ImgGenerada,Circuloinit,0);
    figure(10)
    if NRMSEopt>NRMSE
        kkkopt=kkk;
%         yyREgopt=yyREg;
    end
    imshow([ImgGenerada ImgGENFilt;yy yyREg],[],'InitialMagnification',1024)
    title("K="+kkk+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    pause(0.05)
end
%%
figure(2001)
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]==ErroresIMG(yyREg,ImgGenerada,Circuloinit,0);
valores=[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM];
imshow(yyREgopt,[],'InitialMagnification',1024)
title("K="+kkkopt+" |NRMSE="+valores(5)+" |SNR="+valores(6)+" |PSNR="+valores(7)+" |SSIM="+valores(8));
%% Wiener

% yyWie2=normalizar(deconvwnr(ImgGENR,PSFmod));
% figure
% imshow([ImgGENFilt, yyWie2],[],'InitialMagnification',1024)

signal_var = var(ImgGENFilt(:));
NSR = 0.0007 / signal_var;
NRMSEopt=456789;
for kkk=0.001:0.005:0.6
    NSR =kkk;
    yyWie=deconvwnr(ImgGENR,PSFmod,NSR,0.712);
    yyWie=normalizar(yyWie);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyWie,ImgGenerada,Circuloinit,0);
    if NRMSEopt>NRMSE
        NSRopt=kkk;
        yyWiegopt=yyWie;
    end
    figure(20)
    imshow([ImgGenerada ImgGENFilt;yyWie2 yyWie],[],'InitialMagnification',1024)
    title("NSR_{input}="+kkk+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    pause(0.05)
end
%%
figure(2002)
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]==ErroresIMG(yyWiegopt,ImgGenerada,Circuloinit,0);
valoreswiener=[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM];
imshow(yyWiegopt,[],'InitialMagnification',1024)
title("NSR_{input}="+kkk+" |NRMSE="+valoreswiener(5)+" |SNR="+valoreswiener(6)+" |PSNR="+valoreswiener(7)+" |SSIM="+valoreswiener(8));
%% Richardson Lucy
pause(0.5)
NRMSEopt=456789;
for kkk=1:1:35
    yyRL=deconvlucy(ImgGENFilt,PSFmod,kkk,0.0005,Circuloinit,0.002);
    yyRLmod=normalizar(yyRL.*Circuloinit);
    yyRL=deconvlucy(ImgGENFilt,PSFFinal,kkk,0.0005,Circuloinit,0.00002);
    yyRLfinal=normalizar(yyRL.*Circuloinit);
    yyRL=deconvlucy(ImgGENFilt,gaussCirc,kkk,0.0005,Circuloinit,0.00002);
    yyRLgauss=normalizar(yyRL.*Circuloinit);
    
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yyRLmod,ImgGenerada,Circuloinit,0);
     if NRMSEopt>NRMSE
        NSRopt=kkk;
        yyWiegopt=yyRLmod;
    end
    figure(21)
    imshow([ImgGENFilt, yyRLmod;yyRLfinal,yyRLgauss],[],'InitialMagnification',1024)
    title("#_{iter}="+kkk+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    pause(0.05)
end
%%
figure(2003)
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]==ErroresIMG(yyWiegopt,ImgGenerada,Circuloinit,0);
valoreswiener=[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM];
imshow(yyWiegopt,[],'InitialMagnification',1024)
title("#_{iter}="+kkk+" |NRMSE="+valoreswiener(5)+" |SNR="+valoreswiener(6)+" |PSNR="+valoreswiener(7)+" |SSIM="+valoreswiener(8));
%% SVD
g = vec(fftnc(ImgGENFilt));
time0 = clock;
% gruido = vec(fftnc(Iruidoc));
A = diag(vec(fftnc(PSFmod)));
[U,s,V] = csvd(A);
%% Picard Graph
g = vec(fftnc(ImgGENFilt));
sigma = s;
figure(12345)
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
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(35)
    imshow([ImgGenerada ImgGENFilt; x_ka x_ka2],[],"InitialMagnification",400)
    title("\Delta_{DPF}="+lambda+" \Delta_{DPR}="+lambda2+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
    pause(0.005)
end

%% K TSVD
g = vec(fftnc(ImgGENFilt));
g2 = vec(fftnc(ImgGENR));
NRMSEopt=12345678;
for K=[1200:100:4000]
    time0 = clock;
    [x_k,rho,eta] = tsvd(U,s,V,g,K);
    [x_k2,rho,eta] = tsvd(U,s,V,g2,K);
    x_ka=ifftnc(reshape(x_k,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    x_ka2=ifftnc(reshape(x_k2,79,79));
    x_ka2=normalizar(abs(x_ka2).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    if NRMSEopt>NRMSE
        NSRopt=kkk;
        yyWiegopt=yyRLmod;
    end
    figure(36)
    imshow([ImgGenerada ImgGENFilt; x_ka x_ka2],[],"InitialMagnification",400)
     title("K_{TSVD}="+K+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
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
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    [SNR_ind2,BRISQUE2,NIQE2,PIQE2,NRMSE2,SNR2,PSNR2,SSIM2]=ErroresIMG(x_ka2,ImgGenerada,Circuloinit,0);
    figure(37)
    
    imshow([ImgGenerada ImgGENFilt; x_ka x_ka2],[],"InitialMagnification",400)
    title({["TSVD K_{LCurveF}="+reg_corner+"=>|NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM];...
        ["TSVD K_{LCurveR}="+reg_corner2+"=>|NRMSE="+NRMSE2+" |SNR="+SNR2+" |PSNR="+PSNR2+" |SSIM="+SSIM2]});
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%%
[x_kLcurv,rho,eta] = tsvd(U,s,V,g,2430);
    [x_kLcurv2,rho,eta] = tsvd(U,s,V,g2,reg_corner2);
    x_ka=ifftnc(reshape(x_kLcurv,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
figure(2004)
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]==ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
valoreswiener=[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM];
imshow(x_ka,[],'InitialMagnification',1024)
title("TSVD K_{LCurveF}="+reg_corner+" |NRMSE="+valoreswiener(5)+" |SNR="+valoreswiener(6)+" |PSNR="+valoreswiener(7)+" |SSIM="+valoreswiener(8));

%% GCV TSVD Generalized cross-validation 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_min,G,reg_param] = gcv(U,s,g,'tsvd');
[x_kGCV,rho,eta] = tsvd(U,s,V,g,reg_min);
    x_ka=ifftnc(reshape(x_kGCV,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(38)
    imshow([ImgGenerada ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("K_{GCV}="+reg_min+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% NCP TSVD Normalized cumulative periodogram
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_minNCP,G,reg_param] = ncp2(U,s,g,'tsvd');
[x_kNCP,rho,eta] = tsvd(U,s,V,g,reg_minNCP);
    x_ka=ifftnc(reshape(x_kNCP,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(329)
    imshow([ImgGenerada ImgGENFilt x_ka],[],"InitialMagnification",400)
     title("K_{NCP}="+reg_minNCP+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% QOC Quasi-optimality criterion 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_minQOC,Q,reg_param] = quasiopt(U,s,g,'tsvd');
[x_kQOC,rho,eta] = tsvd(U,s,V,g,reg_minQOC);
    x_ka=ifftnc(reshape(x_kQOC,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(352)
    imshow([ImgGenerada ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("K_{QOC}="+reg_minQOC+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
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
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_tik,ImgGenerada,Circuloinit,0);
figure(402)
imshow([ImgGenerada ImgGENFilt; x_tik x_tik2],[],"InitialMagnification",400)
title("\lambda_{TIKHONOV}="+lambda+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
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
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(422)
    imshow([ImgGenerada ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{Lcurve}="+lambdaTK+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% GCV Tikhonov Generalized cross-validation 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintk,G,reg_param] = gcv(U,s,g,'tikh');
[x_tkGCV,rho,eta] = tikhonov(U,s,V,g,reg_mintk);
    x_ka=ifftnc(reshape(x_tkGCV,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(400)
    imshow([x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{GCV}="+reg_mintk+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% NCP Tikhonov Normalized cumulative periodogram
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintkNCP,G,reg_param] = ncp2(U,s,g,'tikh');
[x_tkNCP,rho,eta] = tikhonov(U,s,V,g,reg_mintkNCP);
    x_ka=ifftnc(reshape(x_tkNCP,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(403)
    imshow([ImgGenerada ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{NCP}="+reg_mintkNCP+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
%% QOC Quasi-optimality criterion 
g = vec(fftnc(ImgGENFilt));
time0 = clock;
figure
[reg_mintkQOC,Q,reg_param] = quasiopt(U,s,g,'tikh');
[x_tkQOC,rho,eta] = tikhonov(U,s,V,g,reg_mintkQOC);
    x_ka=ifftnc(reshape(x_tkQOC,79,79));
    x_ka=normalizar(abs(x_ka).*Circuloinit);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_ka,ImgGenerada,Circuloinit,0);
    figure(404)
    imshow([ImgGENFilt x_ka],[],"InitialMagnification",400)
    title("TK\lambda_{QOC}="+reg_mintkQOC+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
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
x1=ImgGenerada.*Circuloinit;
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
x00=normalizar(x0);
%%
for taus=[0.000001:0.000025:0.0005]
    tau = taus;
    lam1=1e-10;    
    tolA = 1e-8;
    [x_twistTV]=TwIST(y,A,tau,'AT', AT,'lambda',lam1,'True_x', x1,'Psi',...
             Psi,'Phi',Phi,'Monotone',1,'Initialization',x0,'StopCriterion',1,'ToleranceA',tolA,'Verbose', 0);
    x_twistTV=normalizar(x_twistTV);
    [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(x_twistTV,x1,Circuloinit,0);
    if NRMSEopt>NRMSE
        NSRopt=kkk;
        yyWiegopt=x_twistTV;
    end
    figure(505);
    imshow([x1 y; x00  x_twistTV],[],'InitialMagnification',1024);
    title("\tau="+tau+" |NRMSE="+NRMSE+" |SNR="+SNR+" |PSNR="+PSNR+" |SSIM="+SSIM);
    pause(0.005)
end
%%
figure(2009)
[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]==ErroresIMG(yyWiegopt,ImgGenerada,Circuloinit,0);
valoreswiener=[SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM];
imshow([x_twistTVak],[],'InitialMagnification',1024)
title("\tau="+tau+" |NRMSE="+valoreswiener(5)+" |SNR="+valoreswiener(6)+" |PSNR="+valoreswiener(7)+" |SSIM="+valoreswiener(8));