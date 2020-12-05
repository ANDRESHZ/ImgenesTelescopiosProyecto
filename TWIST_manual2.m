function x_twist = TWIST_manual2(x,y,alp,bt,iterTWIST,lamb,iterTV,PSF)
%x=imagen original
%y=imagen ruidosa
%alp=factor Alpha de Twist (respecto de datos anteriores xt-1)
%bt=factor beta de Twist (respecto de datos anteriores Xt)
%iterTWIST=iteraciones de TWIST
%lamb=factor asosiado a la regularizacion  lambda*TV(f)
%iterTV=iteraciones en TVdenoising
%% OPERADORES funciones y operadores
fftnc = @(x) fftshift(fftn(fftshift(x)));
ifftnc = @(x) ifftshift(ifftn(ifftshift(x)));
vec=@(data) data(:);

K = @(x1) fftnc(PSF).*fftnc(x1);%transformada
KT = @(x1) ifftnc(x1./fftnc(PSF));%tranformada inversa
proxf = @(x1,l,TViter) (max(abs(x1)-l,0).*sign(x1));%tvdenoise(x1,l,TViter);%regularizacion (tvdenoising obsequida por libreria TWIST original)
proxg = @(x1,y1) x1 - KT(K(x1)-y1); %consitencia de datos
Gamma = @(x1,y1,l,TViter) proxf(proxg(x1,y1),l,TViter);%Factor GAMMA asociado a TWIST
%% DATOS
% yin=K(y);
yin=fftnc(y);
%% ITERACION
% alp = .6;
% bt = .5;
xt = zeros(size(x)); %xt 
xt1 = zeros(size(x)); %xt-1 incializacion
%iteraTV=3; %iteraciones
% lamb = (4+(7*2))*1.3^(-7)
xt = Gamma(xt,yin,lamb,iterTV);%incializamos xt
for k = 1:iterTWIST
    xt1 = xt;
    xt = (1-alp)*xt1 + (alp-bt)*xt + bt*Gamma(xt,yin,lamb,iterTV);
end
% imshow([abs(x) abs(KT(yin)) abs(xt)]);
% title("1)Original              2)Ruidosa            3)recuperada TwIST")
x_twist=abs(xt);