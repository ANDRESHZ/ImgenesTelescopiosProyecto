function [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yanalizar,xorg,mascara,imprimir)
%yanalizar=imagen o matriz a analizar
%xorg=imagen de referencia u original
%mascara=mascara de valores para ponderar o discrimisnar las medisiones.

%%%%%TOMADO DE https://la.mathworks.com/help/images/image-quality-metrics.html
%Un modelo CALIDAD ESPACIAL BRISQUE está formado en una base de datos de imágenes con distorsiones conocidas, 
%y BRISQUE se limita a evaluar la calidad de las imágenes con el mismo tipo de distorsión.
%BRISQUE es , lo que significa que las puntuaciones de calidad subjetivas acompañan las imágenes
%de entrenamiento.opinion-aware
%%%%%
%Evaluador de Calidad de Imagen Natural (NIQE).niqe Aunque un modelo NIQE está entrenado en una base de datos 
%de imágenes prístinas, NIQE puede medir la calidad de las imágenes con distorsión arbitraria. NIQE es , y no 
%utiliza puntuaciones de calidad subjetivas.opinion-unaware La contrapartida es que la puntuación NIQE de una 
%imagen podría no correlacionarse, así como la puntuación BRISQUE con la percepción humana de la calidad.
%%%%%
%Evaluador de calidad de imagen basado en la percepción (PIQE).piqe El algoritmo PIQE no reconoce la opinión, 
%lo que significa que no requiere un modelo entrenado.Unsupervised PIQE puede medir la calidad de las imágenes
%con distorsión arbitraria y en la mayoría de los casos realiza un rendimiento similar al NIQE. PIQE estima la
%distorsión en bloque y mide la varianza local de bloques distorsionados perceptiblemente para calcular la 
%puntuación de calidad.
%%%%%
%% AJUSTA Y PONDERA SEGUN LA MASCARA
    if mean(size(yanalizar)==size(mascara))==1
        yanalizar=yanalizar.*mascara;
    else
        if imprimir>0
            disp('Se continua sin mascara');
        end
    end
     if mean(size(xorg)==size(mascara))==1
        xorg=xorg.*mascara;
    end
%% Calcula los errores sin referencia
    SNR_ind=mean2(yanalizar)/var(yanalizar(:));
    BRISQUE=brisque(yanalizar);
    NIQE=niqe(yanalizar);
    PIQE=piqe(yanalizar);
    if imprimir>0
        disp("SNR_ind="+SNR_ind+"    BRISQUE="+BRISQUE+"    NIQE="+NIQE+"    PIQE="+PIQE);
    end
%%Calcula errores con Referencias si se envia una matriz de iguaL TAMAÑO
    if mean(size(xorg)==size(yanalizar))==1
        NRMSE=norm(xorg-yanalizar,'fro')^2/norm(xorg,'fro')^2;
        [PSNR SNR]=psnr(yanalizar,xorg);
        SSIM=ssim(yanalizar,xorg);
        if imprimir>0
            disp("NRMSE="+NRMSE+"    SNR="+SNR+"    PSNR="+PSNR+"    SSIM="+SSIM);
        end
    else
        if imprimir>0
            disp("NO se calcula [NRMSE,SNR, PSNR,SSIM]");
        end
        NRMSE=-123456789;
        PSNR=-123456789;
        SNR=-123456789;
        SSIM=-123456789;
    end
end