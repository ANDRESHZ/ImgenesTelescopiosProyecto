function [SNR_ind,BRISQUE,NIQE,PIQE,NRMSE,SNR,PSNR,SSIM]=ErroresIMG(yanalizar,xorg,mascara,imprimir)
%yanalizar=imagen o matriz a analizar
%xorg=imagen de referencia u original
%mascara=mascara de valores para ponderar o discrimisnar las medisiones.

%%%%%TOMADO DE https://la.mathworks.com/help/images/image-quality-metrics.html
%Un modelo CALIDAD ESPACIAL BRISQUE est� formado en una base de datos de im�genes con distorsiones conocidas, 
%y BRISQUE se limita a evaluar la calidad de las im�genes con el mismo tipo de distorsi�n.
%BRISQUE es , lo que significa que las puntuaciones de calidad subjetivas acompa�an las im�genes
%de entrenamiento.opinion-aware
%%%%%
%Evaluador de Calidad de Imagen Natural (NIQE).niqe Aunque un modelo NIQE est� entrenado en una base de datos 
%de im�genes pr�stinas, NIQE puede medir la calidad de las im�genes con distorsi�n arbitraria. NIQE es , y no 
%utiliza puntuaciones de calidad subjetivas.opinion-unaware La contrapartida es que la puntuaci�n NIQE de una 
%imagen podr�a no correlacionarse, as� como la puntuaci�n BRISQUE con la percepci�n humana de la calidad.
%%%%%
%Evaluador de calidad de imagen basado en la percepci�n (PIQE).piqe El algoritmo PIQE no reconoce la opini�n, 
%lo que significa que no requiere un modelo entrenado.Unsupervised PIQE puede medir la calidad de las im�genes
%con distorsi�n arbitraria y en la mayor�a de los casos realiza un rendimiento similar al NIQE. PIQE estima la
%distorsi�n en bloque y mide la varianza local de bloques distorsionados perceptiblemente para calcular la 
%puntuaci�n de calidad.
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
%%Calcula errores con Referencias si se envia una matriz de iguaL TAMA�O
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