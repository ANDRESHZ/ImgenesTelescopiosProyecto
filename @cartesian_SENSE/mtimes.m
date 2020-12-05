function [res] = mtimes(a,b)
% G. Cruz, C. Prieto, Medical Image Reconstruction
% a = encoding operator, 
% b = full_k_dUa : (Ny,Nx,Nc) -> (Ny,Nx) when applying EH
% b = image_space : (Ny,Nx) -> (Ny,Nx,Nc) when applying E

if a.adjoint % EH operation
    res_coils = zeros(size(b,1),size(b,2),size(a.C,3)); % init
    for coil = 1:size(a.C,3) % number of coils
        %% apply Fourier transform
        aux_b = itok(b(:,:,coil));
        % apply sampling operator
        aux_b = aux_b.* a.U;
        %% TO DO: apply inverse Fourier transform
        aux_b = ktoi(aux_b);
        %% TO DO apply coil sensitivities
        aux_b = conj(a.C(:,:,coil)).*aux_b;
        % Save output...
        res_coils(:,:,coil) = aux_b; 
    end

    % Sum over channels
    res = sum(res_coils,3);    
    % Applying intensity correction
    res = res ./ a.coil_rss;
    res(isnan(res)) = 0;
    res(isinf(res)) = 0;

else % E operation
    res = zeros(size(b,1),size(b,2),size(a.C,3)); % init
    % Applying intensity correction
%     b = b ./ a.coil_rss;
%     b(isnan(b)) = 0; 
%     b(isinf(b)) = 0;
    for coil = 1:size(a.C,3) % number of coils
        %% TO DO: apply coil maps 
        aux_b = b;
        %% TO DO: apply Fourier transform
        aux_b = itok(aux_b);
        %% TO DO: apply sampling operator
        aux_b = (a.U).*aux_b;
        %% To DO: Apply Inverse Fourier transform
        aux_b=ktoi(aux_b);
        % Save output...
        res(:,:,coil) = aux_b;
    end

end