function res = ctranspose(a)
% G. Cruz, C. Prieto, Medical Image Reconstruction

a.adjoint = xor(a.adjoint,1);
res = a;

