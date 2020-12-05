function [Data] = ktoi(Data, Dims)
% ITOK - G. Cruz, C. Prieto, Medical Image Reconstruction
% 

if nargin<2,
  Data = fftshift(ifftn(ifftshift(Data)));
else
  for d=[1:length(Dims)],
    Data = fftshift(ifft(ifftshift(Data),[],Dims(d)));
  end
end
return
