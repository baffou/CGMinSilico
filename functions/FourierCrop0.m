function [SImout, SRfout] = FourierCrop0(FIm,FRf,cropParams)
% Compute the cropped images, according to the FcropParams object cropParams
% SImout: Centered crop in the Fourier space
% SRfout: Centered crop in the Fourier space (for the reference)
% cropParams: the parameters of the crop (FcropParams object)

if size(FIm)~=size(FRf)
    size(FIm)
    size(FRf)
    error('Interferogram and reference images do not have the same size')
end
[Ny,Nx] = size(FIm);

R = cropParams.R;

if length(R)==1
    rx = R;
    ry = R;
else
    rx = R(1);
    ry = R(2);
end

[xx,yy] = meshgrid(1:Nx, 1:Ny);

R2C = (xx  -Nx/2-1-cropParams.shiftx).^2/rx^2 + (yy - Ny/2-1-cropParams.shifty).^2/ry^2;
circle = (R2C < 1); % mask

% crop and demodulation in the Fourier space
Imout = FIm.*circle;
Rfout = FRf.*circle;
SImout= circshift(Imout,[-cropParams.shifty -cropParams.shiftx]);
SRfout= circshift(Rfout,[-cropParams.shifty -cropParams.shiftx]);


