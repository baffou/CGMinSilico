function [W,Int] = CGMprocess(Itf,ItfRef,alpha,crop)
% Function that processes the intensity and OPD images from interferograms
% G. Baffou, CNRS, Jan 2022

FIm = fftshift(fft2(Itf));
FRf = fftshift(fft2(ItfRef));

%    fprintf('ZERO ORDER\n')
[FImc,FRfc] = FourierCrop0(FIm,FRf,crop{1});
Im_int = ifft2(ifftshift(FImc));
Rf_int = ifft2(ifftshift(FRfc));

%    fprintf('FIRST ORDER ALONG XD\n')
[FImc,FRfc] = FourierCrop0(FIm,FRf,crop{2});
%figure('Name','QLSIprocess.m/ TF Fourier crop 1st order'),imagetf(FImc)

Im_DW1 = ifft2(ifftshift(FImc));
Rf_DW1 = ifft2(ifftshift(FRfc));

% fprintf('FIRST ORDER ALONG YD\n')
[FImc,FRfc] = FourierCrop0(FIm,FRf,crop{3});

Im_DW2 = ifft2(ifftshift(FImc));
Rf_DW2 = ifft2(ifftshift(FRfc));

%% computation of the I & OPD images

% intensity image
Int=abs(Im_int)./abs(Rf_int);

% phase gradient images

DW1 = angle(Im_DW1.*conj(Rf_DW1))* alpha;
DW2 = angle(Im_DW2.*conj(Rf_DW2))* alpha;

DWx = crop{2}.angle.cos*DW1-crop{2}.angle.sin*DW2;
DWy = crop{2}.angle.sin*DW1+crop{2}.angle.cos*DW2;

DWx = DWx-mean( DWx(:));
DWy = DWy-mean( DWy(:));

W = intgrad(DWx,DWy);
