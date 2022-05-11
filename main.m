%% InSilex ALGORITHM
%% Program that simulates experimental cross-grating phase microscopy images
% G. Baffou
% CNRS
% Jan 2022

% Associated with the article:
%  Cross-grating phase microscopy (CGM):
%  In-silico experiment (insilex) algorithm, noise and accuracy
%  Optics Communication 2022
%  B. Marthy and G. Baffou
%  Github repository: https://github.com/baffou/CGMinSilico

% This program can either image noise (GaussianOPD=0), or image a Gaussian
% phase profile (GaussianOPD=1).

clear
close all

addpath(genpath(pwd))

%% Parameters

lambda = 625e-9;  % actual wavelength sent to the grating [m]
eD = lambda;      % etching depth of the checkerboard [m]. Must equal lambda to set a pi phase shift.
camPxSize = 6.5e-6;   % camera chip pixel size [m]
d0 = 1e-3;        % distance of the grating from the image plane [m]
zeta = 3;         % zeta factor: Gamma/(2*camPxSize)
Nim = 25;         % Number of averaged interferograms
w = 40000;        % full well capacity of the camera
Npx = 600;        % Desired final image size [px], must be square
zoom0 = 1;        % Relay lens zoom
beta = acos(3/5); % tilt of the cross-grating [rad]. Possible values to maintain periodicity: 0 or acos(3/5)

GaussianOPD = 1;  % false: flat OPD profile ; true: Gaussian OPD profile
Gampl = 600e-9;   % Gaussian profile amplitude [m]
Gradius = Npx/10; % Gaussian profile radius [m]

%% Definition of some parameters
nCell = 20;     % Initial overdimensioned size[px] of the grating unit cell: 6*nCellx6*nCell
Gamma = 2*zoom0*zeta*camPxSize;       % size [m] of the unit cell of the grating

%% construction of the interferograms (Itf & Ref) according to Fig 2.

% (Fig 2a) Build the unit cell :
    grexel = QLSIunitCell(nCell,pi*lambda/eD,Gamma); 
% (Fig 2b) 5x5-tile and rotate by acos(4/5) the unit cell to form the superunit cell
% image, enabling periodicity when further tiling :
    superUnit = grexel.TileRot5(beta);
% (Fig 2c) Resizing of the superunit cell by imresize to match the pixel size of the camera :
    superUnitPixelized = superUnit.redimension(camPxSize,zeta*zoom0/(3*nCell));
% (Fig2d) generation of the final E field by tiling until reaching the desired Npx
% px number of the image :
    Grating = tile(superUnitPixelized,Npx);

Grating.figure
title('real part of the grating transmittance')

if GaussianOPD == 1 % Gaussian OPD profile of amplitude Gampl
    [X,Y] = meshgrid(1:Npx,1:Npx);
    X = X-mean(X(:));
    Y = Y-mean(Y(:));
    OPD = Gampl*exp(-(X.^2+Y.^2)/Gradius.^2); % (Fig 2f)
    Pha = 2*pi/lambda*OPD;
else % uniform phase profile, just to study the image noise
    Pha = zeros(Npx,Npx); % (Fig 2f)
end

Emodel = Grating;
Emodel.im = exp(-1i*Pha);
EmodelRef = Grating;
EmodelRef.im = ones(Npx);

Emodelb = Emodel.propagation(lambda,-d0);% Backpropagation of the light from the camera chip
Egrating = Grating.*Emodelb;
EmodelRefb = EmodelRef.propagation(lambda,-d0);
ErefGrating = Grating.*EmodelRefb;

% (Fig 2e,g) Propagate the light after the unit cell:
E    = Egrating.propagation(lambda,d0);
Eref = ErefGrating.propagation(lambda,d0);

% (Fig 2h,i) compute the intensity images:
E20 = E.square();
[E2,fac] = E20.setI0(w*Nim); % set the max counts of the interferogram
E2Ref0 = Eref.square();
E2Ref = E2Ref0.timesC(fac);% Apply the same correction factor to the reference interferogram

hfig=figure;
subplot(1,2,1)
E.figure(hfig)
title(sprintf('Real part of E-field at the camera plane\n Interfero'))
subplot(1,2,2)
Eref.figure(hfig)
title(sprintf('Real part of E-field at the camera plane\n Reference'))
drawnow

cut = 1;% =1 or 2. to make artificial crosses in the Fourier space and better
        % match what is observed experimentally. Otherwise, set to 0.
        % Note it reduces the image size from Npx to Npx-cut.

% Add the shoit noise on the interferograms:  
Itf    = poissrnd(   E2.im(1:end-cut,1:end-cut));
ItfRef = poissrnd(E2Ref.im(1:end-cut,1:end-cut));

% plot of the interferogram and its Fourier transform:
figure
subplot(1,2,1)
imagebw(Itf)
title('interferogram')
subplot(1,2,2)
imagetf(Itf)
title('Fourier transform of the interferogram')

%% definition of the crops of Fig2j:
crop = cell(3,1);
Rcrop = Npx/(2*zeta*zoom0);
crop{1} = FcropParams(Npx/2+1                   ,Npx/2+1                   ,Rcrop,Npx);
crop{2} = FcropParams(Npx/2+1+2*Rcrop*cos(beta),Npx/2+1+2*Rcrop*sin(beta),Rcrop,Npx);
crop{3} = FcropParams(Npx/2+1-2*Rcrop*sin(beta),Npx/2+1+2*Rcrop*cos(beta),Rcrop,Npx);

%% CGM process to get the OPD and Intensity images (Fig2m)
alpha = camPxSize*Gamma/(4*pi*d0);
[W,Int] = CGMprocess(Itf,ItfRef,alpha,crop);

%%
figure
subplot(2,2,1)
imagebw(Int)
title('Intensity image')
subplot(2,2,2)
imageph(W)
subplot(2,2,3)
plot(Int(round(end/2),:),'Color','#444444','lineWidth',1.2)
axis square
subplot(2,2,4)
hold on
plot(W(round(end/2),:)*1e9,'Color','#D95319','lineWidth',1.5)
plot(OPD(round(end/2),:)*1e9,'k-.','lineWidth',0.8)
ylabel('[nm]')
axis square
legend({'in silico','model'})
