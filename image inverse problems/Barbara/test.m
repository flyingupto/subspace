clc
clear
close all;
addpath('../src/'); 

% 1.1. Load original 512 x 512 image 
name = 'Barbara.bmp';
Img =  double(imread(name));

% 1.2. Define the regularization 
psf = [[0,-1,0];[-1,4,-1];[0,-1,0]];
[FL,FLC,F2L,~] = HXconv(Img,psf,'Hx');
delta = 1e-1;
FL = -FL + delta;
FLC = -FLC + delta;
F2L = abs(FL).^2;

% 1.3. Define the blurring kernel and its associated Fourier matrices 
mask_type = 'gaussian';
mask_size = 3;
mask_deviation =1;
standard_deviation=2;
mask=fspecial(mask_type,mask_size,mask_deviation);
[FB,FBC,F2B,Bx] = HXconv(Img,mask,'Hx');

% 1.4. Apply the blurring operator on the original image
blured_Img=imfilter(Img,mask, 'circular');
N = size(blured_Img,1);
blured_Img=blured_Img+standard_deviation* randn(N,N);

% 1.5. Define the parameters of SPA
rho =10;
alpha = 1;%改变这个没用

% 1.6. Define MCMC parameters
N_MC =110; % total number of MCMC iterations
N_bi =90; % number of burn-in iterations

% 1.7. Other parameters and precomputing
gamma = 6e-3; % regularization parameter (fixed here)
D =standard_deviation * ones(N,N);
D = D.^(-2); % precision matrix associated to the likelihood
mu1 = 0.99 / max(D(:)); % parameter used in AuxV1 method embedded in SPA

var_signal=var(Img(:));
NSR=standard_deviation^2/var_signal;
lambda=NSR;
mu=0.1;
types_frame=1;
level=1;
p=N*N;
k=round(p*1/10);
[X_SPA,Z_SPA,U_SPA] = SPA(D,mu1,FB,F2B,rho,alpha,blured_Img,FBC,gamma,F2L,N,N_MC);
[X_SGSA_FFT,Z_SGSA_FFT,U_SGSA_FFT]=SGSA_FFT(FB,F2B,rho,alpha,blured_Img,FBC,N,N_MC,types_frame,level,k);
x_SB=SB(Img,standard_deviation, blured_Img, mask, lambda, mu, types_frame,level, N_MC);

% Take out the last sample
x_SPA=X_SPA(:,:,end);
x_SGSA_FFT=X_SGSA_FFT(:,:,end);
% Plot the results
p1=200;
p2=270;
PLOT(p1,p2,N,mask_size,mask_deviation,standard_deviation,Img,blured_Img,x_SPA,x_SGSA_FFT,x_SB,name);
