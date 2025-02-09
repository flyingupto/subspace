
function [X_MC,Z_MC,U_MC] =SGSA_FFT(FB,F2B,rho,alpha,y,FBC,N,N_MC,types_frame,level,k)
%-------------------------------------------------------------------------%
% This function computes the SPA algorithm to solve the linear inverse
% problem y = A*x + n associated to the image deconvolution problem.

% INPUTS:
% FB: counterpart of the blur operator in the Fourier domain.
% F2B: same as FB with coefficients equal to |FB|.^2.
% rho: user-defined standard deviation of the variable of
%      interest x.
% alpha: user-defined hyperparameter of the prior p(u).
% y: observations (2D-array).
% FBC: conjugate of FB.
% N: dimension of x (2-D array).
% N_MC : total number of MCMC iterations.
%types_frame--- reconstruciton filters of tight frame's types with
%level--DC or RE level by tight frame\

% OUTPUT:
% X_MC,Z_MC,U_MC: samples (2D-array) from the joint
% posterior.
%%%%%%%%%%%%22-04-2023 by Wang Shangfei%%%%%%%%%%%%%%%%%%%%%%%
tic;
disp(' ');
disp('BEGINNING OF THE SAMPLING');

%-------------------------------------------------------------------------
% Initialization
N2=N^2;
n=N2;
p=N2;
if types_frame==1
    Length_TF=3;
elseif types_frame==2
    Length_TF=5;
end

beta=cell(Length_TF,Length_TF);
beta(:)={rand(p,1)};
lambda2 =cell(Length_TF,Length_TF);
lambda2(:)={rand(p,1)};
tau2 =ones(Length_TF,Length_TF);
sigma2=1;
nu=cell(Length_TF,Length_TF);
nu(:)={rand(p,1)};
xi=ones(Length_TF,Length_TF);
% define matrices to store the iterates
X_MC = zeros(N,N,N_MC);
Z_MC = zeros(N,N,N_MC);
U_MC = zeros(N,N,N_MC);
% initialize the latter matrices
X_MC(:,:,1) = rand(N,N)*255;
Z_MC(:,:,1) = rand(N,N)*255;
U_MC(:,:,1) = rand(N,N)*255;
%-------------------------------------------------------------------------
% Gibbs sampling
for t = 1:(N_MC-1)
    % 2.2.1. Sampling the variable of interest x
    zeta= sqrt(0.5) * (randn(N,N) + sqrt(-1)*randn(N,N));
    delta=sqrt(0.5) * (randn(N,N) + sqrt(-1)*randn(N,N));
    Lambda_F= F2B./sigma2  + (1./ rho^2);
    b= FBC .* fft2(y)./sigma2+ fft2(Z_MC(:,:,t) - U_MC(:,:,t))./rho^2+FBC.*zeta./sqrt(sigma2)+delta./rho;
    
    x_F= b./ Lambda_F;
    X_MC(:,:,t+1) = real(ifft2(x_F));
    e = y - real(ifft2(FB.*x_F));
    scale = sum(sum(e.*e))/2+1e-4;
    sigma2 = 1./ gamrnd(n/ 2+1, 1./scale);
    
    % 2.2.2. Sampling z=WTbeta
    f=QASDCML(X_MC(:,:,t+1)+U_MC(:,:,t),level,types_frame);
    for l=1:level
        for i=1:Length_TF
            for j=1:Length_TF
                if(i==1&&j==1)
                else
                    y_tilde=reshape(f{l}{i,j},[p,1]);
                    if t>0
                        omega=sort(lambda2{i,j},'descend');
                        m=find(lambda2{i,j}>=omega(k));
                        y_tilde=y_tilde(m);
                        D=tau2(i,j) * lambda2{i,j}(m);
                        zeta=randn(k,1)./sqrt(D);
                        delta=randn(k,1)./rho;
                        beta{i,j}(m)=(y_tilde./rho^2 +zeta+delta)./(1./D  + (1 / rho^2));
                        shape=(k+1)/2;
                    else
                        D=tau2(i,j) * lambda2{i,j};
                        zeta=randn(k,1)./sqrt(D);
                        delta=randn(k,1)./rho;
                        beta{i,j}=(y_tilde./rho^2 +zeta+delta)./(1./D  + (1 / rho^2));
                        shape=(p+1)/2;
                    end
                    f{l}{i,j}= reshape(beta{i,j},[N,N]);
                    %% Sample lambda2
                    scale = 1./nu{i,j}(m) + beta{i,j}(m).^2./2./tau2(i,j);
                    lambda2{i,j}(m)= 1 ./ exprnd(1./scale);
                    
                    %% Sample tau2
                    
                    scale = 1/xi(i,j) + sum(beta{i,j}(m).^2./lambda2{i,j}(m))/2;
                    tau2(i,j)= 1 / gamrnd(shape, 1/scale);
                    
                    %% Sample nu
                    scale = 1 + 1./lambda2{i,j}(m);
                    nu{i,j}(m) = 1 ./ exprnd(1./scale);
                    
                    %% Sample xi
                    scale = 1 + 1/tau2(i,j);
                    xi(i,j) = 1 / exprnd(1/scale);
                    
                end
            end
        end
    end
    Z_MC(:,:,t+1) =QASREML(f,level,types_frame);
    clear eps x0 u0 z0;
    
    % 2.2.3. Sampling u
    cov = (alpha^2 * rho^2) / (alpha^2 + rho^2);
    moy = (Z_MC(:,:,t+1)-X_MC(:,:,t+1)) * alpha^2 / (rho^2 + alpha^2);
    moy = reshape(moy,[N^2,1]);
    u0 = moy + randn(N^2,1).* sqrt(cov);
    U_MC(:,:,t+1)= reshape(u0,[N,N]);
    clear moy cov eps u0;
    %    U_MC(:,:,t+1)= zeros(N,N);
    
    
    waitbar(t/N_MC);
    
end
%-------------------------------------------------------------------------
t_1 = toc;
%close(h);
disp('END OF THE GIBBS SAMPLING');
disp(['Execution time of the Gibbs sampling: ' num2str(t_1) ' sec']);

end