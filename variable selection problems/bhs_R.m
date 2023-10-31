
function [beta,iter] = bhs_R(Xorg, yorg,biao)
[n, p] = size(Xorg);
%% Normalise data
[X, normX, y]=standardise(Xorg, yorg);
%% Initial values
sigma2 = 100;
lambda2 = rand(p, 1);
tau2 = 1;
xi = 1;
nu = ones(p,1);
shape_sigma2= (n + p) / 2;
%shape_sigma2= (n + 1) / 2;
shape_tau2 = (p + 1)/2;
sigma = sqrt(sigma2);
%% Gibbs sampler
iter = 0;
%Phi n*p
%sigma2 1*1
%lambda2 p*1
%D p*1
%u p*1
%b p*1
%X n*p
while(sigma>biao)
    iter = iter + 1;
    %% Sample from the conditional posterior dist. for beta
    sigma = sqrt(sigma2);
    Phi=X ./ sigma;
    D=sigma2 * tau2 * lambda2;
    u = randn(p,1) .* sqrt(D);
    v = Phi*u + randn(n,1);
    
    Dpt = bsxfun(@times, Phi', D);
    w = (Phi*Dpt + eye(n)) \ (y ./ sigma - v);
    b = u + Dpt*w;
    
    
    
    scale = 1 + 1./lambda2;
    nu = 1 ./ exprnd(1./scale);
    
    
    e = y - X*b;
    
    b2=b.^2./2;
    ssu= sum(b2./lambda2);
    
    scale = 1/xi +ssu/sigma2;
    tau2 = 1 / gamrnd(shape_tau2, 1/scale);
    
    scale = e'*e/2 + ssu/tau2;
    sigma2 = 1 / gamrnd(shape_sigma2, 1/scale);
    
    %% Sample lambda2
    scale = 1./nu + b2./tau2./sigma2;
    lambda2 = 1 ./ exprnd(1./scale);
    
    scale = 1 + 1/tau2;
    xi = 1 / exprnd(1/scale);
end
%% Re-scale coefficients
beta = bsxfun(@rdivide, b, normX');
end


%% Standardise the covariates to have zero mean and x_i'x_i = 1
function [X,stdX,y] = standardise(X, y)
%% params
n = size(X, 1);
meanX = mean(X);
stdX = std(X,1) * sqrt(n);
%% Standardise Xs
X = bsxfun(@minus,X,meanX);
X = bsxfun(@rdivide,X,stdX);
%% Standardise ys (if neccessary)
if(nargin == 2)
    meany = mean(y);
    y = y - meany;
end
%% done
end
