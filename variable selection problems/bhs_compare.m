%%%%%%%%%%%%08-06-2023 by Wang Shangfei%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters:
%   A        = regressor matrix [n * p]
%   yorg     = response vector  [n * 1]
%   nsamples = number of samples for the Gibbs sampler (nsamples > 0)
%   burnin   = number of burnin (burnin >= 0)
%   thin     = thinning (thin >= 1)
%
% Returns:
%   x1       = regression parameters of the subspace algorithm  [p * nsamples]
%   x2       = regression parameters of the non-subspace algorithm  [p * nsamples]
%   time     = time cost of two algorithms [2 * 1]
%  
% Example:
% % Load a dataset:
% % Run horseshoe sampler. Normalising the data is not required.
% [x1,x2,time] = bhs_compare(A, y_noisy,nsamples,burnin, thin);
%
%
% References:
% A simple sampler for the horseshoe estimator
% E. Makalic and D. F. Schmidt
% arAiv:1508.03884
%
% Fast sampling with gaussian scale mixture priors in high-dimensional regression
% Anirban Bhattacharya, Antik Chakraborty, and Bani K Mallick
% arAiv:1506.04778
%
% Scalable approximate mcmc algorithms for the horseshoe prior.
% James Johndrow, Paulo Orenstein, and Anirban Bhattacharya.
% arAiv:1705.00841
%
% (c) Copyright Wang Shangfei and Shen Zhengwei
function [x1,x2,time] = bhs_compare(A, y_noisy,nsamples,burnin, thin)
[n, p] = size(A);
%% Normalise data
[A, normA, y]=standardise(A, y_noisy);
%% Initial values
sigmx2 = 1;
lambdx2 = 100*rand(p, 1);
tau2 = 1;
xi = 1;
nu = ones(p,1);
%% Gibbs sampler
k = 0;
iter = 0;
x1=zeros(p,nsamples);
x2=zeros(p,nsamples);
while(iter < burnin)
    iter = iter + 1;
    %% Sample from the conditional posterior dist. for x
    sigma = sqrt(sigmx2);
    
    D=tau2 * lambdx2;
    u = randn(p,1).*sigma .* sqrt(D);
    v = A*u + randn(n,1)*sigma;
    Dpt = bsxfun(@times, A', D);
    w = (A*Dpt + eye(n)) \ (y - v);
    b= u + Dpt*w;
    
    scale = 1 + 1./lambdx2;
    nu = 1 ./ exprnd(1./scale);
    
    e = y - A*b;
    
    b2=b.^2./2;
    ssu= sum(b2./lambdx2);
    
    scale = 1/xi +ssu/sigmx2;
    tau2 = 1 / gamrnd((p + 1)/2, 1/scale);
    
    scale = e'*e/2 + ssu/tau2;
    sigmx2 = 1 / gamrnd((n + p) / 2, 1/scale);
    %% Sample lambdx2
    scale = 1./nu + b2./tau2./sigmx2;
    lambdx2 = 1 ./ exprnd(1./scale);
    
    scale = 1 + 1/tau2;
    xi = 1 / exprnd(1/scale);
end
s2=sigmx2;
t2=tau2;
l2=lambdx2;
xx=xi;
tic
while(k<nsamples)
    iter = iter + 1;
    sigma = sqrt(sigmx2);
    omega=sort(lambdx2);
    m=find(lambdx2>=omega(p*9/10));
    r=length(m);
    D=tau2 * lambdx2(m,:);
    u = randn(r,1).*sigma .* sqrt(D);
    v = A(:,m)*u + randn(n,1)*sigma;
    Dpt = bsxfun(@times, A(:,m)', D);
    w = (A(:,m)*Dpt + eye(n)) \ (y - v);
%   if(~isempty(find(isnan(w), 1)))
%      break;
%   end
%   b=zeros(p,1);
    b(m,:)= u + Dpt*w;
    
    scale = 1 + 1./lambdx2;
    nu = 1 ./ exprnd(1./scale);
    
    e = y - A*b;
    
    b2=b.^2./2;
    ssu= sum(b2(m,:)./lambdx2(m,:));
    scale = 1/xi +ssu/sigmx2;
    tau2 = 1 / gamrnd((r + 1)/2, 1/scale);
    
    scale = e'*e/2 + ssu/tau2;
    sigmx2 = 1 / gamrnd((n + r) / 2, 1/scale);
    %% Sample lambdx2
    scale = 1./nu + b2./tau2./sigmx2;
    lambdx2 = 1 ./ exprnd(1./scale);
    
    scale = 1 + 1/tau2;
    xi = 1 / exprnd(1/scale);
        % thinning
        if(mod(iter,thin) == 0)
            k = k + 1;
            x1(:,k) = b;
        end
end
%% Re-scale coefficients
x1 = bsxfun(@rdivide, x1, normA');
toc
time(1)=toc;
k = 0;
iter= burnin;
sigmx2=s2;
tau2=t2;
lambdx2=l2;
xi=xx;
tic
while(k<nsamples)
    iter = iter + 1;
    %% Sample from the conditional posterior dist. for x
    sigma = sqrt(sigmx2);
    D=tau2 * lambdx2;
    u = randn(p,1).*sigma .* sqrt(D);
    v = A*u + randn(n,1)*sigma;
    Dpt = bsxfun(@times, A', D);
    w = (A*Dpt + eye(n)) \ (y - v);
    b= u + Dpt*w;
    
    scale = 1 + 1./lambdx2;
    nu = 1 ./ exprnd(1./scale);
    
    
    e = y - A*b;
    
    b2=b.^2./2;
    ssu= sum(b2./lambdx2);
    
    scale = 1/xi +ssu/sigmx2;
    tau2 = 1 / gamrnd((p + 1)/2, 1/scale);
    
    scale = e'*e/2 + ssu/tau2;
    sigmx2 = 1 / gamrnd((n + p) / 2, 1/scale);
    %% Sample lambdx2
    scale = 1./nu + b2./tau2./sigmx2;
    lambdx2 = 1 ./ exprnd(1./scale);
    
    scale = 1 + 1/tau2;
    xi = 1 / exprnd(1/scale);
        % thinning
        if(mod(iter,thin) == 0)
            k = k + 1;
            x2(:,k) = b;
        end
end
x2 = bsxfun(@rdivide, x2, normA');
toc
time(2)=toc;
end


%% Standardise the covariates to have zero mean and a_i'a_i = 1
function [A,stdA,y] = standardise(A, y)
%% params
n = size(A, 1);
meanA = mean(A);
stdA = std(A,1) * sqrt(n);
%% Standardise As
A = bsxfun(@minus,A,meanA);
A = bsxfun(@rdivide,A,stdA);
%% Standardise ys (if neccessary)
if(nargin == 2)
    meany = mean(y);
    y = y - meany;
end
%% done
end
