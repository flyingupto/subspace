function deblurimg=SB(img,standard_deviation,blured_img, mask, lambda, mu, types_frame, level, no_iteration)

%solv the image restoration f=A*u+n, where * is convolution and A is the
%known bluring mask.
% min_u || lambda * d||+1/2||Au-f|| +mu/2||Au-d||        (1)
%blured_img---is f in f=A*u+n;
%blured_mask---is A of f=A*u+n;
%lambda-- in (1) is a vector with length silmar to the level, for example
%lambda=[2,1,1;1,1,1;1,1,1] if level=3;
%mu--in (1)
%types_frame--- reconstruciton filters of tight frame's types with
%level--DC or RE level by tight frame\
%no_iteration--number of iteration for program
% set initial guesses d_0 and b_0. Choose an appropriate set of parameters
% (lambda, mu)
%%%%%%%%%%%%20-07-2013 by Shen ZhengWei%%%%%%%%%%%%%%%%%%%%%%%
var_signal=var(img(:));
NSR=standard_deviation^2/var_signal;



imgf=fft2(img);
%use the frequency information of each point to estimate the SD of the image

 local_simga_x=sqrt(abs(imgf));
 

local_NSR=standard_deviation^2./local_simga_x;



if types_frame==1
    %dc filter
    h{1}=(1/4)*[1,2,1];
    h{2}=(sqrt(2)/4)*[1,0,-1];
    h{3}=(1/4)*[-1,2,-1];
elseif types_frame==2
    h{1}=(1/16)*[1,4,6,4,1];
    h{2}=(1/8)*[-1,-2,0,2,1];
    h{3}=(sqrt(6)/16)*[1,0,-2,0,1];
    h{4}=(1/8)*[-1,2,0,-2,1];
    h{5}=(1/16)*[1,-4,6,-4,1];
end
%the length of different tight frame, for exampe, if types_frame==1, then
%Length_TF
Length_TF=length(h);
%size of input blured_image, also is the size of output image deblurimg
[m,n]=size(blured_img);

tic;

%%%%Another compute the eignvalue method 18, April, 2016
% conj_mask=mask';
% ATf=imfilter(blured_img,conj_mask,'circular');
% 
% eigen_mask=Eigen_Mask(mask,conj_mask,mu,m,n) ;



%%%%%  compute eigen_value diagonal matrix of mask %%%%%
 Sbig = psf2otf(mask, [m, n]); % psf2otf is first convert mask to Forier domain.
 ATA = conj(Sbig).*Sbig;
ATy = conj(Sbig).*fft2(blured_img);
% I_matrix=eye(m);
wLevel=-1/2;
%Compute the weighted thresholding parameters.这一步是计算针对不同分解的小波的频道，给予不同的
%阈值化处理的权重，当然，后边在函数wthresh中使用时，又除以了lambda.
lambdaLevel=getwThresh(lambda,wLevel,Length_TF,h);

%%%%%%%%%%%%%%
u=zeros(m,n);
b=QASDCML(u,level,types_frame);
d=b;
%cc=QASDCML(blured_img,level,types_frame);

%start iteration
for no_i=1:no_iteration   
    % solve u.
        for l=1:level
            for i=1:Length_TF
                for j=1:Length_TF
                   % b{l}{i,j}=b{l}{i,j}./mu;%标准化处理，但是效果不好
                    C{l}{i,j}=d{l}{i,j}-b{l}{i,j};
                end
            end
        end
        TempC=QASREML(C,level,types_frame);
       
   
         
         %u =ifft2(fft2(mu*TempC+ATf)./(eigen_mask+local_NSR));%都是在fourier域内进行操作
          %u=ifft2((fft2(mu*TempC)+ATy)./(ATA+mu));
           u =ifft2((fft2(mu*TempC)+ATy)./(ATA+mu+local_NSR));
         u=real(u);
       %  u1=real(u1)
         
    %second slove d. get the decompostion of u from the above step     
      C=QASDCML(u,level,types_frame);
      %C1=QASDCML(u1,level,types_frame);
      
      for l=1:level
            for i=1:Length_TF
                for j=1:Length_TF
                    temp{l}{i,j}=C{l}{i,j}+b{l}{i,j};%% how to define the norm of d{l}{i,j}?
                    d{l}{i,j}=wthresh(temp{l}{i,j},'s', lambdaLevel{l}{i,j}./mu); 
                    %by the Matlab inbulit function wthresh to threshold .这一步去软阈值去噪处理。 
                end
            end
      end
      
      
      
      
      
      
% d=GetwThresh(coeffs, mu,lambda,level,types_frame, types_norm)
    %compute the threshold 
    
    %d=GetwThresh(temp,mu,lambda,level,types_frame, types_norm);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % slove b %%%%%%%%%%%%%%%   
        for l=1:level
            for i=1:Length_TF
                for j=1:Length_TF
                    b{l}{i,j}=b{l}{i,j}+(C{l}{i,j}-d{l}{i,j});
                end
            end 
        end
end
%d{1}{1,1}=cc{1}{1,1};
%u=QASREML(d,level,types_frame);

u(u>255)=255;u(u<0)=0;
deblurimg=u;

t_1 = toc;
disp(['Execution time of SplitBregman: ' num2str(t_1) ' sec']);

function lambdaLevel=getwThresh(lambda,wLevel,Level,D)
%Level is the decomposition level, Level=1
nfilter=1;
nD=length(D);
if wLevel<=0
    for ki=1:Level
        for ji=1:nD
            for jj=1:nD
                lambdaLevel{ki}{ji,jj}=lambda*nfilter*norm(D{ji})*norm(D{jj});
            end
        end
        nfilter=nfilter*norm(D{1});
    end
else
    for ki=1:Level
        for ji=1:nD
            for jj=1:nD
                if ji==1 && jj==1
                    lambdaLevel{ki}{ji,jj}=0;
                else
                    lambdaLevel{ki}{ji,jj}=lambda*nfilter;
                end
            end
        end
        nfilter=nfilter*wLevel;
    end 
    
%     %%%%%%%
%     for ki=1:Level
%         for ji=1:nD
%             for jj=1:nD
%                 muLevel{ki}{ji,jj}=mu*nfilter;
%             end
%         end
%         nfilter=nfilter*wLevel;
%     end
end


