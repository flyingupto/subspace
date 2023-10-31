clc
clear
p=10000;    % 特征维度
n =1800;    % 样本数量
mu=10*rand(1,p);
sigma=20*rand(1,p);
Xorg=zeros(n,p);
for i=1:p
    Xorg(:,i)=mu(i)+sigma(i)*randn(1,n);
end

btrue=zeros(p,1);
k=p/100;
a=randperm(p);
a=sort(a(1:k));
byou=zeros(k,1);

mean=5;
for i=1:k
    byou(i)=mean+2*k*rand;
    btrue(a(i))=byou(i);
end

% 计算对应的一组Y值
yorg = Xorg* btrue ;

noise =randn(n,1);
biao=1;
y_noisy = yorg +biao*noise;

data2 =byou;
%子空间
l=5;
data1 = btrue;
tong=zeros(8,4*l);


s=1;
for i=1:l
    tic
    [b,iter,m,sigma1]= bhs_zu(Xorg, y_noisy,biao,s);
    toc
    disp(['运行时间: ', num2str(toc)]);
    tong(1,i)=toc;
    tong(2,i)=iter;
    betasi=find(b>mean);
    kf=intersect(a,betasi);
    tong(3,i)=length(m);
    tong(4,i)=length(kf);
    tong(5,i)= sqrt(sum((yorg- Xorg(:,betasi)*b(betasi)).^2)) ./ k;
    tong(6,i)=norm(b(kf)-btrue(kf), 1)/norm(btrue(kf), 1);
    tong(7,i)= sqrt(sum((yorg- Xorg*b).^2)) ./ p;
    tong(8,i)=norm(b-btrue, 1)/norm(btrue, 1);
    data1 = [data1 b];
end
s=1e1;
for i=l+1:2*l
    tic
    [b,iter,m,sigma1]= bhs_zu(Xorg, y_noisy,biao,s);
    toc
    disp(['运行时间: ', num2str(toc)]);
    tong(1,i)=toc;
    tong(2,i)=iter;
    betasi=find(b>mean);
    kf=intersect(a,betasi);
    tong(3,i)=length(m);
    tong(4,i)=length(kf);
    tong(5,i)= sqrt(sum((yorg- Xorg(:,betasi)*b(betasi)).^2)) ./ k;
    tong(6,i)=norm(b(kf)-btrue(kf), 1)/norm(btrue(kf), 1);
    tong(7,i)= sqrt(sum((yorg- Xorg*b).^2)) ./ p;
    tong(8,i)=norm(b-btrue, 1)/norm(btrue, 1);
    data1 = [data1 b];
end
s=1e2;
for i=2*l+1:3*l
    tic
    [b,iter,m,sigma1]= bhs_zu(Xorg, y_noisy,biao,s);
    toc
    disp(['运行时间: ', num2str(toc)]);
    tong(1,i)=toc;
    tong(2,i)=iter;
    betasi=find(b>mean);
    kf=intersect(a,betasi);
    tong(3,i)=length(m);
    tong(4,i)=length(kf);
    tong(5,i)= sqrt(sum((yorg- Xorg(:,betasi)*b(betasi)).^2)) ./ k;
    tong(6,i)=norm(b(kf)-btrue(kf), 1)/norm(btrue(kf), 1);
    tong(7,i)= sqrt(sum((yorg- Xorg*b).^2)) ./ p;
    tong(8,i)=norm(b-btrue, 1)/norm(btrue, 1);
    data1 = [data1 b];
end




%非子空间
for i=3*l+1:4*l
    tic
    [b,iter]= bhs_R(Xorg, y_noisy,biao);
    toc
    disp(['运行时间: ', num2str(toc)]);
    tong(1,i)=toc;
    tong(2,i)=iter;
    betasi=find(b>mean);
    kf=intersect(a,betasi);
    tong(3,i)=length(betasi);
    tong(4,i)=length(kf);
    tong(5,i)= sqrt(sum((yorg- Xorg(:,betasi)*b(betasi)).^2)) ./ k;
    tong(6,i)=norm(b(kf)-btrue(kf), 1)/norm(btrue(kf), 1);
    tong(7,i)= sqrt(sum((yorg- Xorg*b).^2)) ./ p;
    tong(8,i)=norm(b-btrue, 1)/norm(btrue, 1);
    data1 = [data1 b];
end

