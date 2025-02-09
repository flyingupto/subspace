function PLOTRGB(p1,p2,N,mask_size,mask_deviation,standard_deviation,Img,blured_Img,x_SPA,x_SGSA_FFT,x_SB,name)
Img=uint8(Img);
blured_Img=uint8(blured_Img);
x_SPA=uint8(x_SPA);
x_SGSA_FFT=uint8(x_SGSA_FFT);
x_SB=uint8(x_SB);

figure;
subplot(2,3,1.5);
imshow(Rectangle(Img,p1,p2,N),[0 255]);
axis equal; axis off;
title(['Original image',num2str(mask_size),' ',num2str(mask_deviation),' ',num2str(standard_deviation)]);

subplot(2,3,2.5);
imshow(Rectangle(blured_Img,p1,p2,N), [0 255]);
axis equal; axis off;
title(['Blurred and noisy observation','  PSNR: ',num2str(psnr(blured_Img, Img)),'  SSIM:',num2str(ssim(blured_Img, Img))]);
imwrite(blured_Img,['Blurred',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);


subplot(2,3,4);
imagesc(Rectangle(x_SPA,p1,p2,N),[0 255]);
axis equal; axis off;
title(['xSPA','  PSNR: ',num2str(psnr(x_SPA, Img)),'  SSIM:',num2str(ssim(x_SPA, Img))]);
imwrite(x_SPA,['x_SPA',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

subplot(2,3,5);
imagesc(Rectangle(x_SGSA_FFT,p1,p2,N),[0 255]);
axis equal; axis off;
title(['xSGSAFFT','  PSNR: ',num2str(psnr(x_SGSA_FFT, Img)),'  SSIM:',num2str(ssim(x_SGSA_FFT, Img))]);
imwrite(x_SGSA_FFT,['x_SGSA_FFT',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

subplot(2,3,6);
imagesc(Rectangle(x_SB,p1,p2,N),[0 255]);
axis equal; axis off;
title(['xSB','  PSNR: ',num2str(psnr(x_SB, Img)),'  SSIM:',num2str(ssim(x_SB, Img))]);
imwrite(x_SB,['x_SB',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

figure;
imshow(Rectangle(Img,p1,p2,N));
print(['RectangleImg',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(a)1.eps'],'-depsc','-r1200');%

imshow(Rectangle(blured_Img,p1,p2,N));
print(['RectangleBlurred',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(b)1.eps'],'-depsc','-r1200');

imshow(Rectangle(x_SPA,p1,p2,N));
print(['RectangleSPA',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(c)1.eps'],'-depsc','-r1200');

imshow(Rectangle(x_SGSA_FFT,p1,p2,N));
print(['RectangleSGSA_FFT',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(d)1.eps'],'-depsc','-r1200');

imshow(Rectangle(x_SB,p1,p2,N));
print(['RectangleSB',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(e)1.eps'],'-depsc','-r1200');
%-------------------------------------------------------------------------%
function img=Rectangle(img,p1,p2,N)
size1 = 60;
size2 = N - 270;

% 选取放大部分
part = img(p1+1:p1+size1, p2+1:p2+size1, :);

% 放大后的大小
mask = imresize(part, [size2, size2]);

% 放置位置
img(1:size2, 1:size2, :) = mask;

% 画框
img = insertShape(img, 'Rectangle', [p2, p1, size1, size1], 'Color', 'red', 'LineWidth', 1);
img = insertShape(img, 'Rectangle', [1, 1, size2, size2], 'Color', 'black', 'LineWidth', 1);
