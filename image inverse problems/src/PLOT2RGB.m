function PLOT2RGB(mask_size,mask_deviation,standard_deviation,Img,blured_Img,x_SPA,x_SGSA_FFT,x_SB,name)
Img=uint8(Img);
blured_Img=uint8(blured_Img);
x_SPA=uint8(x_SPA);
x_SGSA_FFT=uint8(x_SGSA_FFT);
x_SB=uint8(x_SB);
figure;
subplot(2,3,1.5);
imshow(Img,[0 255]);
axis equal; axis off;
title(['Original image',num2str(mask_size),' ',num2str(mask_deviation),' ',num2str(standard_deviation)]);

subplot(2,3,2.5);
imshow(blured_Img, [0 255]);
axis equal; axis off;
title(['Blurred and noisy observation','  PSNR: ',num2str(psnr(blured_Img, Img)),'  SSIM:',num2str(ssim(blured_Img, Img))]);
imwrite(blured_Img,['Blurred',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);


subplot(2,3,4);
imagesc(x_SPA,[0 255]);
axis equal; axis off;
title(['xSPA','  PSNR: ',num2str(psnr(x_SPA, Img)),'  SSIM:',num2str(ssim(x_SPA, Img))]);
imwrite(x_SPA,['x_SPA',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

subplot(2,3,5);
imagesc(x_SGSA_FFT,[0 255]);
axis equal; axis off;
title(['xSGSAFFT','  PSNR: ',num2str(psnr(x_SGSA_FFT, Img)),'  SSIM:',num2str(ssim(x_SGSA_FFT, Img))]);
imwrite(x_SGSA_FFT,['x_SGSA_FFT',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

subplot(2,3,6);
imagesc(x_SB,[0 255]);
axis equal; axis off;
title(['xSB','  PSNR: ',num2str(psnr(x_SB, Img)),'  SSIM:',num2str(ssim(x_SB, Img))]);
imwrite(x_SB,['x_SB',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_',name]);

figure;
imshow(Img);
print(['RectangleImg',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(a)1.eps'],'-depsc','-r1200');%

imshow(blured_Img);
print(['RectangleBlurred',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(b)1.eps'],'-depsc','-r1200');

imshow(x_SPA);
print(['RectangleSPA',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(c)1.eps'],'-depsc','-r1200');

imshow(x_SGSA_FFT);
print(['RectangleSGSA_FFT',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(d)1.eps'],'-depsc','-r1200');

imshow(x_SB);
print(['RectangleSB',num2str(mask_size),'_',num2str(mask_deviation),'_',num2str(standard_deviation),'_Fig.9(e)1.eps'],'-depsc','-r1200');