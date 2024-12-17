function [labels]= LSCseg(data,cc)

[M,N,B]=size(data);
Y_scale=scaleForSVM(reshape(data,M*N,B));
p = 1;
[Y_pca] = pca(Y_scale, p);
% 
%  figure;
 Y_Pca_show = reshape(Y_pca,M,N);
%  imshow(Y_Pca_show, []);
% print('C:\Users\20355\Desktop\view\Trento-pca.png', '-dpng', '-r300');
img = im2uint8(mat2gray(reshape(Y_pca', M, N, p)));
img(:,:,1) = data(:,:,9);
img(:,:,2) = data(:,:,19);
img(:,:,3) = data(:,:,29);
K=cc;
% gaus=fspecial('gaussian',3);
% I=imfilter(img,gaus);
ratio=0.075;
[labels] = LSC_mex(img,K,ratio);
%labels = labels + 1; % Because the first segmentation label is zero.

end





