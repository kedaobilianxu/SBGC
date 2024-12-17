function [labels]= SLICseg(data,cc)

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
K=cc;

[labels, ~] = superpixels(double(img),K);
%labels = labels + 1; % Because the first segmentation label is zero.

end