function [labels]= Kmeans_seg(data,cc)

[M,N,B]=size(data);
Y_scale=scaleForSVM(reshape(data,M*N,B));
p = 1;
[Y_pca] = pca(Y_scale, p);
K=cc;

labels = kmeans(Y_pca',K, 'MaxIter', 1000);
labels = labels  + 1;
labels = reshape(labels,M,N);