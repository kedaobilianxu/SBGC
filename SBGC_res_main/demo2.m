clear; close all; clc;

addpath('./funs');
addpath('./Data/HSI-LiDAR-Trento');


dataType = 'Trento';   
%dataType = 'MUUFL'; 
%dataType = 'Houston'; 

%beta = [0.5:0.1:1];
%w_d = [10:10:50];
%M = [6:6:120];

 beta = [1];
 M = [24];
 w_d = [50];


%  beta = [0.5];
%  M = [44];
%  w_d = [40];
% 
% beta = [1];
% M = [15];
% w_d = [40];

row_index = 2;
[data3D, lidar2D, gt, ind, c] = loadHSI(dataType);

[nr,nc,~] = size(data3D);
for betai = 1:size(beta,2)
    fprintf('Beta : %f\n', beta(betai));
    for M_i = 1:size(M,2)
        for d_j = 1:size(w_d,2)
            tic;
            [y_pred,Z] = TBD2(data3D, lidar2D, gt, ind, beta(betai), M(M_i),c,w_d(d_j));
            running_time = toc;
            [result,PA, UA, AA, OA, Kappa] = HSI_ClusteringMeasure(gt(ind),y_pred(ind));
            results = [beta(betai), M(M_i), w_d(d_j),AA, OA, Kappa,running_time];
            label = ['A',num2str(row_index)];
            row_index = row_index + 1;
        end
    end
end

resall = y_pred;
res = resall(ind);
new_res = bestMap(gt(ind),y_pred(ind));
resallmap = zeros(nr*nc,1);
for i = 1:c
    temp = new_res(find(res==i));
    if ~isempty(temp)
        resallmap(find(resall==i)) = max(temp);
    else
        resallmap(find(resall==i)) = setdiff(gt(ind),new_res);
    end
end
%map = mat2rgb_Trento(resallmap,nr,nc); 
%map = mat2rgb_MUUFL(resallmap,nr,nc); 
%map = mat2rgb_Houston(resallmap,nr,nc);  % ground truth 
%figure(1)
%imshow(map,'border','tight','initialmagnification','fit')
%imshow(map,'border','tight','initialmagnification','fit')
%set(gcf,'Position',[0,0,1*nr,1*nc]);
%axis off
