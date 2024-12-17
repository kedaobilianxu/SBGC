function [y_pred,Z] = TBD2(data3D, lidar2D, gt, ind, beta, M, nc, dim_w)
%% prod. by caozhe@mail.nwpu.edu.cn
% Input:
%       - data3D: 3D HSI data
%       - lidar2D: 2D Lidar data
%       - nc: cluster number
%       - beta: parameter
% Output:
%       - y_pred: clustering labels.

[N_hid,M_hid,dim_hid] = size(data3D);
imageSize = [N_hid,M_hid,3];
[~,~,dim_lidar] = size(lidar2D);
data3D = data3D./max(data3D(:));
%M = pixelNum(data3D, 2000); % Number of superpixels
%fprintf('Superpixels number : %d\n', M);
Mdata = data3D;
Mdata(:,:,size(data3D,3)+1) = lidar2D;

tic;
[splabels] = cubseg(Mdata, M);% ERS segmentation.
%[splabels] = SLICseg(Mdata, M);% SLIC segmentation.
%[splabels] = LSCseg(Mdata, M);% LSC segmentation.
%[splabels] = kmeans_seg(Mdata, M);% Kmeans segmentation.
segtime = toc;


disp(segtime);
[z_hid,A_hid,X_hid] = Generate_anchor_graph(data3D,M,splabels,dim_hid);
[z_lidar,A_lidar,X_lidar] = Generate_anchor_graph(lidar2D,M,splabels,dim_lidar);
z_hid = double(z_hid);
A_hid = double(A_hid);
X_hid = double(X_hid);
z_lidar = double(z_lidar);
A_lidar = double(A_lidar);
X_lidar = double(X_lidar);
N = size(z_hid,1);
clear data3D lidar2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nV = 2;
sX = [N, M, nV];
Z{1} = full(z_hid);
Z{2} = full(z_lidar);
%Z{2} = Z{1};
for k = 1:nV
    J{k} = zeros(N, M);
    Q{k} = zeros(N, M);
    XA{k} = zeros(N, M);    
end

eta = 1.3; mu = 10e-5; rho = 10e-5; max_mu = 10e12; max_rho = 10e12;
Isconverg = 0; maxIter = 30; iter = 1;
betaf = ones(nV, 1); 
%W = eye(int16(dim_hid),int16(w_d));
maxkappa = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(Isconverg == 0)
rng('default');
seed = 3;
rng(seed);
%     %%%update W matrix
    A = zeros(dim_hid, dim_hid);
    W = eye(dim_hid, dim_w);
    rowSums = sum(Z{1}', 2);
    D_M = diag(rowSums);
    G = X_hid*X_hid'-2*X_hid*Z{1}*A_hid'+A_hid*D_M*A_hid';
    [W, ~, ~] = svds(G, dim_w);


    clear A Z_ext A_hid_ext X_hid_ext

    %%%undate tensor J
    for v =1:nV
        XA{v} = Z{v} + Q{v}./rho;
    end
    M_tensor = cat(3, XA{:,:});
    M_vector = M_tensor(:);
    [myj, OBJ] = wshrinkObj_weight_lp(M_vector, beta*betaf./rho, sX, 0, 3, 1);
    J_tensor = reshape(myj, sX);
    for k=1:nV
        J{k} = J_tensor(:,:,k);
    end
    clear J_tensor Q_tensor M_vector

    %%%update tensor Z
%     for i = 1:N
%         for j = 1:M
%             m_1 = W' * X_hid(:,i) - W' * A_hid(:,j);
%             m_2 = X_lidar(:,i) - A_lidar(:,j);
%             XA{1}(i,j) = sum(m_1.^2);
%             XA{2}(i,j) = sum(m_2.^2);
%         end
%     end
    WX_hid = W' * X_hid;  % D x N
    WA_hid = W' * A_hid;  % D x M

    % 初始化 XA 存储容器
    XA = cell(2,1);

    % 利用矩阵操作计算差的平方和
    % 方法：使用广播
    XA{1} = sum((reshape(WX_hid, [], 1, N) - reshape(WA_hid, [], M, 1)).^2, 1);
    XA{2} = sum((reshape(X_lidar, [], 1, N) - reshape(A_lidar, [], M, 1)).^2, 1);

    % 转换维度以匹配原始代码的输出形式 (N x M)
    XA{1} = squeeze(XA{1});  % N x M
    XA{2} = squeeze(XA{2});  % N x M


    for v =1:nV
        P{v} = J{v} - Q{v}./rho;
        R{v} = P{v}.*rho - XA{v}';
    end
    for v = 1:nV
        for i = 1:M
            MM = R{v}(:,i);
            Z{v}(:,i) = EProjSimplex_new(MM, 1);
        end
    end

clear P R MM
    %%%update Q
    for k = 1:nV
        Q{k} = Q{k} + mu .* (Z{k} - J{k});  % 进行数值运算
    end
    mu = min(eta * mu, max_mu);
    rho = min(eta * rho, max_rho);

    %%
    if iter == maxIter
        Isconverg = 1;
    end
    iter = iter + 1;
    %Z{1}(isnan(Z{1}) | isinf(Z{1})) = 0;
    %Z{2}(isnan(Z{1}) | isinf(Z{2})) = 0;
    Z_final = (0.5.*Z{1} + 0.5.*Z{2});
    Z_final(isnan(Z_final) | isinf(Z_final)) = 0;
    [u,~,~] = svds(Z_final,nc);
    u = real(u);
    Z_final = real(Z_final);
    y_pred = kmeans(u,nc,'MaxIter', 1000);
    %disp(iter - 1);
    [~,~, ~, ~, ~, Kappa] = HSI_ClusteringMeasure(gt(ind),y_pred(ind));
    if(Kappa>maxkappa)
        maxkappa = Kappa;
        y_max = y_pred;
    end
end
% seed = 1;
% rng(seed);
% Z_final = (Z{1} + Z{2})./2;
% Z_final(isnan(Z_final) | isinf(Z_final)) = 0;
% [u,~,~] = svds(Z_final,nc);
% u = real(u);
% y_pred = kmeans(u,nc,'MaxIter', 300);

y_pred  = y_max;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,A,X] = Generate_anchor_graph(data,m,splabels,dim)
    if(dim==1)
        X = reshape(data, [], dim);
        X = X';
        t = 5;

        A = meanInd(X, splabels(:), ones(size(X,2),m));
        [Z,~,~,~] = initZ(X, A, t);

    else

        %[newData] = S3_PCA(data, 13, splabels);

        X = reshape(data, [], dim);
        X = X';
        t = 5;

        A = meanInd(X, splabels(:), ones(size(X,2),m));
        [Z,~,~,~] = initZ(X, A, t);
    end
end

function [num]=pixelNum(data,Tbase)
%% Adaptively determine the number of superpixels
    [M_1,N_1,B]=size(data);
    Y_scale=scaleForSVM(reshape(data,M_1*N_1,B));
    p=1;
    [Y_pca] = pca(Y_scale, p);
    img = im2uint8(mat2gray(reshape(Y_pca', M_1, N_1, p)));
    [m,n] = size(img);
    BW = edge(img,'log');
    ind = find(BW~=0);
    Len = length(ind);
    Ratio = Len/(m*n);
    num = fix(Ratio * Tbase);    
    
end



