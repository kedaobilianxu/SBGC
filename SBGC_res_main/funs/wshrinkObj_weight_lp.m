function [x,objV] = wshrinkObj_weight_lp(x, rho, sX, isWeight, mode, p)
%����ԭ�������⣬�ļ�: wshrinkObj_weight_lp.m ��: 12 ��: 1
%rho��һ�������������ĵ�i��Ԫ�ض�Ӧ���ǵ�i������ֵ��Ȩ��
% if find(isnan(x(:)))
%     disp('error')
% end
if isWeight == 1
    % mode = 3�ǲ���top slice�ķ���
    mode = 1;
end
%%�˴������޸�
if ~exist('mode','var')
    % mode = 1�ǲ���lateral slice�ķ���
    % mode = 2�ǲ���front slice�ķ���
    %     C = 2*sqrt(2)*sqrt(sX(3)*sX(2));
    C = sqrt(sX(3)*sX(2));
end
X=reshape(x,sX);
% if find(isnan(X))
%     disp('error')
% end
% �����׵ľ�����һ����ģ���������ת

if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);
    %Y=permute(X, [1, 3, 2]);
else
    Y = X;
end
% if find(isnan(Y))
%     disp('error')
% end
Yhat = fft(Y,[],3);
objV = 0;

% ȡһ��ģ�ķ���ĸ���
if mode == 1
    n3 = sX(2);
elseif mode == 3
    n3 = sX(1);
else
    n3 = sX(3);
end

if isinteger(n3/2)
    endValue = int16(n3/2+1);
    for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        
        if isWeight
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau = rho;
            shat=diag(shat);
            shat = solve_Lp_w(shat, tau, p);
            shat=diag(shat);
        end
        
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
            objV = objV + sum(shat(:));
        end
    end
    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');
    if isWeight
        weight = C./(diag(shat) + eps);
        tau = rho*weight;
        shat = soft(shat,diag(tau));
    else
        tau = rho;
        shat=diag(shat);
        shat = solve_Lp_w(shat, tau, p);
        shat=diag(shat);
    end
    objV = objV + sum(shat(:));
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
else
    %L = 10;
    %Y_ini = zeros(size(Yhat));
    %Y_ini(:,:,1:L) = Yhat(:,:,1:L);
    %Y_ini(:,:,n3-L+1:n3) = Yhat(:,:,n3-L+1:n3);
    %Yhat = Y_ini;
    %J = zeros(size(Yhat));
     endValue = int16(n3/2+1);
     for i = 1:endValue
        [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
        if isWeight
            weight = C./(diag(shat) + eps);
            tau = rho*weight;
            shat = soft(shat,diag(tau));
        else
            tau=rho;
            shat=diag(shat);
            shat = solve_Lp_w(shat, tau, p);
            shat=diag(shat);
        end
        objV = objV + sum(shat(:));
        Yhat(:,:,i) = uhat*shat*vhat';
        if i > 1
            Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
            objV = objV + sum(shat(:));
        end
%         if(i == 1)
%             J(:,:,1) = Yhat(:,:,1);
%         else
%             J(:,:,i) = Yhat(:,:,i);
%             J(:,:,n3-i+2) = conj(J(:,:,i));
    
    end
end

%% ԭ������仯
Y = ifft(Yhat,[],3);
%Y = ifft(J,[],3);
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X=shiftdim(Y, 2);
    %X=permute(Y, [1, 3, 2]);
else
    X = Y;
end

x = X(:);

end
