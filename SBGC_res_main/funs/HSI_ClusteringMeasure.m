function [result,PA, UA, AA, OA, Kappa] = HSI_ClusteringMeasure(Y,predY)
%   ��������ڷ������Բ����������з��࣬������
%   �������� cf_matrix, ƽ�����ྫ�� Average classification Accuracy (AA)
%   ������ྫ�� Overall classification Accuracy (OA), Kappa coefficient (kappa)

predY = bestMap(Y,predY);

cf_matrix = confusionmat(Y,predY); % ��������

N = sum(sum(cf_matrix));    % �������������
col_num = sum(cf_matrix,1); % �������������������������
row_num = sum(cf_matrix,2); % �������������������������

PA = diag(cf_matrix)./col_num';  % �����߾��� 

UA = diag(cf_matrix)./row_num;  % ʹ���߾��� 

AA = sum(UA)/size(UA,1);

OA = trace(cf_matrix)/N;  % �ܷ��ྫ��

Kappa = (N*trace(cf_matrix)-col_num*row_num)/(N^2-col_num*row_num);

result = [PA;UA;AA;OA;Kappa];

end



