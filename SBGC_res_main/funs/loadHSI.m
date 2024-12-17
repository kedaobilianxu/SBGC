function [data3D, lidar2D, gt, ind, c] = loadHSI(dataType)

switch dataType
    case 'Salinas'
        HSIdata = load('Salinas.mat');
        data3D = HSIdata.salinas;
        HSIlabel = load('Salinas_gt.mat');
        gt2D = HSIlabel.salinas_gt;
    case 'PaviaU'
       HSIdata = load('PaviaU.mat');
        data3D = HSIdata.paviaU;
        HSIlabel = load('PaviaU_gt.mat');
        gt2D = HSIlabel.paviaU_gt;
    case 'PaviaC'
        HSIdata = load('Pavia.mat');
        data3D = HSIdata.pavia;
        HSIlabel = load('Pavia_gt.mat');
        gt2D = HSIlabel.pavia_gt;
    case 'Trento'
        HSIdata = load('Trento-HSI.mat');
        Lidardata = load('Trento-Lidar.mat');
        data3D = HSIdata.HSI;
        lidar2D = Lidardata.Lidar(:,:,1);
        HSIlabel = load('Trento-GT.mat');
        gt2D = HSIlabel.GT;
    case 'MUUFL'
        HSIdata = load('HSI.mat');
        Lidardata = load('LiDAR_data_first_return.mat');
        data3D = HSIdata.HSI_data;
        lidar2D = Lidardata.LiDAR_data_first_return(:,:,1);
        HSIlabel = load('GT.mat');
        gt2D = HSIlabel.GT;
    case 'Houston'
        HSIdata = load('data_HS_LR.mat');
        Lidardata = load('data_LiDAR.mat');
        data3D = HSIdata.data_HS_LR;
        lidar2D = Lidardata.LiDAR;
        HSIlabel = load('GT-ALL.mat');
        gt2D = HSIlabel.gt;
    case 'Augsburg'
        HSIdata = load('data_HS_LR.mat');
        SARdata = load('data_SAR_HR.mat');
        DSMdata = load('data_DSM.mat');
        data3D = HSIdata.data_HS_LR;
        sar3D = SARdata.data_SAR_HR;
        dsm2D = DSMdata.data_DSM;
        HSIlabel = load('GT-ALL.mat');
        gt2D = HSIlabel.gt;
    otherwise
        error('No corresponding HSI dataset!');
end

gt = double(gt2D(:));
ind = find(gt);
c = length(unique(gt(ind)));

end