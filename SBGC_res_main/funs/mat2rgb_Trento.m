% %------------------------Rong Wang 2015.09.30------------------------------
% function g = mat2rgb_Trento(result,m,n)
% 
% %---------------------------Background-------------------------------------
% cc(result==0,1,1)=0;  
% cc(result==0,1,2)=0;
% cc(result==0,1,3)=0;
% %---------------------------Asphalt(6631)----------------------------------
% cc(result==1,1,1)=196;  
% cc(result==1,1,2)=196;
% cc(result==1,1,3)=196;
% %-----------------------Meadow(18649)--------------------------------------
% cc(result==2,1,1)=0;
% cc(result==2,1,2)=255;
% cc(result==2,1,3)=0;
% %------------------------Gravel(2099)--------------------------------------
% cc(result==3,1,1)=0;
% cc(result==3,1,2)=255;
% cc(result==3,1,3)=255;
% %---------------------------Trees(3064)------------------------------------
% cc(result==4,1,1)=128;
% cc(result==4,1,2)=0;
% cc(result==4,1,3)=128;
% %-------------------Metal Sheets(1345)-------------------------------------
% cc(result==5,1,1)=0;
% cc(result==5,1,2)=128;
% cc(result==5,1,3)=0;
% %--------------------------Soil(5029)--------------------------------------
% cc(result==6,1,1)=165;
% cc(result==6,1,2)=82;
% cc(result==6,1,3)=40;
% %-----------------------------Bitumen(1330)--------------------------------
% cc(result==7,1,1)=128;
% cc(result==7,1,2)=0;
% cc(result==7,1,3)=128;
% %-------------------------Bricks(3682)-------------------------------------
% cc(result==8,1,1)=255;
% cc(result==8,1,2)=0;
% cc(result==8,1,3)=255;
% %----------------------------Shadows(947)----------------------------------
% cc(result==9,1,1)=255;
% cc(result==9,1,2)=255;
% cc(result==9,1,3)=0;
% 
% cc2=reshape(cc,m,n,3);
% cc2=uint8(cc2);
% g=[];g=uint8(g);
% g(:,:,1)=cc2(:,:,1);
% g(:,:,2)=cc2(:,:,2);
% g(:,:,3)=cc2(:,:,3);
%------------------------Zhe Cao 2024.08.25------------------------------
%color of Trento
function g = mat2rgb_Trento(result,m,n)

%---------------------------Background-------------------------------------
cc(result==0,1,1)=0;  
cc(result==0,1,2)=0;
cc(result==0,1,3)=0;
%----------------------------Apple trees(4034)-----------------------------
cc(result==1,1,1)=0;
cc(result==1,1,2)=100;
cc(result==1,1,3)=0;
%----------------------------Buildings(2903)-------------------------------
cc(result==2,1,1)=196;  
cc(result==2,1,2)=196;
cc(result==2,1,3)=196;
%----------------------------Ground(479)-----------------------------------
cc(result==3,1,1)=255;
cc(result==3,1,2)=255;
cc(result==3,1,3)=0;
%----------------------------Wood(9123)------------------------------------
cc(result==4,1,1)=165;
cc(result==4,1,2)=82;
cc(result==4,1,3)=40;
%----------------------------Vineyard(10501)-------------------------------
cc(result==5,1,1)=128;
cc(result==5,1,2)=0;
cc(result==5,1,3)=128;
%----------------------------Roads(3174)------------------------------------
cc(result==6,1,1)=0;
cc(result==6,1,2)=255;
cc(result==6,1,3)=0;

cc2=reshape(cc,m,n,3);
cc2=uint8(cc2);
g=[];g=uint8(g);
g(:,:,1)=cc2(:,:,1);
g(:,:,2)=cc2(:,:,2);
g(:,:,3)=cc2(:,:,3);
