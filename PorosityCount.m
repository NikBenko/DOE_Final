L=499;
manualCount=zeros(L,1);
auto2DCount=zeros(L,1);
auto3DCount=zeros(L,1);
erode2DCount=zeros(L,1);
erode3DCount=zeros(L,1);

topImage=imread('CGN_Pyrolysis_offcen_0978 _1.tif');
for i=1:L
    image1=imread('Manual_Segmentation_Stack.tif',i);
    image2=imread('Weka_2D_1_stack.tif',i);
    image3=imread('Weka_3D_1_stack.tif',i);
    bw1=imbinarize(image1,0);
    bw2=~imbinarize(image2,0);
    bw3=~imbinarize(image3,0); 
    
    se=strel('disk',1,4);
    erode1=imerode(bw2,se); 
    erode2=imerode(bw3,se); 
    
    manualCount(i)=sum(sum(bw1));
    auto2DCount(i)=sum(sum(bw2));
    auto3DCount(i)=sum(sum(bw3));
    erode2DCount(i)=sum(sum(erode1));
    erode3DCount(i)=sum(sum(erode2));
    
    while i==1
        topEdges=bw3-erode2;
        break
    end
    
        
end

% figure; imagesc(bw1)
% figure; imagesc(bw2)



manualPorosity=manualCount./(L^2);
auto2DPorosity=auto2DCount./(L^2);
auto3DPorosity=auto3DCount./(L^2);
erode2DPorosity=erode2DCount./(L^2);
erode3DPorosity=erode3DCount./(L^2);


%% Plots
figure;
hold on
plot(1:L,manualPorosity,1:L,auto2DPorosity,1:L,auto3DPorosity,1:L,erode2DPorosity,1:L,erode3DPorosity)
hold off
legend('Manual','Auto2D','Auto3D','Erode2D','Erode3D')
%% Stats

Data=[manualPorosity,auto2DPorosity,auto3DPorosity,erode2DPorosity,erode3DPorosity];
[p,tbl,stats] = anova1(Data); 
[c,m,h,gnames] = multcompare(stats);

Diff=(Data(:,2:end)-Data(:,1))./Data(:,1)*100;

%% Image Overlay
overlay=zeros(500,500,3);
overlay(:,:,1)=topImage;
overlay(:,:,2)=topImage;
overlay(:,:,3)=topImage;
red=topImage;
blue=topImage;
green=topImage;
red(topEdges==1)=255;
blue(topEdges==1)=0;
green(topEdges==1)=0;
overlay(:,:,1)=red;
overlay(:,:,2)=green;
overlay(:,:,3)=blue;
figure; imshow(uint8(overlay))

%%

% diff=autoPorosity-manualPorosity;
% pDiff=(autoPorosity-manualPorosity)./manualPorosity*100;
% 
% meanM=mean(manualPorosity);
% stdM=std(manualPorosity);
% meanA=mean(autoPorosity);
% stdA=std(autoPorosity);
% meadD=mean(diff);
% stdD=std(diff);
% 
% figure; plot(1:length(diff),pDiff)
% xlabel('Slice Number')
% ylabel('Percent Difference')
% set(gca,'FontSize',16)
% 
% figure; hist(diff)
