close all
clear all

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
    bw1=im2bw(image1,0);
    bw2=~im2bw(image2,0);
    bw3=~im2bw(image3,0); 
    
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
        topSeg=bw3;
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
%% One Way Anova and Tukey's HSD

Data=[manualPorosity,auto2DPorosity,auto3DPorosity,erode2DPorosity,erode3DPorosity];


Diff=zeros(length(Data),4);
for i=1:4
    Diff(:,i)=(Data(:,i+1)-Data(:,1))./Data(:,1)*100;
end

[p,tbl,stats] = anova1(Data); 
[c,m,h,gnames] = multcompare(stats);
%% 2 Way Anova
group1=cell(4*L,1);
group2=cell(4*L,1);
for ii=1:L
    group1{ii}='2D';
    group1{ii+L}='3D';
    group1{ii+2*L}='2D';
    group1{ii+3*L}='3D';
    group2{ii}='Lo';
    group2{ii+L}='Lo';
    group2{ii+2*L}='Hi';
    group2{ii+3*L}='Hi';    
end
anovan(Diff(:)/100,{group1,group2},'model','full')

%% 2by2 Factorial
A=[1 -1 -1 1;
    1 1 -1 -1;
    1 -1 1 -1;
    1 1 1 1];
B=mean(Diff)';
c=A\B;

x=-1:0.1:1; y=x;
[X,Y]=meshgrid(x,y);
f=c(1)+c(2)*X+c(3)*Y+c(4)*X.*Y;
figure;
%[C,h]=contourf(X,Y,f);
% surf(X,Y,f,'FaceColor','interp',...
%     'EdgeColor','none',...
%     'FaceLighting','phong')
surf(X,Y,f)
daspect([1 1 100])
axis tight
view(+115,30)
camlight left)
xlabel('Analysis Demension')
ylabel('Erosion')
d=colorbar;
d.Label.String='% Error';
d.FontSize=18;
set(gca,'FontSize',18)
surfax=gca;
surfax.XTick=[-1 1];
surfax.YTick=[-1 1];
surfax.YTickLabel={'None','1 Pixel'};
surfax.XTickLabel={'2D','3D'};
%anovan(Diff,{'2D','3D','2D erode','3D erode'})
%% Image Overlay
overlay=zeros(500,500,3);
overlay(:,:,1)=topImage;
overlay(:,:,2)=topImage;
overlay(:,:,3)=topImage;
red=topImage;
blue=topImage;
green=topImage;
red(topEdges==1)=255;
blue(topSeg==1)=255;
green(topEdges==1)=0;
overlay(:,:,1)=red;
overlay(:,:,2)=green;
overlay(:,:,3)=blue;
figure; imshow(uint8(overlay))

%% Histograms
figure;
h1=histogram(Data(:,1));
hold on
h2=histogram(Data(:,2));
h3=histogram(Data(:,3));
h4=histogram(Data(:,4));
h5=histogram(Data(:,5));

h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h4.Normalization = 'probability';
h5.Normalization = 'probability';

xlabel('Porosity')
ylabel('Probability')
set(gca,'FontSize',16)
legend('Manual','Auto2D','Auto3D','Erode2D','Erode3D')
%% Diff Histograms
figure;
h1=histogram(Diff(:,1));
hold on
h2=histogram(Diff(:,2));
h3=histogram(Diff(:,3));
h4=histogram(Diff(:,4));


h1.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
h4.Normalization = 'probability';


xlabel('% Error')
ylabel('Probability')
set(gca,'FontSize',16)
legend('Auto2D','Auto3D','Erode2D','Erode3D')
title('Difference Histograms')
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
