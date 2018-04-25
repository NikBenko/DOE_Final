function [ denoised ] = OCT_Denoise( image_filename )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Load image, check size, convert to greyscale

if ischar(image_filename)==1
    image=imread(image_filename);
else
    image=image_filename;
end

[L,W,D]=size(image);

if D>1
    input=im2double(image(:,:,1));
else
    input=im2double(image);
end


% create a colormap
map = parula(1000); 
map = map(1:999,:);
map = [map;[1,1,1]];

%==========================================================================
% Statistical parameter estimation
%==========================================================================

%alpha = estimatePar(input);

%==========================================================================
% Run Proposed Method
%==========================================================================
if ~exist('alpha','var')
    % default value
    alpha = 0.35;
    %alpha = 1.0;
end
% compute distribution coefficient
c1 = (1-alpha^2/2)^(1/4);
c2 = 1-(1-alpha^2/2)^(1/2);

% parameter selection
par.lambda = 0.4;
par.gamma = 2;
par.delta = 0.002;
par.theta = 0.98;
par.c1 = c1;
par.c2 = c2;
par.maxIter = 20;

tic;
[ denoised ] = ladexp_huberTV( input, par );
toc;

%==========================================================================
% Display images
%==========================================================================
% figure; imshow(input,[],'Colormap',map); colorbar;
% title('Noisy');
% figure; imshow(denoised,[],'Colormap',map); colorbar;
% title('Denoised');
% 
figure;
subplot(1,2,1); imshow(input)
subplot(1,2,2); imshow(denoised)

%BW3=imbinarize(input,'adaptive');
BW1=imbinarize(denoised,'adaptive');
BW2=bwareaopen(BW1,150);
BW3=imfill(BW2,'holes');
figure;
subplot(1,3,1); imshow(BW1)
subplot(1,3,2); imshow(BW2)
subplot(1,3,3); imshow(BW3)

disp(alpha)

end
