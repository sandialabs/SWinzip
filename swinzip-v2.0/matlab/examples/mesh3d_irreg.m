% This script compresses and reconstructs a dataset defined on a 3D complex
% geometry and unstructured mesh. Full wavelet details are considered

clear all

% Load the dataset
load('../../data/mesh3D_irreg.mat');
p(3,:)=p(3,:)-0.5; % Centering the domain around the origin

% p contains the nodes coordinates in x and y
x=p(1,:);
y=p(2,:);
z=p(3,:);
N=length(x);

% function to be compressed/reconstructed
% We have two functions

% The first highly oscillatory
% f = (4.0*sin(6.0*pi*x) - 3.0*sin(5.0*pi*z))*2.0.*(sin(4.0*pi*y) - 2.0*sin(3.0*pi*z));
% The second is pretty smooth
f = (3.0*sin(2.0*pi*x) - 4.0*sin(2.0*pi*z))*2.0.*(2.0*sin(2.0*pi*y) - sin(2.0*pi*z));

% R is the compression ratio
% R=15;  % For the highly oscillatory function
R=100;  % For the pretty smooth function

% setting the options
options.vm=3;
options.mt='bernoulli';
options.j=-1;
options.seed=123;
options.part_type='axis';

% Perform all compression, wavelet computation and reconstruction
cd ../src
[y,fw,fs,ew,es,inc,x_I, activeSet,ts,wc,D]=CompressReconstruct(p, f, R ,options, [], []);
cd ../examples

% Plotting
I=find(p(3,:)>0);
p(3,I)=p(3,I)+0.25;
I=find(p(3,:)<=0);
p(3,I)=p(3,I)-0.25;

hFig = figure(1);
set(hFig, 'Position', [100 100 1800 800]);

S = repmat([5],size(p,2),1);
subplot(1,2,1)
scatter3(p(1,:),p(2,:),p(3,:),S,f,'filled');
caxis([min(f) max(f)])
xlim([min(p(1,:)) max(p(1,:))])
ylim([min(p(2,:)) max(p(2,:))])
zlim([min(p(3,:)) max(p(3,:))])
title('Original data field');

subplot(1,2,2)
scatter3(p(1,:),p(2,:),p(3,:),S,fs,'filled');
caxis([min(f) max(f)])
xlim([min(p(1,:)) max(p(1,:))])
ylim([min(p(2,:)) max(p(2,:))])
zlim([min(p(3,:)) max(p(3,:))])
title('Data reconstructed from compressed samples');

% output the error metrics
disp('The error metrics of Kuan-Liu Ma:');
D

fprintf('The normalized relative mean square error (NRMSE) \n of reconstruction from compressed samples: %f \n', es);
fprintf('Compression ratio: %f \n', R);