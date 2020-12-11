% This script compresses and reconstructs a dataset defined on a 2D complex
% geometry and unstructured mesh. Full wavelet details are considered

clear all

% Load the dataset
load('../../data/mesh2D_irreg.mat');

% t is the connectivity matrix: nodes of each element
t=t(1:3,:);

% p contains the nodes coordinates in x and y
x=p(1,:);
y=p(2,:);
N=length(x);

% function to be compressed/reconstructed
% We have two functions

% The first highly oscillatory
f = (4.*sin(8*pi.*x) ).*(4.*sin(7*pi.*y) ).*3.*sin(6*pi.*x);

% The second is pretty smooth
% f = (4.0*sin(2.0*pi*x) - 4.0*sin(2.0*pi*y))*3.0.*sin(2.0*pi*x);

% R is the compression ratio
R=10;  % For the highly oscillatory function
% R=70;  % For the pretty smooth function (affords a much larger compression)

% setting the options
options.vm=5;                % wavelet order
options.mt='bernoulli';      % sampling matrix type
options.j=-1;                % wavelet detail level, -1 means full details
options.seed=123;            % random sampling seed
options.part_type='axis';    % domain subdivision method

% Perform all compression, wavelet computation and reconstruction
cd ../src
[y,fw,fs,ew,es,inc,x_I, activeSet,ts,wc,D]=CompressReconstruct(p, f, R ,options, [], []);
cd ../examples

% Plotting
subplot(1,2,1)
trisurf(t',p(1,:),p(2,:),f)
view(2)
shading interp
caxis([min(f) max(f)])
axis square
title('Original data field');

subplot(1,2,2)
trisurf(t',p(1,:),p(2,:),fs)
view(2)
shading interp
caxis([min(f) max(f)])
axis square
title('Data reconstructed from compressed samples');

% output the error metrics
disp('The error metrics of Kuan-Liu Ma:');
D

fprintf('The normalized relative mean square error (NRMSE) \n of reconstruction from compressed samples: %f \n', es);
fprintf('Compression ratio: %f \n', R);