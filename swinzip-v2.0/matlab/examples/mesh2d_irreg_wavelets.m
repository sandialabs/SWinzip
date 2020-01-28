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
vm=5;                % wavelet order

% Perform all compression, wavelet computation and reconstruction
cd ../src

% Precompute mesh hierarchy and corresponding moment matrices
[Uj,part,k2] = compute_Uj(p,vm);

% Foreward wavelet transform
w = perform_alpert_transform(f,Uj,part,k2,1);

% Sort wavelet amplitudes and getting indices
[ws,i]=sort(abs(w));

% NR is the number of remaining values after compression
NR=floor((1-1/R)*N);

% Truncating coefficients of low magnitude
w(i(1:NR))=0;

% Inverse wavelet transform
fw = perform_alpert_transform(w,Uj,part,k2,-1);

% error
ew=sqrt(sum((f'-fw).^2)/N)/(max(f)-min(f));
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
trisurf(t',p(1,:),p(2,:),fw)
view(2)
shading interp
caxis([min(f) max(f)])
axis square
title('Data reconstructed from compressed samples');

fprintf('The normalized relative mean square error (NRMSE) \n of reconstruction from compressed samples: %f \n', ew);
fprintf('Compression ratio: %f \n', R);