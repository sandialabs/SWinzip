% This script compresses and reconstructs a dataset defined on a 2D complex
% geometry and unstructured mesh. Full wavelet details are considered but
% the reconstruction is performed in a continuous level of details (CLOD)
% fashion

clear all
x_I=[];activeSet=[];jj=1;

% Load the dataset
load('../../data/mesh2D_irreg.mat');
ne=size(t,2);
t=t(1:3,:);
x=p(1,:);
y=p(2,:);
N=length(x);
d=size(p,1);

% CLOD works best on data exhibiting many oscillations and features
f = (4.*sin(8*pi.*x) ).*(4.*sin(7*pi.*y) ).*3.*sin(6*pi.*x);
R=10;
vm=5;

cd ../src
rng(123);
fprintf('Computing sampling matrix\n');
P = constructSamplingMatrix( ceil(N/R), N, 'bernoulli');
y=P*f';


fprintf('Computing wavelet matrix product at level %i\n',1);
part_type = 'axis';
[V,Jmax,Ui,Mi] = build_alpert_matrix(p,vm,1,part_type);
cd ../examples

Ujj=V';
T=P;

if length(vm)==1
    vm = repmat(vm, 1, d);
end
cd ../src
[part,B,G] = dichotomic_grouping(p, prod(vm), [num2str(d) part_type]);
Jmax = log2(length(G))+1;

for kk=1:Jmax

    disp('Computing sampling-wavelet matrix product');
    T=T*Ujj';

    disp('Running StOMP');    
    tic
    [sol, numIters, I, x_I, activeSet] = SolveStOMP(T, y, 50, 0, 1e-8, x_I, activeSet);
    ts(kk)=toc;    
    
    fs=V*sol;
    es(kk)=sqrt(sum((f'-fs).^2)/N)/(max(f)-min(f));
    
    if (kk<Jmax)
        fprintf('Computing wavelet matrix product at level %i\n',kk+1);
        
        % Adding detail level to the wavelet matrix
        [V,Ui,Mi,Ujj] = addlevel_alpert_matrix(Ui,Mi,V',p,kk,1,prod(vm),Jmax); 
        w2=V'*fs;
        w1=zeros(N,1);w1(activeSet)=x_I;
        ma=sort(abs(x_I));        
        
        % Limiting the number of active wavelet coefficients to be used in
        % the next level reconstruction
        nn=ceil(0.0425*length(ma));
        activeSet=find(abs(w2)>ma(nn));
        x_I=w2(activeSet);
    end
    
end
cd ../examples

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

figure
subplot(1,2,1)
plot([1:Jmax],es)
xlabel('Detail level');
ylabel('NRMSE');

subplot(1,2,2)
plot([1:Jmax],ts,[1:Jmax],cumsum(ts))
xlabel('Detail level');
ylabel('Wall clock time (s)');
legend('Relative','Cumulative')