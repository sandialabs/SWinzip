% This code computes the global error norm between two functions
% represented on two different meshes while accounting for mode switching

clear all

load('../../data/mesh2D_irreg.mat');

vm=5;

% p1 is the larger mesh
p1=p;
x1=p1(1,:);
y1=p1(2,:);
N1=length(x1);

% p2 is smaller, scaled and shifted
p2=p1(:,1:3:end)*2+1.1;
x2=p2(1,:);
y2=p2(2,:);
N2=length(x2);

% computing wavelet matrices and normalized
cd ../src
V1 = build_alpert_matrix(p1,vm,-1)/sqrt(N1);
V2 = build_alpert_matrix(p2,vm,-1)/sqrt(N2);
cd ../examples

% computing the functions
f1 = 48*sin(6*pi.*x1).*sin(5*pi.*y1).*sin(4*pi.*x1);
g1 = 48*(sin(2.0*pi*x1) - sin(2.0*pi*y1)).*sin(2.0*pi*x1)-52;
h1 = 5*(9*sin(2.0*pi*x1) - 7*sin(2.0*pi*y1)).*sin(2.0*pi*x1)-52;

f2=f1(1:3:end);
g2=g1(1:3:end);
h2=h1(1:3:end);

% wavelet transforms
wf1=V1'*f1';
wg1=V1'*g1';
wh1=V1'*h1';

wf2=V2'*f2';
wg2=V2'*g2';
wh2=V2'*h2';

% sorting amplitudes
[wf1s,iif1]=sort(abs(wf1),'descend');
[wg1s,iig1]=sort(abs(wg1),'descend');
[wh1s,iih1]=sort(abs(wh1),'descend');

[wf2s,iif2]=sort(abs(wf2),'descend');
[wg2s,iig2]=sort(abs(wg2),'descend');
[wh2s,iih2]=sort(abs(wh2),'descend');

M=20;

subplot(1,2,1)
scatter(p1(1,:),p1(2,:),5,g1,'filled')
hold
scatter(p2(1,:),p2(2,:),5,g2,'filled')
xlabel('x');
ylabel('y');
xlim([0 3.1])
ylim([0 3.1])
subplot(1,2,2)
semilogy(1:M,wf1s(1:M),1:M,wf2s(1:M),1:M,wg1s(1:M),1:M,wg2s(1:M))
ylim([0.9 40])
xlabel('Wavelet mode');
ylabel('Wavelet transform amplitude');
legend('\phi_{f_1}','\phi_{f_2}','\phi_{g_1}','\phi_{g_2}','Interpreter','tex');

% if you want to compare g1 with h2 uncomment the following two lines
% iig2=iih2;
% wg2=wh2;

% the sorting indices will reveal mode switching
iigerr=iig2(1:M)+N1-N2-iig1(1:M);
iigerr2=iig2(1:M)+N1-N2;
iigerr1=iig1(1:M);

% finding the mode switching
im=1:M;
jj=find(iigerr~=0);
for i=1:length(jj)
    j=find(iigerr2(jj)-iigerr1(jj(i))==0);
    if isempty(j)
        break;
    else
        im(jj(i))=jj(j);
    end
end

% calculating norms with and without mode switching and the corresponding
% relative error
norm_not_fixed=norm(abs(wg1(iig1(1:M)))-abs(wg2(iig2(1:M))));
norm_fixed=norm(abs(wg1(iig1(1:M)))-abs(wg2(iig2(im))));
rel_err=abs(norm_not_fixed-norm_fixed)/norm_fixed
