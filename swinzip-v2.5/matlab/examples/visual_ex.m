% This code compares three functions f, g and h represented on the same
% mesh by computing a wavelet-based global error norms

clear all

load('../../data/mesh2D_irreg.mat');

% t is the connectivity matrix: nodes of each element
t=t(1:3,:);

M=20;

% p contains the nodes coordinates in x and y
x=p(1,:);
y=p(2,:);
N=length(x);

f = 48*sin(6*pi.*x).*sin(5*pi.*y).*sin(4*pi.*x);
g = 48*(sin(2.0*pi*x) - sin(2.0*pi*y)).*sin(2.0*pi*x)-52;
h = 5*(9*sin(2.0*pi*x) - 7*sin(2.0*pi*y)).*sin(2.0*pi*x)-52;

% computing the wavelet matrix
cd ../src
V = build_alpert_matrix(p,5,-1);
cd ../examples

mxmnf=max(f)-min(f);
mxmng=max(g)-min(g);
mxmnh=max(h)-min(h);

% performing the wavelet transforms
wf=V'*f';
wg=V'*g';
wh=V'*h';

% sorting the amplitudes
[wfs,iif]=sort(abs(wf),'descend');
[wgs,iig]=sort(abs(wg),'descend');
[whs,iih]=sort(abs(wh),'descend');

% plotting the functions
figure
subplot(1,3,1)
trisurf(t',p(1,:),p(2,:),f);cc=[min(g) max(g)];
view(2)
shading interp
axis square
title('f');
caxis(cc);

subplot(1,3,2)
trisurf(t',p(1,:),p(2,:),g)
view(2)
shading interp
axis square
title('g');
caxis(cc);

subplot(1,3,3)
trisurf(t',p(1,:),p(2,:),h)
view(2)
shading interp
axis square
title('h');
caxis(cc);
colorbar

% plotting the wavelet amplitudes and coefficients
figure
subplot(1,2,2)
semilogy(wfs(1:M))
hold
semilogy(wgs(1:M))
semilogy(whs(1:M))
ylim([1e2 1e4])
xlim([0 20])
xlabel('Wavelet mode');
ylabel('Wavelet transform amplitude \phi','Interpreter','tex');
legend('\phi_f','\phi_g','\phi_h');

subplot(1,2,1)
plot(wf(iif(1:M)))
hold
plot(wg(iig(1:M)))
plot(wh(iih(1:M)))
ylim([-3000 6000])
xlim([0 20])
xlabel('Wavelet mode');
ylabel('Wavelet transform w');
legend('w_f','w_g','w_h');

MM=N;
ef=zeros(MM,1);eg=ef;eh=ef;
ewfg=ef;ewgh=ef;
lewfg=ef;lewgh=ef;
lrewfg=ef;lrewgh=ef;
for M=1:MM
    % WARNING: computing ef, eg and eh for all values of M in this loop
    % is time consuming (see below)
%     ef(M)=norm(f'-V(:,iif(1:M))*wf(iif(1:M)))/sqrt(N)/mxmnf;
%     eg(M)=norm(g'-V(:,iig(1:M))*wg(iig(1:M)))/sqrt(N)/mxmng;
%     eh(M)=norm(h'-V(:,iih(1:M))*wh(iih(1:M)))/sqrt(N)/mxmnh;

    ewfg(M)=norm(wfs(1:M)-wgs(1:M))/sqrt(M)/wgs(1);
    ewgh(M)=norm(wgs(1:M)-whs(1:M))/sqrt(M)/wgs(1);
    
    lewfg(M)=norm(log(wfs(1:M)/sqrt(N))-log(wgs(1:M)/sqrt(N)))/sqrt(M)/log(wgs(1)/sqrt(N));
    lewgh(M)=norm(log(wgs(1:M)/sqrt(N))-log(whs(1:M)/sqrt(N)))/sqrt(M)/log(wgs(1)/sqrt(N));
    
    lrewfg(M)=mean(abs(log(wfs(1:M)/sqrt(N))-log(wgs(1:M)/sqrt(N)))./abs(log(wgs(1:M)/sqrt(N))));
    lrewgh(M)=mean(abs(log(wgs(1:M)/sqrt(N))-log(whs(1:M)/sqrt(N)))./abs(log(wgs(1:M)/sqrt(N))));
end

% Instead of computing ef, eg and eh in the previous loop, they can be read
% as pre-computed from this file 
load ('res_visual_ex.mat','ef','eg','eh');

% plotting global error norms
figure
subplot(1,3,1)
semilogy(1:MM,ef,1:MM,eg,1:MM,eh)
ylabel('NRMSE');
legend('f','g','h');
ylim([1e-8 1])
subplot(1,3,2)
semilogy(1:MM,lewfg,1:MM,lewgh)
ylim([1e-3 100])
ylabel('\epsilon');
subplot(1,3,3)
semilogy(1:MM,lrewfg,1:MM,lrewgh)
ylim([1e-3 100])
ylabel('\rhi');
