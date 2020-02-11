% This code compares two functions q1 and q2 represented on two different
% meshes and plots wavelet coefficients and amplitudes

clear all

load('../../data/g3d.mat');

q1 = 48*(sin(2.0*pi*x) - sin(2.0*pi*y) - sin(2.5*pi*z)).*sin(2.0*pi*x).*sin(2.5*pi*z)-52;
q2 = 38*(1.1*sin(2.0*pi*x1) - 0.9*sin(2.0*pi*y1) - 1.05*sin(2.5*pi*z1)).*sin(2.0*pi*x1).*sin(2.5*pi*z1)-52;

cc=[min(g) max(g)];

% Plotting the fields g and g1 on the meshes (x,y,z) and (x1,y1,z1), respectively
figure
subplot(2,2,1)
ii=find(p(1,:)>=p(3,:));scatter3(p(1,ii),p(2,ii),p(3,ii),25*ones(length(ii),1),q1(ii),'filled');axis equal;caxis(cc);
hold
ii=find(p(1,b)>=p(3,b));scatter3(p(1,b(ii)),p(2,b(ii)),p(3,b(ii)),'k.');
title('q_1 plotted on Mesh 1','Interpreter','tex');
subplot(2,2,2)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),q2(ii),'filled');axis equal;caxis(cc);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');
colorbar
title('q_2 plotted on Mesh 2','Interpreter','tex');

N=length(x);
N1=length(x1);

% M is the number of wavelet coefficients to retain when plotting the wavelet
% coefficients below
M=100;


% vm is the wavelet order
vm=3;

% Building wavelet matrices for the two different meshes
disp('Building wavelet matrices');
cd ../src
V = build_alpert_matrix(p,vm,-1);
V1 = build_alpert_matrix(p1,vm,-1);
cd ../examples

mxmn=max(g)-min(g);
mxmn1=max(g1)-min(g1);

w=V'*q1';
w1=V1'*q2';
 
[ws,ii]=sort(abs(w),'descend');
[w1s,ii1]=sort(abs(w1),'descend');

% plotting the first M sorted wavelet coefficient amplitudes
subplot(2,2,4)
semilogy(ws(1:M)/sqrt(N))
hold
semilogy(w1s(1:M)/sqrt(N1))
xlabel('Wavelet mode');
ylabel('Wavelet transform amplitude');
ylim([0.6 100])
legend('\phi_{q_1}','\phi_{q_2}');

% plotting the first M sorted wavelet coefficient in raw form
subplot(2,2,3)
plot(w(ii(1:M))/sqrt(N))
hold
plot(w1(ii1(1:M))/sqrt(N1))
xlabel('Wavelet mode');
ylabel('Wavelet transform');
legend('w_{q_1}','w_{q_2}');

