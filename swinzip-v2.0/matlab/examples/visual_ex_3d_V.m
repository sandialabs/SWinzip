% This code compares two functions g and g1 represented on two different
% meshes and reconstructs the error field using wavelet transforms.

clear all

load('../../data/g3d.mat');

g = 48*(sin(2.0*pi*x) - sin(2.0*pi*y) - sin(2.5*pi*z)).*sin(2.0*pi*x).*sin(2.5*pi*z)-52;
g1 = 38*(1.1*sin(2.0*pi*x1) - 0.9*sin(2.0*pi*y1) - 1.05*sin(2.5*pi*z1)).*sin(2.0*pi*x1).*sin(2.5*pi*z1)-52;
cc=[min(g) max(g)];

N=length(x);
N1=length(x1);

% computing the error fields using the iterative and direct methods
Nd=8;
cd ../src
ee1=compute_error_field_alpert_dichotomic(p',p1',g',g1',Nd,0.01);
ee1V=compute_error_field_alpert_dichotomic_V(p',p1',g',g1',Nd);
cd ../examples

% computing the true error
eet1=-g1+(48*(sin(2.0*pi*x1) - sin(2.0*pi*y1) - sin(2.5*pi*z1)).*sin(2.0*pi*x1).*sin(2.5*pi*z1)-52);
eet1=eet1';
cce=[min(eet1) max(eet1)];

% plotting the true and estimated error fields
figure
subplot(1,3,1)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),eet1(ii),'filled');axis equal;caxis(cce);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');
xlabel('$$x$$','Interpreter','latex');
ylabel('$$y$$','Interpreter','latex');
zlabel('$$z$$','Interpreter','latex');
title('$$q_1 - q_2 $$ analytical','Interpreter','latex');
subplot(1,3,2)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),ee1V(ii),'filled');axis equal;caxis(cce);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');
xlabel('$$x$$','Interpreter','latex');
ylabel('$$y$$','Interpreter','latex');
zlabel('$$z$$','Interpreter','latex');
title('$$q_1 - q_2 $$ computed by the direct method','Interpreter','latex');
subplot(1,3,3)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),ee1(ii),'filled');axis equal;caxis(cce);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');
xlabel('$$x$$','Interpreter','latex');
ylabel('$$y$$','Interpreter','latex');
zlabel('$$z$$','Interpreter','latex');
title('$$q_1 - q_2 $$ computed by the iterative method','Interpreter','latex');
colorbar

norm(eet1-ee1)
norm(eet1-ee1V)
