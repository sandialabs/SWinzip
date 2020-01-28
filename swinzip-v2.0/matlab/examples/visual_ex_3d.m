% This code compares two functions g and g1 represented on two different
% meshes and reconstructs the error field usign wavelet transforms.
% It also computes a wavlet-based global norm of the error

clear all

load('../../data/g3d.mat');

g = 48*(sin(2.0*pi*x) - sin(2.0*pi*y) - sin(2.5*pi*z)).*sin(2.0*pi*x).*sin(2.5*pi*z)-52;
g1 = 38*(1.1*sin(2.0*pi*x1) - 0.9*sin(2.0*pi*y1) - 1.05*sin(2.5*pi*z1)).*sin(2.0*pi*x1).*sin(2.5*pi*z1)-52;

cc=[min(g) max(g)];

% Plotting the fields g and g1 on the meshes (x,y,z) and (x1,y1,z1), respectively
figure
subplot(1,2,1)
ii=find(p(1,:)>=p(3,:));scatter3(p(1,ii),p(2,ii),p(3,ii),25*ones(length(ii),1),g(ii),'filled');axis equal;caxis(cc);
hold
ii=find(p(1,b)>=p(3,b));scatter3(p(1,b(ii)),p(2,b(ii)),p(3,b(ii)),'k.');
subplot(1,2,2)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),g1(ii),'filled');axis equal;caxis(cc);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');

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

w=V'*g';
w1=V1'*g1';
 
[ws,ii]=sort(abs(w),'descend');
[w1s,ii1]=sort(abs(w1),'descend');

% plotting the first M sorted wavelet coefficient amplitudes
figure
semilogy(ws(1:M)/sqrt(N))
hold
semilogy(w1s(1:M)/sqrt(N1))
xlabel('Wavelet mode');
ylabel('Wavelet transform amplitude');

% plotting the first M sorted wavelet coefficient in raw form
figure
plot(w(ii(1:M))/sqrt(N))
hold
plot(w1(ii1(1:M))/sqrt(N1))
xlabel('Wavelet mode');
ylabel('Wavelet transform');


MM=min(N1,N);
e=zeros(MM,1);e1=e;
ewgg1=e;
lewgg1=e;
lrewgg1=e;

% Computing global error metrics as a function of the wavelet modes
disp('Computing global error metrics');
for M=1:MM
%     if mod(M,100)==0
%        M
%     end

    % WARNING: computing e and e1 for all values of M in this loop
    % is time consuming (see below)
%     e(M)=norm(g'-V(:,ii(1:M))*w(ii(1:M)))/sqrt(N)/mxmn;
%     e1(M)=norm(g1'-V1(:,ii1(1:M))*w1(ii1(1:M)))/sqrt(N1)/mxmn1;
    
    ewgg1(M)=norm(ws(1:M)-w1s(1:M)*sqrt(N/N1))/sqrt(M)/ws(1);
        
    lewgg1(M)=norm(log(ws(1:M))-log(w1s(1:M))*sqrt(N/N1))/sqrt(M)/log(ws(1));
        
    lrewgg1(M)=mean(abs(log(ws(1:M))-log(w1s(1:M))*sqrt(N/N1))./abs(log(ws(1:M))));    
end

% Instead of computing e and e1 in the previous loop, they can be read
% as pre-computed from this file 
load ('res_visual_ex_3d.mat','e','e1');

% plotting global error metrics as a function of the retained number of
% wavelet coefficients
figure
subplot(1,4,1)
semilogy(1:MM,e,1:MM,e1)
title({'NRMSE','Reconstruction error'});
subplot(1,4,2)
semilogy(1:MM,ewgg1)
title({'$$L_2$$ error between','the wavelet modes'});
subplot(1,4,3)
semilogy(1:MM,lewgg1)
title({'$$L_2$$ error between the','log of the wavelet modes'});
subplot(1,4,4)
semilogy(1:MM,lrewgg1)
title({'Average relative error between',' the log of the wavelet modes'});

% Computing the error field
disp('Computing error field');
cd ../src
ee1=compute_error_field_alpert(p',p1',g',g1',4,4,4,0.001,vm);
cd ../examples

% this is the true analytical error field
eet1=-g1+(48*(sin(2.0*pi*x1) - sin(2.0*pi*y1) - sin(2.5*pi*z1)).*sin(2.0*pi*x1).*sin(2.5*pi*z1)-52);
cce=[min(eet1) max(eet1)];

% plotting the true and computed error fields
figure
subplot(1,2,1)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),eet1(ii),'filled');axis equal;caxis(cce);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');
subplot(1,2,2)
ii=find(p1(1,:)>=p1(3,:));scatter3(p1(1,ii),p1(2,ii),p1(3,ii),25*ones(length(ii),1),ee1(ii),'filled');axis equal;caxis(cce);
hold
ii=find(p1(1,b1)>=p1(3,b1));scatter3(p1(1,b1(ii)),p1(2,b1(ii)),p1(3,b1(ii)),'k.');