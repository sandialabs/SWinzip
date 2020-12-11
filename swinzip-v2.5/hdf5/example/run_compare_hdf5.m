clear all

d=importdata('out.txt');
N=length(d);
d=reshape(d,N/2,2);

norm(d(:,1)-d(:,2))/sqrt(N)/(max(d(:,1))-min(d(:,1)))

subplot(1,2,1)
slice(reshape(d(:,1),200,180,240),[1 0 0],[0 1 0],[0 0 1]);shading interp;axis equal

subplot(1,2,2)
slice(reshape(d(:,2),200,180,240),[1 0 0],[0 1 0],[0 0 1]);shading interp;axis equal;caxis([min(d(:,1)) max(d(:,1))])