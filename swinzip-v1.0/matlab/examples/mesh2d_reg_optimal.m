% This script compresses and reconstructs a large dataset defined on a 2D 
% structured mesh. The datasets is composed of different subdomains where
% the optimal compression ratio is estimated before performing the
% compression.
% This examples shows how such datasets can be processed in parallel

% The datasets are obtained from a 2D RCCI trubulent combustion simulation

% MANY thanks to Jackie Chen and Ankit Bhagatwala from Sandia Labs for
% providing the data!

% Details on each dataset are found in SWinzip_ROOT/data/README.txt

clear all
fi=2;
ff=['rcci_Y_CO2_000' num2str(fi) '.dat'];

% aggregation levels in each direction
% The data chunk size increases with s1 and s2
% s1 and s2 have to be equal unless the code below is tweaked to prevent
% errors
% For best results the maximum value for s1 and s2 is 4.
s1=4;s2=4; 

f=fopen(['../../data/' ff]);

% n is the size of data in the x and y directions
n=[1600 1600];
% np is the way it is divided into subdomains in the x and y directions
np=[40/s1 40/s2];

if (fi>2)
    n=n*2;
    np=np*2;
end

% setting the options
options.vm=5;
options.mt='bernoulli';
options.j=-1;
options.seed=123;
options.part_type='axis';

% Setting a minimum limit to the compression ratio
Rmin=10;
A=fread(f,n,'float32');

% data skipped during plotting, if sk=1 then no data is skipped
sk=1;

% Setting data plot range
plotr=[min(min(A)) max(max(A))];

% Filling the coordinates and data for each subdomain into matrices
npt=prod(np);
nl=[n(1)/np(1) n(2)/np(2)];
B=zeros(nl(1),nl(2),npt);
xyl=zeros(2,2,prod(np));

for i=1:np(1)
    for j=1:np(2)
        ii=j+np(1)*(i-1);
        xyl(:,:,ii)=[i/np(1) (i-1)/np(1);j/np(2) (j-1)/np(2)];
    end
end

for j=1:np(1)
    for i=1:np(2)
        m=(j-1)*nl(1)+1;
        p=(i-1)*nl(2)+1;
        ii=j+np(1)*(i-1);
        B(:,:,ii)=A(m:m+nl(1)-1,p:p+nl(2)-1);
    end
end
dx=1/n(1);
dy=1/n(2);

R=zeros(1,npt);
nl2=prod(nl);
AA=1/prod(n);

% a is the vector of coefficients for the optimal compression estimation
a = [-2.8341 ;  0.3863 ;  0.0889 ;  0.1601 ;  0.0430 ; -0.4153];

% computing the compression ratio for all subdomains
for i=1:npt
   C=B(:,:,i);
   
   % normalize the data by the max-min
   Cn=C/(max(reshape(C,1,size(C,1)*size(C,2)))-min(reshape(C,1,size(C,1)*size(C,2))));
   dCngx=(Cn(2:end,2:end)+Cn(1:end-1,2:end)-Cn(1:end-1,1:end-1)-Cn(2:end,1:end-1))/dx/2;
   dCngy=(Cn(2:end,2:end)-Cn(1:end-1,2:end)-Cn(1:end-1,1:end-1)+Cn(2:end,1:end-1))/dy/2;
   
   % normalize the gradient by the mesh size
   dCng1=sqrt(dCngx.^2+dCngy.^2)*sqrt(AA);

   mdCng=mean(reshape(dCng1,1,size(dCng1,1)*size(dCng1,2)));                  % mean
   sdCng=std(reshape(dCng1,1,size(dCng1,1)*size(dCng1,2)));                   % standard deviation
   mmed=abs(mdCng - median(reshape(dCng1,1,size(dCng1,1)*size(dCng1,2))));    % mean minus the median
   mx=max(abs(reshape(dCng1,1,size(dCng1,1)*size(dCng1,2))));                 % maximum
      
    RR=exp([1 -log(mx) -log(sdCng) -log(mdCng) -log(mmed) -log(prod(n./np))]*a);    
    R(i)=RR;

    % limiting the compression ratio
    if (R(i)<=Rmin)
        R(i)=Rmin;
    end
end

% Plotting the original data
A=A(1:sk:end,1:sk:end);
[x,y]=ndgrid(0:1/(n(1)/sk-1):1,0:1/(n(2)/sk-1):1);
subplot(1,2,1)
pcolor(x,y,A');shading interp
caxis(plotr)
xlim([0 1])
ylim([0 1])
% colorbar
axis square
title('Original data field');

s=zeros(npt,1);
inc=s;
ew=s;
es=s;
NN=s;thr=s;normya=s;
xw=zeros(nl2,npt);
xs=zeros(nl2,npt);
p=zeros(2,nl2,npt);
d=xw';

x_I0=[];
activeSet0=[];
x_I=[];
activeSet=[];

cd ../src

% Since the mesh is the same in all subdomains, we compute the matrices
% only once prior to compression/decompression
disp('Building large matrices');
V = build_alpert_matrix(get_mesh(xyl(:,:,1),nl),options.vm,options.j,options.part_type);
P = constructSamplingMatrix(ceil(nl2/Rmin),nl2,'bernoulli');
disp('Building large matrices product');
T = P*V;
disp('done building large matrices');

% loop all over subdomains. This can be performed in parallel by replacing
% 'for' with 'parfor'
x_I=[];
activeSet=[];
tic
for i=1:npt
    i
   p(:,:,i)=get_mesh(xyl(:,:,i),nl);
   d(i,:)=reshape(B(:,:,i),1,nl2);
   [z,xw(:,i),xs(:,i),ew(i),es(i),inc(i)]=CompressReconstruct(p(:,:,i),d(i,:), R(i), options, x_I0, activeSet0, V, P, T);
   s(i)=length(d(i,:))/R(i);
end
toc
esg=sqrt(sum(sum((d-xs').^2))/prod(n))/(max(max(d))-min(min(d)));
ewg=sqrt(sum(sum((d-xw').^2))/prod(n))/(max(max(d))-min(min(d)));

cd ../examples

% Plotting the reconstructed data
subplot(1,2,2)
hold on
for i=1:npt
   x=reshape(p(1,:,i),nl(1),nl(2));
   y=reshape(p(2,:,i),nl(1),nl(2));
   x=x(1:sk:end,1:sk:end);
   y=y(1:sk:end,1:sk:end);
   dd=reshape(xs(:,i),nl(1),nl(2))';
   dd=dd(1:sk:end,1:sk:end);
   pcolor(x,y,dd)
end
caxis(plotr)
xlim([0 1])
ylim([0 1])
shading interp
% colorbar
axis square
title('Data reconstructed from compressed samples');

RR=prod(n)/sum(s);
fprintf('The normalized relative mean square error (NRMSE) \n of reconstruction from compressed samples: %f \n', esg);
fprintf('Average Compression ratio: %f \n', RR);
