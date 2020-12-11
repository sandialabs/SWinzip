function D=KLMmetrics(w,wc,J,w1,wc1,Rmin)

% A function to compute metrics from a paper by the group of Kuan-Liu Ma
% C. Wang and K-L. Ma, Statistical Approach to Volume Data Quality
% Assessment, IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 14, NO. 3, MAY/JUNE 2008

% w: The "true" wavelet representation
% wc: The approximate wavelet representation such as the one
%     obtained by StOMP
% J: The total number of levels in the wavelet transform.
%    It can be obtained when building the wavelets
% w1: The "true" level 1 wavelet representation of the data
% wc1: The level 1 wavelet representation of the data
%      recovered by StOMP
% Rmin: The minimum compression ratio that decides the minimum threshold
%        parameter to compute D1 abd D2

N=length(w);
[ws,i]=sort(abs(w));

% Rmin is a minimum compression ratio that sets
% a threshold below which the wavelets coefficients in
% the "true" wavelet representation w are considered small
% and can be set to zero
% Rmin=N/length(find(abs(w)>ws(end)/2^(J+4)));

% L is the number of coeffs considered small
L=N-floor(N/Rmin);

% D1
wwm=w;

% Setting the small coefficients to zero
wwm(i(1:L))=0;
wwm=setdiff(wwm,0);

% The pdf of the "true" wavelet coefficients is
% ugly in the sense that most of the coeffs are centered
% around zero while the large ones are VERY few and form
% a set of outliers. These cannot be neglected because they
% contain most of the info in the data
% Here we separate the large coeffs on the positive and
% negative sides
iin=find(wwm<prctile(wwm,5));
iip=find(wwm>prctile(wwm,95));
pxp=sort(wwm(iip));
pxn=sort(wwm(iin));

% px is the support of the pdf of the "true" data 
px=min(pxn):(max(pxp)-min(pxn))/100000:max(pxp);

% ff is the pdf of the "true" data
% ksdensity is part of the statistics toolbox
ff=ksdensity(wwm,px);

% ffc is the pdf of the approximate data
ii=find(wc~=0);
if isempty(ii)
    ffc=ksdensity([0],px);
else
    ffc=ksdensity(wc(ii),px);
end

% Make sure that only positive values exist in the pdf
ii=find(ff>0);
ff=ff(ii);ffc=ffc(ii);px=px(ii);
ii=find(ffc>0);
ff=ff(ii);ffc=ffc(ii);px=px(ii);

% Compute D1
ff=ff.*log(ff./ffc);
ffc=diff(px);
if (length(ffc)==0)
    D1 = 0;
else
    D1=log(1+sum(ff(1:end-1).*ffc));
end

% D2
wmax=ws(end);
D2=log(sqrt((sum((w(i(L:end))-wc(i(L:end))).^2)+sum(wc(i(1:L-1)).^2))/N/wmax^2)+1);

% D3
bi=wc1;
bj=w1;
mi=mean(bi);mj=mean(bj);
si=std(bi);sj=std(bj);
sij=norm(cov(bi,bj))/2;
D3 = sqrt(abs(1-(4*sij*mi*mj/(si^2+sj^2)/(mi^2+mj^2)+1)/2));

D=[D1 D2 D3];