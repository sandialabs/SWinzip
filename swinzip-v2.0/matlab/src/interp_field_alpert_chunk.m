function [f0,w]=interp_field_alpert_chunk(p,f,p0,th0,vm,mxmn)

% This function uses wavletes to interpolate a field f represented on a
% mesh p, onto another mesh p0.

% p: is the mesh matrix of size N x d where d is the number of dimensions
%    and N is the number of points
% f: is the function represented on the mesh p and to be interpolated on p0
% p0: mesh where f will be interpolated
% th0: wavelet threshold multiplier (typical value is 0.0075)
% vm: wavelet order (typical value is 5)
% mxmn: reference maximum of the global f field minus its minimum

% f0: interpolated field on p0
% w: wavelet transform of f0

N=size(p,1);
pa=[p;p0];
N0=size(p0,1);

mxmnf=prctile(f,99)-prctile(f,1);

[Uj,part,k2] = compute_Uj(p',vm);
w = perform_alpert_transform(f,Uj,part,k2,1);
aw=abs(w)/mxmn;
th=th0/sqrt(mxmn/mxmnf)/mean(aw)*mean(abs(aw-mean(aw)))/0.6745*sqrt(2*log(N));
ii=find(aw>th);
M=length(ii);

[Uja,parta,k2a] = compute_Uj(pa',vm);

wa=zeros(size(pa,1),1);
Na=length(wa);
wa(ii+N0)=w(ii)/sqrt(N)*sqrt(N+N0);
fw=perform_alpert_transform(wa,Uja,parta,k2a,-1);
ew=sqrt(sum((f-fw(1:N)).^2));

[aw,jj]=sort(abs(wa),'descend');

for i=1:M
    wa(jj(i))=-wa(jj(i));
    fw = perform_alpert_transform(wa,Uja,parta,k2a,-1);
    ew1=sqrt(sum((f-fw(1:N)).^2));
    if (ew1<ew)
        ew=ew1;
    else
        wa(jj(i))=-wa(jj(i));
    end
end

f0=fw(N+1:end);

w=wa(Na-M:end);
