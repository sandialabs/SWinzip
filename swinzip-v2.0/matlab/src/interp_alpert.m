function [f0,ww,ew]=interp_alpert(p,f,p0,th0,vm,mxmn)

N=size(p,1);
pa=[p;p0];
N0=size(p0,1);

ar=max(pa,[],1)-min(pa,[],1);
ar=max(ar)/min(ar);

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

% iterating over wavelet coefficients and adjusting signs
[aw,jj]=sort(abs(wa),'descend');
for i=1:M
    % flip the sign of a coefficient
    wa(jj(i))=-wa(jj(i));
    fw = perform_alpert_transform(wa,Uja,parta,k2a,-1);
    ew1=sqrt(sum((f-fw(1:N)).^2));
    
    % if the error does not decrease, flip the sign back
    if (ew1<ew)
        ew=ew1;
    else
        wa(jj(i))=-wa(jj(i));
    end
end

f0=fw(N+1:end);
ww=wa(Na-M:end)/sqrt(Na)*sqrt(N);
