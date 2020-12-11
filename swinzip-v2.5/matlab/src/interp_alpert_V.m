function [e,e0]=interp_alpert_V(pre,pr,fre,fr,vm)

V = build_alpert_matrix(pr,vm,-1);
Ve = build_alpert_matrix(pre,vm,-1);

N=size(pr,2);
Ne=size(pre,2);

pr=bsxfun(@minus,pr,min(pr,[],2));
pr=bsxfun(@times,pr,1./max(pr,[],2));
pre=bsxfun(@minus,pre,min(pre,[],2));
pre=bsxfun(@times,pre,1./max(pre,[],2));

% computing signs of wavelet non-vanishing moments
s=sign(sum(bsxfun(@times,V,prod(pr.^vm)')));
se=sign(sum(bsxfun(@times,Ve,prod(pre.^vm)')));

if (Ne>N)
    ss=se(Ne-N+1:Ne).*s;
    VV=bsxfun(@times,V,ss);    

    we=Ve'*fre/sqrt(Ne);

    fre_pr=VV*we(Ne-N+1:Ne)*sqrt(N);

    e=fre_pr - fr;

    [wes,jj]=sort(abs(we),'descend');
    wwe=V'*fre_pr/sqrt(N);
    [wwes,jj]=sort(abs(wwe),'descend');

    e0=norm(wes(1:N)-wwes)/wes(1)/sqrt(N);
else
    ss=se.*s(N-Ne+1:N);
    VV=bsxfun(@times,V(:,N-Ne+1:N),ss);    

    we=Ve'*fre/sqrt(Ne);

    fre_pr=VV*we*sqrt(N);

    e=fre_pr - fr;

    [wes,jj]=sort(abs(we),'descend');
    wwe=V'*fr/sqrt(N);
    [wwes,jj]=sort(abs(wwe),'descend');

    e0=norm(wes-wwes(1:Ne))/wes(1)/sqrt(Ne);    
end
    
