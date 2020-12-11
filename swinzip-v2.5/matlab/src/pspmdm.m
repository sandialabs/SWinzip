function T = pspmdm(P,V)

% Computes the product of a dense matrix P and a sparse square matrix V in
% parallel, it is sometimes faster than the direct product especially when
% the number of columns in P (and T) is large.

T=zeros(size(P));
parfor i=1:size(P,1)
    T(i,:)=P(i,:)*V;
end

end

