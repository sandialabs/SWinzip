function y=CompressBernoulli(x, M, seed)

% This code performs a matrix-free Bernoulli sampling on a vector x of length N.
% The result is a column vector of length M < N

x=x(:);
N=length(x);
y=zeros(M, 1);

% Default seed
if (nargin<3)
    seed=123;
end

rng(seed);
for j=1:M
    y(j)=(mod(randperm(N),2)*2-1)*x;
end

y=y/sqrt(N);