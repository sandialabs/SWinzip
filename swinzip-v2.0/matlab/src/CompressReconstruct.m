function [y,xw,xs,ew,es,inc, x_I, activeSet,ts, wc, D]=CompressReconstruct(p, x, R, options, x_I0, activeSet0, V, P, T)

% This function is a pipeline takes data as input, compresses it, performs reconstruction
% then compute error metrics.

% inputs:
% p: mesh coordinates matrix (d x N)

% x: data vector (1 x N)

% R: compression ratio (double scalar)

% options: a struct with the folloeing entities:
%          .vm: wavelet order aka number of vanishing moments (integer)
%          .mt: type of sampling matrix (string)
%          .j: number of detail levels requested (integer). If j=-1
%              then the wavelet matrix is computed at full details
%          .seed: seed for sampling matrix computation (integer). If this number is
%                 less or equal than zero then a random number is generated
%          .part_type: the partition type
%                      'axis': for a partioning based on the corrdinates directions
%                      'kmeans': for a partioning based on the k-means clustering

% R: compression ratio (double scalar)

% x_I0: Initial solution vector for StOMP (vecor of double)

% activeSet0: Initial locations of active coefficients in the solution
%             vector for StOMP (vector of integer)

% V, P and T are wavelet, sampling and measurement matrices (T=P*V).
% They are optional parameters used if they are precomputed which is
% useful for different datasets based on the same mesh

% outputs:
% y: wavelet coefficients vector recontructed by StOMP (N x 1)
% xw: solution vector recontructed by inverse wavelet transform (N x 1)
% xs: solution vector recontructed by StOMP (N x 1)
% ew: error obtained by inverse wavelet transform (double scalar)
% es: error obtained by StOMP (double scalar)
% inc: incoherence between sampling and wavelet matrices (double scalar)
% x_I0: solution vector for StOMP (vecor of double)
% activeSet0: locations of active coefficients in the solution
%             vector for StOMP (vector of integer)
% ts: time taken by StOMP in seconds (double scalar)
% wc: wavelet coefficients solved by StOMP (N x 1)
% D: the Kuan-Liu Ma error metrics

M=size(p);
N=M(2);
d=M(1);

vm = get_option(options, 'vm', 5);
j = get_option(options, 'j', -1);
mt = get_option(options, 'mt', 'bernoulli');
seed = get_option(options, 'seed', 123);
part_type = get_option(options, 'part_type', 'axis');

if (nargin<7)
    disp('Computing wavelet matrix');
    tic
    V = build_alpert_matrix(p,vm,j,part_type);
    fprintf('wavelet matrix construction took %f\n', toc);
end

% Do foreward wavelet transform 
w=V'*x';
% Sort the wavelet vector
[ws,i]=sort(abs(w));

% NR is the number of remaining values after compression
NR=floor((1-1/R)*N);

if (nargin<7)
    % Start random number generator
    if (seed>0)
        rng(seed);
    end
    % Construct sampling matrix 
    P = constructSamplingMatrix( ceil(N/R), N, mt);

    disp('Computing sampling-wavelet matrix product');
    
    % Do the Compression
    tic
    y=P*x';
    % fprintf('compression took %f\n', toc);
    P=P*V;    
else
    P=P(1:ceil(N/R),:);
    
    % Do the Compression
    tic
    y=P*x';
    clear P
    % fprintf('compression took %f\n', toc);
    
    P=T(1:ceil(N/R),:);
end

% Compute incoherence 
inc=sqrt(N)*max(max(P));

% disp('Running StOMP');
tic
[sol, numIters, I, x_I, activeSet] = SolveStOMP(P, y, 50, 1, 1e-8, x_I0, activeSet0);
ts=toc;

% Compute the inverse wavelet tranform and the associated error
xw=V(:,i(NR:end))*w(i(NR:end));
ew=sqrt(sum((x'-V(:,i(NR:end))*w(i(NR:end))).^2)/N)/(max(x)-min(x));

% Compute the inverse wavelet transform based on the StOMP solution and the associated error
xs=V*sol;


es=sqrt(sum((x'-V*sol).^2)/N)/(max(x)-min(x));
wc = sol;

Vb = build_alpert_matrix(p,vm,1,part_type);
if (j==-1)
    if length(vm)==1
        vm = repmat(vm, 1, d);
    end
    [part,B,G] = dichotomic_grouping(p, prod(vm), [num2str(d) part_type]);
    J = log2(length(G))+1;
else
    J=j;
end
D=KLMmetrics(w,wc,J,Vb'*x',Vb'*xs,N/length(find(abs(w)>ws(end)/2^(J+4))));
