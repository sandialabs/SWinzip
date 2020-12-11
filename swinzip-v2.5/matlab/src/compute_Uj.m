function [Uj,part,k2] = compute_Uj(pos,k)

% compute_Uj - computes moment matrices and mesh hierarchy.

% pos: the points coordinates. It is a matrix of size d x N, where d is the
%       number of dimensions and N the number of points.

% k: is the wavelet order in each direction (a 1 x d vector). If only one number is provided
%    then the order will be the same in all directions

% This code has been implemented from the library of Gabriel Peyre found
% at: https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
% It was substantially modified and improved to fit the needs of the SWinzip library

n = size(pos,2);    % number of points
d = size(pos,1);    % number of dimensions

if length(k)==1
    % use same order for X and Y, etc
    k = repmat(k, d, 1);
end

if n==0
    w = [];
    return;
end

if n<prod(k)
    % special case, not enough data.
    % We want prod(k)<n i.e. k^d=n
    k = floor( n^(1/d) );
    k = repmat(k, d, 1);
end
k2 = prod(k);             % equivalent to k^2 in Alpert 2D.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
options.ptmax = k2;
options.part_type = [num2str(d) 'axis'];
% options.part_type = 'kmeans';

part_type = 'axis';
part_type = [num2str(d) part_type];

[part,B,G] = dichotomic_grouping(pos,k2,part_type);


P = length(part);          % number of groups
for i=1:P
    G(i) = length(part{i});    
end
si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

% we have got nbr.packets = 2^(J-1)
J = log2(P)+1;

KK=make_table(k)-1;

% J=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Sparse Matrix construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices (j=1)
for i = 1:P
    [Mi{i}, Ui{i}]=compute_moment_matrix(pos,k,KK,d,part{i},k2);
end

Uj{1} = Ui;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=2:J   % for each scale
    % at this scale, we have nj = P/2^(j-1) groups
    nj = P/2^(j-1) ; % n/(2^j*k);
    mj = nj*2*k2;     % total length of the blocks
    
    % update each sub matrix
    for i = 1:nj
        M = zeros(2*k2, 2*k2);
        M(1:k2,:)           = Ui{2*i-1}(1:k2,:)  *   Mi{2*i-1};         % Ui^U is just k first row
        M((k2+1):2*k2,:)    = Ui{2*i}(1:k2,:)    *   Mi{2*i};
        MMi{i} = M;
        [UUi{i},R] = qr(MMi{i});    
        UUi{i} = transpose(UUi{i});
    end
    Mi = MMi;
    Ui = UUi;    
    Uj{j} = Ui;
end
