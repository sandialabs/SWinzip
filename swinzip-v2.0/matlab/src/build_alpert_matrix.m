function [V,J,Ui,Mi] = build_alpert_matrix(pos,k,J,part_type)

% build_alpert_matrix - builds an Alpert basis adapted to the provided mesh.
%   Extension of build_alpert_matrix_2d to build 3D Alpert basis

% pos: the points coordinates. It is a matrix of size d x N, where d is the
%       number of dimensions and N the number of points.

% k: is the wavelet order in each direction (a 1 x d vector). If only one number is provided
%    then the order will be the same in all directions

% J: is the number of levels required. If J=-1 then the full level details
%    are considered

% part_type: is the partitioning type. It can be either:
%             'axis': for a partioning based on the corrdinates directions
%             'kmeans': for a partioning based on the k-means clustering

% motherFunction: a choice of motherFunction for our tree wavelets

% This code has been implemented from the library of Gabriel Peyre found
% at: https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
% It was substantially modified and improved to fit the needs of the SWinzip library


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n = size(pos,2);
d = size(pos,1);    % number of dimensions

if length(k)==1
    % use same order for X and Y, etc
    k = repmat(k, 1, d);
end

if n==0
    V = [];
    return;
end

k2 = prod(k);

if n<k2
    % special case, not enough data
    k = repmat(floor(sqrt(n)),1,d);
    k2=prod(k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find a regroupement
if (nargin<4)
    part_type = 'axis';
end

part_type = [num2str(d) part_type];

[part,B,G] = dichotomic_grouping(pos,k2,part_type);

P = length(G);          % number of groups

% we have got nbr.packets = 2^(J-1)
if (J<=0)
    J = log2(P)+1;
end

si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

KK=make_table(k)-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of moments matrices
for i = 1:P
    
    seli = part{i};    
    
    % build moment matrix       
    kk = length(seli);                   % nbr of point in this bin, ie. 2k^2 in regular case
    M = zeros( kk, 2*k2 );
    for ii=1:kk    
        for j=1:k2
           M(ii,j)=1;
           for dd=1:d
              M(ii,j)=M(ii,j)*pos(dd,seli(ii))^KK(j,d-dd+1);
           end
        end
    end
        
    Mi{i} = M;    
    % orthogonalize
    [Ui{i},R] = qr(Mi{i});
    Ui{i} = transpose(Ui{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization of first matrix
U = spalloc(n,n,ceil(n^2*0.01));
for i = 1:P
    selj = part{i};                   % selected column
    % to keep : upper part is of size n-P*k^2
    offs = si(i)-(i-1)*k2;   % offset on row
    long = G(i)-k2;          % length on row
    U( offs + (0:long-1), selj )       = Ui{i}((k2+1):end, :);   % we keep the G(i)-k^2 last
    % to retransform : lower part is of size P*k2
    offs = n - P*k2+1 + (i-1)*k2;
    if offs>0
        U( offs+(0:k2-1), selj )   = Ui{i}(1:k2, :);   % the first k^2 doesn't have vanishing moments, keep them   
    else
        warning('Pbm empty bin.');
        % bug ...
    end
end

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
    
    % lower part of the multiplicative matrix
    Uj = speye(n);
    Uj( (end-mj+1):end, (end-mj+1):end ) = 0;
    for i = 1:nj  
        Uj( (n-mj)+k2*(i-1)+(1:k2), (n-mj)+2*k2*(i-1)+(1:2*k2) ) = mlow( Ui{i} );
        Uj( (n-mj)+mj/2+k2*(i-1)+(1:k2), (n-mj)+2*k2*(i-1)+(1:2*k2) ) = mup( Ui{i} );
    end
    
    % update vectors
    U = Uj*U;
end

V = U';


% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);
