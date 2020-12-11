function [part,B,G] = dichotomic_grouping(pos, ptmax, part_type)

% dichotomic_grouping - regroup points into 2^s (s is unknown) bins of nearly equal size.
% it is used when performing a wavelet transform on an unstructured set of
% points
%
% pos: the points coordinates. It is a matrix of size d x N, where d is the
%      number of dimensions and N the number of points.

% ptmax: is the full wavelets order in all dimensions

% part_type: is the partitioning type. It can be either:
%             '2axis': for a partioning based on the corrdinates directions
%                      in a 2D domain
%             '3axis': for a partioning based on the corrdinates directions
%                      in a 3D domain                       
%             'kmeans': for a partioning based on the k-means clustering

% This code has been implemented from the library of Gabriel Peyre found
% at: https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
% It was slightly modified to fit the needs of the SWinzip library

n = size(pos, 2);
d = size(pos, 1);   % number of dimensions

depthmax = Inf;
nb_iter = 10;

part = {1:n};
cpt = 0;
type=part_type;

if length(type)>=5 && strcmp(type(2:end), 'axis')
    % retrieve the number of axis
    type = str2num(type(1));
    type = max(min(type,d), 0);     % should be in [0,d]
else
    % otherwise use kmeans grouping
    type = 0;
end

while cpt<depthmax  % to avoid infinite loop ...
    cpt = cpt+1;
    
    part_prev = part;
    
    % perform binary split
    part1 = part;
    part = {};
    for k=1:length(part1)
        P = part1{k};
        pp = floor(length(P)/2);    % half # points
        % compute sub-partition
        if type>0
            % axis-split clustering
            s = mod(cpt-1,type)+1;
            [tmp,I] = sort( pos(s,P) );
            ppart1 = I(1:pp);
            ppart2 = I(pp+1:end);
        else
            % k-mean clustering
            [ppart1,ppart2] = dist_part( pos(:,P), nb_iter );
        end
        % assign new partition
        part{2*k-1} = P(ppart1);
        part{2*k-0} = P(ppart2);
    end
    
    for i=1:length(part)
        if length(part{i})<ptmax
            part = part_prev;
            B = zeros(n,1);
            for i=1:length(part)
                B(part{i}) = i;
                G(i) = length(part{i});  
            end
            return;
        end
    end
end

B = zeros(n,1);
for i=1:length(part)
    B(part{i}) = i;
    G(i) = length(part{i});    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J] = dist_part(X,nb_iter)

% dist_part - partition a set of points into too components.
%
%   [I,J] = dist_part(X);
%
%   Copyright (c) 2004 Gabriel Peyr�

if nargin<2
    nb_iter = 5;
end

n = size(X,2);
n1 = floor(n/2);
d = size(X,1);

% try to reduce a bit the number of iterations
nb_iter = min(nb_iter,ceil(n/2-1));

% random initialization
seed_num = floor(rand(2,1)*n)+1;
if seed_num(2)==seed_num(1)
    seed_num(2) = seed_num(1) + 1;
    if seed_num(2)>n
        seed_num(2) = 1;
    end
end
seeds = X(:,seed_num);

% compute region of influence
D = sqrt( compute_distance_to_points(X,seeds) );
[tmp,A] = sort( D(1,:)-D(2,:) );
part{1}= A(1:n1); part{2} = A(n1+1:end);

[tmp,B] = min(D);

for i=1:nb_iter
    % compute region center
    for k=1:2
        % geometric barycenter
        if ~isempty(part{k})
            seeds(:,k) = sum( X(:,part{k}), 2)/length( part{k} );
        else
            warning('Empty cluster created.');
        end
    end
    % compute region of influence
    D = sqrt( compute_distance_to_points(X,seeds) );
    [tmp,A] = sort( D(1,:)-D(2,:) );
    part{1}= A(1:n1); part{2} = A(n1+1:end);
end

I = part{1};
J = part{2};


function D = compute_distance_to_points(X,seeds)

% compute_distance_to_points - compute euclidean distance to a set of points.
%
%   D = compute_distance_to_points(X,seeds)
%
%   'X' is a [d,n] matrix, X(:,i) is the ith point living in R^d.
%   'seeds' is a [d,k] matrix.
%   D(i,j) = |X(:,j)-seeds(:,i)|^2.
%
%   Copyright (c) 2004 Gabriel Peyre

%   This code was obtained from:
%   https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
%   published under a BSD license

nbCluster = size(seeds,2);
n = size(X,2);
D = zeros(nbCluster,n);
d = size(X,1);

for k=1:nbCluster
    % distance to seed
    D(k,:) = sum( (X - repmat(seeds(:,k),1,n)).^2 );
end