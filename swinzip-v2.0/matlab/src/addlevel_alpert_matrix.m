function [V,Uo,Mo,Ujj] = addlevel_alpert_matrix(Ui,Mi,U,pos,cj,nJ,k2,Jmax)

%  Written by Maher Salloum

% addlevel_alpert_matrix_2d - adds detail levels to an existing Alpert
% matrix
%
%   Ui and Mi: moments matrices at level cj

%   U: is the existing matrix compute at level cj 

%   pos: is a 2D matrix, pos(:,i) is the ith point.

%   cj: is the current detail level

%   nJ: is the number of levels to add

%   k2: is the full wavelets order in all directions
%       This should be the same as the one used to build U

%   Jmax: is the maximum number of levels
%         This should be the same as the one used to build U

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of points
n = size(pos,2);
d = size(pos,1);

J = Jmax;
P=2^(J-1);

if (cj<=0)
    disp('This function cannot be used to start building detail levels from scratch, use build_alpert_matrix_2d instead.');
    return
end

if (cj>=J)
    disp('Cannot add anymore detail levels');
    return
end

if ((cj+nJ)>=J)
    disp('Maximum detail levels reached');
    nJ=J-cj;
end

Ujj = speye(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=cj:cj+nJ-1   % for each scale
    
    % at this scale, we have nj = P/2^(j-1) groups
    nj = P/2^j ; % n/(2^j*k);
    mj = nj*2*k2;     % total length of the blocks
    
    % update each sub matrix
    for i = 1:nj
        M = zeros(2*k2, 2*k2);
%         M = spalloc(2*k2,2*k2,ceil(4*k2^2*0.1));        
        M(1:k2,:)           = Ui{2*i-1}(1:k2,:)  *   Mi{2*i-1};         % Ui^U is just k first row
        M((k2+1):2*k2,:)    = Ui{2*i}(1:k2,:)    *   Mi{2*i};
        MMi{i} = M;
        [UUi{i},R] = qr(MMi{i});    UUi{i} = transpose(UUi{i});
        
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
    Ujj=Uj*Ujj;
end

V = U';
Uo=Ui;
Mo=Mi;

% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);