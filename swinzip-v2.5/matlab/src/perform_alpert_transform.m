function w = perform_alpert_transform(v,Uj,part,k2,dir)

% perform_alpert_transform - performs forward and inverse Alpert wavelet 
%                            transform based on a provided mesh.

% v: data represented on an unsructured mesh, a vector of length N

% Uj: structure containing moment matrices pre-computed  based on the
%     unsructured mesh at which v is represented.

% part: structure containing the unsructured mesh hirearchy.

% k2: the full wavelet order. It should be the same as the one used to
%     pre-compute Uj and part

% dir: an index equal to:
%      1 for forward wavelet transform
%      -1 for inverse wavelet transform

% This code has been implemented from the library of Gabriel Peyre found
% at: https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_alpert
% It was substantially modified and improved to fit the needs of the SWinzip library

v = v(:);

n=length(v);
J=length(Uj);
P=2^(J-1);

for i=1:P
    G(i) = length(part{i});    
end
si = [0, cumsum(G)]+1;    % si(i) is the index of the 1st point of ith group

if dir==1
    j_list= 1:J;
else    % transpose reverse order
    j_list= J:-1:1;
end

w = v;
for j=j_list
    
    Ui = Uj{j}; 
    if j==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % special treatment for 1st scale 
        r = w;
        for i = 1:P
            selj = part{i};                   % selected column
            % to keep : upper part is of size n-P*k^2
            offs = si(i)-(i-1)*k2;   % offset on row
            long = G(i)-k2;          % length on row
            % we keep the G(i)-k^2 last
            seli = offs + (0:long-1);
            U = Ui{i}((k2+1):end, :);
            if dir==1
                r(seli) = U * w(selj);
            else
                r(selj) = U' * w(seli);    
            end
            
            % to retransform : lower part is of size P*k2
            offs = n - P*k2+1 + (i-1)*k2;
            if offs>0
                seli = offs+(0:k2-1);
                % the first k^2 doesn't have vanishing moments, keep them   
                U = Ui{i}(1:k2, :);
                if dir==1
                    r(seli) = U * w(selj);
                else
                    r(selj) = r(selj)  + U' * w(seli);    
                end
            else    % bug ...
                % warning('Problem empty bin.');
            end
        end
        w = r;
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % at this scale, we have nj = P/2^(j-1) groups
        nj = P/2^(j-1) ; % n/(2^j*k);
        mj = nj*2*k2;     % total length of the blocks
        
        % lower part of the multiplicative matrix
        r = w;
        offs = n-mj;
        for i = 1:nj
            selj = offs + 2*k2*(i-1)+(1:2*k2);
            seli = offs + k2*(i-1)+(1:k2);
            if dir==1
                r(seli) = mlow( Ui{i} ) * w(selj);
            else
                r(selj) = mlow( Ui{i} )' * w(seli);    
            end
            
            seli = offs + mj/2+k2*(i-1)+(1:k2);
            if dir==1
                r(seli) = mup( Ui{i} ) * w(selj);
            else
                r(selj) = r(selj) + mup( Ui{i} )' * w(seli);    
            end
        end
        w = r;
    end
    
end

% extract uper part of the matrix
function MU = mup(M)
MU = M(1:end/2, :);

% extract lower part of the matrix
function ML = mlow(M)
ML = M((end/2+1):end, :);
