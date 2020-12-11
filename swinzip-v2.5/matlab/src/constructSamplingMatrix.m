function [ PsiM ] = constructSamplingMatrix( nrows, ncols, mtype) 
%constructSamplingMatrix Construct a sampling matrix of the specified type
%and size
%   [ PsiM ] = constructSamplingMatrix( nrows, ncols, mtype) creates a
%   sampling matrix, nrows X ncols in size. It is of type mtype.
%
% INPUTS
%   nrows, ncols : The sampling matrix PsiM is nrows X ncols in size. ncols
%   is the size of the signal it will sample. nrows is the number of
%   measurements/samples produced.
%
%   mtype :  type of sampling matrix. 'random'/'Bernoulli'/'Hadamard'/
%
% OUTPUTS :
%   PsiM : nrows x ncols sampling matrix.
%
% This code is written by Jaideep Ray and Jina Lee and published at:
% http://www.sandia.gov/~jairay/kfcs/kfcs.html#Releases
% under a BSD license

% -------------------------------------------------------------------------

switch lower(mtype)
    case 'gaussian'
        fprintf('Generating a %d X %d Gaussian sampling matrix\n', nrows, ncols);
        PsiM = feval(@randomMatrix, nrows, ncols) ;
    case 'bernoulli'
        fprintf('Generating a %d X %d Bernoulli sampling matrix\n', nrows, ncols);
        PsiM = feval(@bernoulliMatrix, nrows, ncols) ;
    case 'hadamard'
        fprintf('Generating a %d X %d Hadamard sampling matrix\n', nrows, ncols);
        PsiM = feval(@hadamardMatrix, nrows, ncols) ;
    case 'sbhadamard'
        fprintf('Generating a %d X %d Scrambled Block Hadamard sampling matrix\n', nrows, ncols);
        PsiM = feval(@sbhadamardMatrix, nrows, ncols) ;
    case 'fourier'
        fprintf('Generating a %d X %d Fourier sampling matrix\n', nrows, ncols);
        PsiM = feval(@fourierMatrix, nrows, ncols) ;
    case 'toeplitz'
        fprintf('Generating a %d X %d Toeplitz sampling matrix\n', nrows, ncols);
        PsiM = feval(@toeplitzMatrix, nrows, ncols) ;
    case 'circulant'
        fprintf('Generating a %d X %d Circulant sampling matrix\n', nrows, ncols);
        PsiM = feval(@circulantMatrix, nrows, ncols) ;
    case 'noiselet'
        fprintf('Generating a %d X %d Noiselet sampling matrix\n', nrows, ncols);
        PsiM = feval(@noiseletMatrix, nrows, ncols) ;
    otherwise
        fprintf('Did not understand sampling matrix type %s\n', mtype) ;
        error('Abort!') ;
end
end

% =========================================================================

function [PsiM] = randomMatrix(nrows, ncols)
%randomMatrix Construct a random matrix to be used in compressive sensing. 
%  PsiM = randomMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. Each row in this matrix is normalized to
%  be a unit vector
%
% Jaideep Ray, 11/9/2012
% -------------------------------------------------------------------------

% ---- Random matrix, independent draws from N(0, 1)
tmp = randn(nrows, ncols) ; 

% ---- Entries should really be from N(mean = 0, var = npixels), where
% npixels is the size of the image = ncols. 
PsiM = 1 / sqrt(ncols) .* tmp ;
end

% =========================================================================

function [PsiM] = bernoulliMatrix(nrows, ncols)
%bernoulliMatrix Construct a bernoulli matrix to be used in compressive sensing. 
%  PsiM = bernoulliMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. Each row in this matrix is normalized to
%  be a unit vector
%
%  Bernoulli matrix has independent draws of -1/sqrt(npixels),
% 1/sqrt(npixels)
%%
% Jaideep Ray, 11/9/2012
% -------------------------------------------------------------------------

% create a vector containing 1/sqrt(npixels) and -1/sqrt(npixels)
rowelems = ones(ncols,1) ; rowelems( (floor(ncols/2)+1) : ncols) = -1 ;
rowelems = 1/sqrt(ncols) * rowelems ;

% create a nrows x ncols matrix containing indices of rowelems
inds = randi( length(rowelems), nrows, ncols ) ;

PsiM = rowelems( inds ) ; 

end

% =========================================================================

function [PsiM] = hadamardMatrix(nrows, ncols)
%hadamardMatrix Construct a partial hadamard matrix to be used in compressive sensing. 
%  PsiM = hadamardMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. nrows are randomly selected from a MxM
%  hadmard matrix.
%
% Jina Lee, 11/15/2012
% -------------------------------------------------------------------------

% create a ncols x ncols Hadamard matrix
hadamardM = hadamard(ncols);

% scramble up the columns
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% cinds = randperm(ncols, ncols);
cinds = randperm(ncols);   % JR_OLD_MATLAB
hadamardM = hadamardM(:, cinds);

% select nrows randomly from the M x M Hadamard matrix, creating a nrows x
% ncols matrix
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% inds = randperm( ncols, nrows );
inds = randperm( ncols ); inds = inds(1 : nrows ) ; % JR_OLD_MATLAB
PsiM = hadamardM(inds, :) ; 

% Normalize each elem by 1/sqrt(ncols)
PsiM = PsiM .* (1/sqrt(ncols));

% arbitrarily select half of the first 
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% rowinds = randperm( nrows, 1/2*nrows );
rowinds = randperm( nrows ); rowinds = rowinds ( 1 : (1/2*nrows) ) ; 
PsiM(rowinds,:) = -1*PsiM(rowinds,:);
end

% =========================================================================

function [PsiM] = sbhadamardMatrix(nrows, ncols)
%sbhadamardMatrix Construct a scrambled block hadamard matrix to be used in
%compressive sensing. 
%  PsiM = sbhadamardMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. nrows are randomly selected from a MxM
%  hadmard matrix.
%  
%  See (Gan,Do&Tran, 2008) for details.
%
% Jina Lee, 12/12/2012
% -------------------------------------------------------------------------

% create a BxB Hadamard matrix. (B = 32, 64, 128,... the smaller the better,
% but harder to achieve good result)
B = 256;
Wb = hadamard(B);

% create a block diagonal matrix W from Wb
numWb = floor(ncols/B)+mod(ncols,B);
W = [];
for i=1:numWb
    W = blkdiag(Wb, W);
end

% scramble the columns of W
cinds = randperm(size(W,2));
WPn = W(:, cinds);

% select nrows randomly from W creating a nrows x ncols matrix
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% rinds = randperm( ncols, nrows );
rinds = randperm( ncols); rinds = rinds(1:nrows) ; % JR_OLD_MATLAB
PsiM = WPn(rinds, :) ; 

% Normalize each elem by 1/sqrt(ncols)
% JR_FIX. 2/6/13. The elements of PsiM are not all +/- 1 (there are zeros)
% so do explicit normalization.
% PsiM = PsiM .* (1/sqrt(ncols));
for i = 1 : size(PsiM, 1)                                     % JR_NORM_FIX
   denom = norm( PsiM(i, :) ) ;
   PsiM(i, :) = PsiM(i, :) / denom ;
end

end

% =========================================================================

function [PsiM] = fourierMatrix(nrows, ncols)
%fourierMatrix Construct a partial fourier matrix to be used in compressive sensing. 
%  PsiM = fourierMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. nrows are randomly selected from a MxM
%  fourier matrix.
%
% Jina Lee, 11/15/2012
% -------------------------------------------------------------------------

% create a ncols x ncols Fourier matrix
fourierM = fft(eye(ncols));

% create a nrows x 1 matrix containing indices of rows to be selected from fourierM
%
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% inds = randperm( ncols, nrows );
inds = randperm( ncols ) ; inds = inds( 1: nrows ) ; % JR_OLD_MATLAB

% select nrows randomly from the M x M discrete Fourier matrix, creating an nrows x
% ncols matrix
PsiM = fourierM(inds, :) ; 

% Normalize each elem by 1/sqrt(ncols)
PsiM = PsiM .* (1/sqrt(ncols));
end

% =========================================================================

function [PsiM] = toeplitzMatrix(nrows, ncols)
%toeplitzMatrix Construct a partial Toeplitz matrix to be used in compressive sensing. 
%  PsiM = toeplitzMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. 

% Note on Toeplitz matrix (source: Wikipedia): a Toeplitz matrix or
% diagonal-constant matrix, is a matrix in which each descending diagonal
% from left to right is constant.
%
% Jina Lee, 12/05/2012
% -------------------------------------------------------------------------

% create the first col and row for constructing Toeplitz matrix
% col and row elements are independently drawn from N(0, 1)
c = randn(nrows, 1);
r = randn(1, ncols);

PsiM = toeplitz(c,r);

% Normalize each elem by 1/sqrt(ncols)
% JR_FIX: 2/6/13. Normalize the sampling matrix. PsiM is not all +/-1
% PsiM = PsiM .* (1/sqrt(ncols));
for i = 1 : size(PsiM, 1)                                     % JR_NORM_FIX
   denom = norm( PsiM(i, :) ) ;
   PsiM(i, :) = PsiM(i, :) / denom ;
end
end

% =========================================================================
function [PsiM] = circulantMatrix(nrows, ncols)
%circulantMatrix Construct a partial circulant matrix to be used in compressive sensing. 
%  PsiM = circulantMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. nrows are randomly selected from a
%  ncols*ncols circulant matrix.

%  Note on Circulant Matrix (source: Wikipedia): A circulant matrix is a
%  special kind of Toeplitz matrix where each row vector is rotated one
%  element to the right relative to the preceding row vector. In numerical
%  analysis, circulant matrices are important because they are diagonalized
%  by a discrete Fourier transform, and hence linear equations that contain
%  them may be quickly solved using a fast Fourier transform.
%
% Jina Lee, 12/05/2012
% -------------------------------------------------------------------------

% create the first row of circulant matrix
firstrow = randn(1, ncols) ;    

% indices of the first row and first col 
rind = 1:ncols;
cind =  ncols + 2 - rind' ;
cind(cind == (ncols+1)) = 1;   

% once the first row and column is given, just call toeplitz
circulantM = firstrow(toeplitz(cind,rind));

% create a nrows x 1 matrix containing indices of rows to be selected from
% circulantM
% JR_OLD_MATLAB. 2/4/13. My old matlab does not understand randperm(a, b).
% Doing it the inelegant way.
% inds = randperm( ncols, nrows );
inds = randperm( ncols ); inds = inds( 1 : nrows ) ;        % JR_OLD_MATLAB

% select nrows randomly from the ncols*ncols discrete Circulant matrix,
% creating an nrows x ncols matrix
PsiM = circulantM(inds, :) ; 

% Normalize each elem by 1/sqrt(ncols)
% JR_FIX 2/6/13. The elements of PsiM are not all +/-1 so normalize
% conventionally
% PsiM = PsiM .* (1/sqrt(ncols));
for i = 1 : size(PsiM, 1)                                     % JR_NORM_FIX
   denom = norm( PsiM(i, :) ) ;
   PsiM(i, :) = PsiM(i, :) / denom ;
end
end
% =========================================================================
function [PsiM] = noiseletMatrix(nrows, ncols)
%noiseletMatrix Construct a noiselet matrix to be used in compressive sensing. 
%  PsiM = noiseletMatrix(nrows, ncols) constructs a nrows X ncols random
%  matrix for compressive sensing. nrows are randomly selected from a
%  ncols*ncols noiselet matrix.
%
% we use Build_Noiselets() written by Laurent Duval licensed under 
% Creative Commons Attribution-Noncommercial-Share Alike 3.0 Unported License.
% http://creativecommons.org/licenses/by-nc-sa/3.0/
%
% Jina Lee, 01/21/2013
% -------------------------------------------------------------------------

% create a ncols x ncols Noiselet matrix
npower = logb(ncols,2);
noiseletM = Build_Noiselets(npower);

% Shuffle the rows of noiseletM and select the first nrows/2 rows. 
% Concatenate the real values of 1 through nrows/2 rows followed by the
% imaginary values. This creates a sensing matrix of 2*(nrows/2) x ncols
% size.   
shuffled_noiseletM = noiseletM(randperm(size(noiseletM,1)),:);
PsiM(1:floor(nrows/2),:) = real(shuffled_noiseletM(1:floor(nrows/2),:)); 

% when nrows is an odd number, match the dimension between two matrices
% if mod(nrows, 2) == 0
%     PsiM(floor(nrows/2)+1:nrows,:) = imag(shuffled_noiseletM(1:floor(nrows/2),:)); 
% else
    PsiM(floor(nrows/2)+1:nrows,:) = imag(shuffled_noiseletM(1:nrows-floor(nrows/2),:)); 
% end

% Normalize each elem by 1/sqrt(ncols)
% JR_FIX, 2/6/13. Normalize conventionally, PsiM may not be +/-1
% PsiM = PsiM .* (1/sqrt(ncols));
for i = 1 : size(PsiM, 1)                                     % JR_NORM_FIX
   denom = norm( PsiM(i, :) ) ;
   PsiM(i, :) = PsiM(i, :) / denom ;
end
end

% =========================================================================

% return the base B logarithm of X for B > 1
function Y = logb(X, B)
    assert(B > 1 , 'The base has to be greater than 1')
    Y = log(X) / log(B);
end

function [noiselet_Matrix] = Build_Noiselets(n_Power)

% [noiselet_Matrix] = Build_Noiselets(n_Power)
% Built an orthonormal noiselet matrix with 2^n_Power
% In:  n_Power (integer); default: n_Power = 5 (n_Sample = 2^n_Power = 32)
% Out: noiselet_Matrix (complex array); default (no outpu): display real and
% imaginery parts
% Example: noiselet_Matrix_1024 = Build_Noiselets(10);
%
% Comments: no optimization at all, suggestions welcome
% Creation: 2008/04/10
% Update:  2008/04/13
%
% Author: Laurent Duval
% URL: http://www.laurent-duval.eu

if nargin == 0
    n_Power = 5;
end
n_Sample = round(2^n_Power);

noiselet_Matrix = zeros(n_Sample,2*n_Sample-1);
noiselet_Matrix(:,1) = 1;
coef1 = 1 - i;
coef2 = 1 + i;
vect_Fill = zeros(floor(n_Sample/2),1);
for i_Col = 1:n_Sample-1
    vect_2x     = [noiselet_Matrix(1:2:n_Sample,i_Col);vect_Fill];
    vect_2x_1   = [vect_Fill;noiselet_Matrix(1:2:n_Sample,i_Col)];
    noiselet_Matrix(:,2*i_Col)   = coef1*vect_2x + coef2*vect_2x_1;
    noiselet_Matrix(:,2*i_Col+1) = coef2*vect_2x + coef1*vect_2x_1;
end
noiselet_Matrix = 1/n_Sample*noiselet_Matrix(:,n_Sample:2*n_Sample-1);

if nargout == 0
    figure(1)
    subplot(1,2,1)
    imagesc(real(noiselet_Matrix))
    xlabel('Real part')
    subplot(1,2,2)
    imagesc(imag(noiselet_Matrix))
    xlabel('Imaginery part')
end

end