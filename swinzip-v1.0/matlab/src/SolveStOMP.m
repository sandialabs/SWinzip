function [sol, iter, I, x_I, activeSet] = SolveStOMP(A, y, maxIters, verbose, OptTol, x_I0, activeSet0)
% SolveStOMP: Implementation of Iterative Threshold-Selective Projection
% algorithm 

% Inputs
%	A           The matrix (of size n x N) equal to the product of the sampling and wavelet
%	            matrices

%	y           vector of samples of length n. 

%	maxIters    maximum number of iterations

%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output

%	OptTol      Error tolerance

%   x_I0        Initial solution vector

%   activeSet0  Initial indices of the solution vector



% Outputs
%	 sol        Solution of StOMP.

%	 iter       Total number of steps taken.

%    I          Number of active coefficients

%    x_I        Solution vector

%   activeSet   Indices of the solution vector

% This is an enhanced version of the original StOMP algortihm developed by Donoho
% et.al. found in the SparseLab library https://sparselab.stanford.edu/

n = length(y);
N = size(A,2);

% Initialize
if (length(x_I0)==0)
    rr=0;
    normy = norm(y);
else
    rr=A(:, activeSet0)*x_I0;
    normy = norm(y);
end

% The thresholding parameter is automatically computed by fitting it to a
% norm metric of the samples vector
thr=3*normy/mean(abs(y))/sqrt(size(A,2))+0.03;
thr=thr/sqrt(n);

if (verbose)
    fprintf('Iteration 1: |I| = %d, ||r||_2 = %g\n', length(activeSet0), normy);
end

iter=1;
done=0;

res=y-rr;
activeSet = activeSet0;
x_I = x_I0;
r=norm(res);

r=norm(res);
while ~done
    % Compute residual correlations
    corr = A'*res/r;

    % Apply hard thresholding
    I = find(abs(corr) > thr);

    % If no new entries, we are done
    J = union(activeSet, I);
    if (length(J) == length(activeSet)) 
        done = 1;
    else
        activeSet = J;

        % Compute current estimate and residual
        Aactive = A(:, activeSet);
        B=Aactive'*Aactive;
        
        % Depending on the condition of the matrix B, we use a suitable
        % method to compute x_I
        if (rcond(B)>eps)
            x_I = inv(B)*(Aactive'*y);
        else
            x_I = Aactive \ y;
        end
      
        res = y - Aactive*x_I;

        if (verbose)
            fprintf('Iteration %d: |I| = %d, ||r||_2 = %g\n', iter, length(activeSet), r);
        end
    end  

    r=norm(res);
    
    iter = iter+1;

    % Check stopping criteria
    if (iter > maxIters)||(r <= OptTol*normy)
        done = 1;
    end
   
end

I=length(activeSet);
sol = zeros(N,1);
sol(activeSet) = x_I;

end


