function npti = make_table(nclp)
% function: uq_table
% purpose:  returns the samples multi-index table 
%           
% usage: npti = uq_order(nclp, ndim) 
%
% input: 
%    nclp: number of samples per dimension
%    ndim: number of stochastic dimensions
%
% output:
%    npti: multi-index matrix

ndim=length(nclp);
nt=prod(nclp);
npti=zeros(nt,ndim);
ind=zeros(1,ndim);
	
for j=1:ndim
  ind(j) = 1;
end

for m=1:nt

  if (m>1)
    
    ind(1) = ind(1)+1;

    ifl=0;
    if (ind(1)>nclp(1))
       ind(1)=1;
       ifl=1;
    end

    for j=2:ndim

       ind(j)=ind(j)+ifl;
       ifl=0;

       if(ind(j)>nclp(j))
          ind(j)=1;
          ifl=1;
       end

     end

  end   

  npti(m,:)=ind;
	  	  
end
