function fr=interp_field_alpert(pre,pr,fre,Nx,Ny,Nz,th0,vm)

% This function interpolates a field fre represented on a mesh pre on another
% mesh pr. Since the fields might be large, the comparison takes place by
% chunks specified by spatially dividing the data in the x, y and z directions.

% pre, fre: mesh representing data field fre
% pr: mesh representing data field fr to be computed
% Nx, Ny, Nz: number of chunks in the x, y and z directions.
% th0: wavelet threshold multiplier (typical value is 0.0075)
% vm: wavelet order (typical value is 5) 

% fr: interpolated data field

N=size(pr,1);
d=size(pr,2);

% We first find the min/max range in the mesh
xl=[min([min(pre(:,1)) min(pr(:,1))]) max([max(pre(:,1)) max(pr(:,1))])];
yl=[min([min(pre(:,2)) min(pr(:,2))]) max([max(pre(:,2)) max(pr(:,2))])];
if d==3
    zl=[min([min(pre(:,3)) min(pr(:,3))]) max([max(pre(:,3)) max(pr(:,3))])];
end

% We divide the domain into chunks
dxl=xl(2)-xl(1);xx=xl(1):dxl/Nx:xl(2);
dyl=yl(2)-yl(1);yy=yl(1):dyl/Ny:yl(2);
if d==3
    dzl=zl(2)-zl(1);zz=zl(1):dzl/Nz:zl(2);
end

fr=zeros(N,1);

% We find the point indices in each chunk
m=1;
if d==2
    for i=1:Nx
        for j=1:Ny            
            ii=find( pr(:,1)>=xx(i) & pr(:,1)<=xx(i+1) & pr(:,2)>=yy(j) & pr(:,2)<=yy(j+1) );
            jj=find( pre(:,1)>=xx(i) & pre(:,1)<=xx(i+1) & pre(:,2)>=yy(j) & pre(:,2)<=yy(j+1) );
            if (length(ii)>0)
                I{m}=ii;
                J{m}=jj;
                m=m+1;
            end
        end
    end
else
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                ii=find( pr(:,1)>=xx(i) & pr(:,1)<=xx(i+1) & pr(:,2)>=yy(j) & pr(:,2)<=yy(j+1) & pr(:,3)>=zz(k) & pr(:,3)<=zz(k+1));
                jj=find( pre(:,1)>=xx(i) & pre(:,1)<=xx(i+1) & pre(:,2)>=yy(j) & pre(:,2)<=yy(j+1) & pre(:,3)>=zz(k) & pre(:,3)<=zz(k+1));
                if (length(ii)>0)
                    I{m}=ii;
                    J{m}=jj;
                    m=m+1;
                end
            end
        end
    end
end

% Adjust some out-of-range indices
m=length(I);
for i=1:m
   if isempty(J{i})
       xll=[min(pr(I{i},1)) max(pr(I{i},1))];
       yll=[min(pr(I{i},2)) max(pr(I{i},2))];
       jj=find( pre(:,1)>=xll(1)-dxl/30 & pre(:,1)<=xll(2)+dxl/30 & pre(:,2)>=yll(1)-dyl/30 & pre(:,2)<=yll(2)+dyl/30);
       if (d==3)
           zll=[min(pr(I{i},3)) max(pr(I{i},3))];
           jj=find( pre(:,1)>=xll(1)-dxl/30 & pre(:,1)<=xll(2)+dxl/30 & pre(:,2)>=yll(1)-dyl/30 & pre(:,2)<=yll(2)+dyl/30 & pre(:,3)>=zll(1)-dzl/30 & pre(:,3)<=zll(2)+dzl/30);
       end
       J{i}=jj;
   end
end

% Interpolate each chunk and for the whole interpolated field
mxmn=max(fre)-min(fre);
parfor i=1:m
    f0{i}=interp_field_alpert_chunk(pre(J{i},:),fre(J{i}),pr(I{i},:),th0,vm,mxmn);
end

for i=1:m
   fr(I{i})=f0{i};
end