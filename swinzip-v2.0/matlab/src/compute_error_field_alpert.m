function e=compute_error_field_alpert(pre,pr,fre,fr,Nx,Ny,Nz,th0,vm)

% This function compares two data fields fre and fr represented on two
% meshes pre and pr, respectively. Since the fields might be large, the
% comparison takes place by chunks specified by spatially dividing the data
% in the x, y and z directions.

% pre, fre: mesh representing data field fre
% pr, fr: mesh representing data field fr
% Nx, Ny, Nz: number of chunks in the x, y and z directions.
% th0: wavelet threshold multiplier (typical value is 0.0075)
% vm: wavelet order (typical value is 5) 

% e: is the error field

N=size(pr,1);
d=size(pr,2);

xl=[min(pr(:,1)) max(pr(:,1))];
yl=[min(pr(:,2)) max(pr(:,2))];
if d==3
    zl=[min(pr(:,3)) max(pr(:,3))];
end

xx=xl(1):(xl(2)-xl(1))/Nx:xl(2);
yy=yl(1):(yl(2)-yl(1))/Ny:yl(2);
if d==3
    zz=zl(1):(zl(2)-zl(1))/Nz:zl(2);
end

e=zeros(N,1);

if d==2
    for i=1:Nx
        for j=1:Ny            
            ii=find( pr(:,1)>xx(i) & pr(:,1)<=xx(i+1) & pr(:,2)>yy(j) & pr(:,2)<=yy(j+1) );
            jj=find( pre(:,1)>xx(i) & pre(:,1)<=xx(i+1) & pre(:,2)>yy(j) & pre(:,2)<=yy(j+1) );

            if (length(ii)>0) && (length(jj)>0)
                f0=interp_field_alpert_chunk(pre(jj,1:2),fre(jj),pr(ii,1:2),th0,vm,max(fr)-min(fr));
                e(ii)=f0-fr(ii);
            end

        end
    end
else
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                ii=find( pr(:,1)>xx(i) & pr(:,1)<=xx(i+1) & pr(:,2)>yy(j) & pr(:,2)<=yy(j+1) & pr(:,3)>zz(k) & pr(:,3)<=zz(k+1));
                jj=find( pre(:,1)>xx(i) & pre(:,1)<=xx(i+1) & pre(:,2)>yy(j) & pre(:,2)<=yy(j+1) & pre(:,3)>zz(k) & pre(:,3)<=zz(k+1));                       
                
                if (length(ii)>0) && (length(jj)>0)
                    f0=interp_field_alpert_chunk(pre(jj,1:3),fre(jj),pr(ii,1:3),th0,vm,max(fr)-min(fr));
                    e(ii)=f0-fr(ii);
                end
            end

        end
    end
end
