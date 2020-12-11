clear all

N=40;
dd=3;
Nv=8;

for i=1:Nv
    ax(1,i)=subplot(2,Nv,i);hold
    ax(2,i)=subplot(2,Nv,i+Nv);hold
end

m=importdata(['../../data/channel/ksgs000' num2str(1) '.txt'],' ',1);m=m.data;
mx=m(1,dd+1:end);
mn=mx;

for i=1:N
    i
    if i<10
        d=importdata(['out_r000' num2str(i) '.txt'],' ',1);d=d.data;
        m=importdata(['../../data/channel/ksgs000' num2str(i) '.txt'],' ',1);m=m.data;
    else
        d=importdata(['out_r00' num2str(i) '.txt'],' ',1);d=d.data;
        m=importdata(['../../data/channel/ksgs00' num2str(i) '.txt'],' ',1);m=m.data;
    end
    
    Np=size(m,1);
    S = repmat([5],1,Np);
    
    for j=1:Nv
        if prctile(m(:,j+dd),80)>=mx(j)
            mx(j)=prctile(m(:,j+dd),80);
        end

        if prctile(m(:,j+dd),20)<=mn(j)
            mn(j)=prctile(m(:,j+dd),20);
        end               
    end
    
    for j=1:Nv
        set(gcf,'CurrentAxes',ax(1,j))
        scatter3(m(:,2),m(:,1),m(:,3),S,m(:,j+dd),'filled')

        set(gcf,'CurrentAxes',ax(2,j))
        scatter3(m(:,2),m(:,1),m(:,3),S,d(:,j),'filled')
    end
        
end

for j=1:Nv
    set(gcf,'CurrentAxes',ax(1,j))
    caxis([mn(j) mx(j)])
    axis equal;
    axis off
    view(3)
    set(gcf,'CurrentAxes',ax(2,j))
    caxis([mn(j) mx(j)])
    axis equal;
    axis off
    view(3)
end
