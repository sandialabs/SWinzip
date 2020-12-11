function e=compute_error_field_alpert_dichotomic(pre,pr,fre,fr,Nc,th0)

% This function compares two data fields fre and fr represented on two
% meshes pre and pr, respectively. Since the fields might be large, the
% comparison takes place by chunks specified by spatially dividing the data
% using a dichotomic grouping technique.

% This method uses an iterative method to enforce wavelet coefficient signs

% pre, fre: mesh representing data field fre
% pr, fr: mesh representing data field fr
% Nc: number of chunks.
% th0: wavelet threshold multiplier (typical value is 0.001)

% e: is the error field

N=size(pr,1);
Ne=size(pre,1);
d=size(pr,2);

for j=1:400
    [parte,B,G] = dichotomic_grouping(pre', j^2, [num2str(d) 'axis']);
    if length(parte)==Nc
        break;
    end
end

for j=1:400
    [part,B,G] = dichotomic_grouping(pr', j^2, [num2str(d) 'axis']);
    if length(part)==Nc
        break;
    end
end

e=zeros(N,1);
mnx=prctile(fr,99)-prctile(fr,1);

parfor j=1:Nc
    if (length(parte{j})>0) && (length(part{j})>0)
        
        % we try several wavelet orders and pick the one resulting with
        % lowest error
            f0=zeros(length(part{j}),4);
            err=zeros(1,4);
            [f0(:,1),ww,err(1)]=interp_alpert(pre(parte{j},:),fre(parte{j}),pr(part{j},:),th0,2,mnx);
            [f0(:,2),ww,err(2)]=interp_alpert(pre(parte{j},:),fre(parte{j}),pr(part{j},:),th0,3,mnx);
            [f0(:,3),ww,err(3)]=interp_alpert(pre(parte{j},:),fre(parte{j}),pr(part{j},:),th0,4,mnx);
            [f0(:,4),ww,err(4)]=interp_alpert(pre(parte{j},:),fre(parte{j}),pr(part{j},:),th0,5,mnx);
            [em,m]=min(err);
            ej{j}=f0(:,m)-fr(part{j});
    end
end

for j=1:Nc
    e(part{j})=ej{j};
end

end
