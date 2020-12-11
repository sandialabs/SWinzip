function e=compute_error_field_alpert_dichotomic_V(pre,pr,fre,fr,Nc)

% This function compares two data fields fre and fr represented on two
% meshes pre and pr, respectively. Since the fields might be large, the
% comparison takes place by chunks specified by spatially dividing the data
% using a dichotomic grouping technique.

% This method uses a direct method to enforce wavelet coefficient signs

% pre, fre: mesh representing data field fre
% pr, fr: mesh representing data field fr
% Nc: number of chunks.

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
        e0=zeros(length(part{j}),4);
        err=zeros(1,4);
        [e0(:,1),err(1)]=interp_alpert_V(pre(parte{j},:)',pr(part{j},:)',fre(parte{j}),fr(part{j}),2);
        [e0(:,2),err(2)]=interp_alpert_V(pre(parte{j},:)',pr(part{j},:)',fre(parte{j}),fr(part{j}),3);
        [e0(:,3),err(3)]=interp_alpert_V(pre(parte{j},:)',pr(part{j},:)',fre(parte{j}),fr(part{j}),4);
        [e0(:,4),err(4)]=interp_alpert_V(pre(parte{j},:)',pr(part{j},:)',fre(parte{j}),fr(part{j}),5);

        [em,m]=min(err);
        ej{j}=e0(:,m);
    end
end

for j=1:Nc
    e(part{j})=ej{j};
end

end
