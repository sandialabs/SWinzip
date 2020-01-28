function p=get_mesh(xyl,nl)
nl2=prod(nl);
[x,y]=ndgrid(fliplr(linspace(xyl(1,1),xyl(1,2),nl(1))),fliplr(linspace(xyl(2,1),xyl(2,2),nl(2))));
p=[reshape(x,1,nl2);reshape(y,1,nl2)];