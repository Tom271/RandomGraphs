%Find degree of a node

function deg=deg(G,i)
%Input graph and the node i for which you
%wish to find the degree
n=max(size(G));
k=spalloc(n,1,10);
for j=1:n
k(j)=G(i,j);
end
deg=full(sum(k));
end

%%Next work out exp deg, triangles, exp tri
