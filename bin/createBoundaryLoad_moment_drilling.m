function fVec=createBoundaryLoad_moment_drilling(me,bEdges,dir)
fVec=zeros(3*me.noNd,1);
for i=1:size(bEdges,1)
    cEdge=bEdges(i,:);
    nodes=me.nodes(cEdge,:);
    len=norm(nodes(2,:)-nodes(1,:));
    dof=(dir*me.noNd)+cEdge;
    load=0.5*len*[1;1];
    fVec(dof)= fVec(dof)+load;
end