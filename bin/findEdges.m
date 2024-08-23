function G = findEdges(elements,noNd)
el= size(elements);
% if el(1,2)==4
% gEdge=[me.elements(:,1) me.elements(:,2);
%     me.elements(:,2) me.elements(:,3);
%     me.elements(:,3) me.elements(:,4);
%     me.elements(:,4) me.elements(:,1);];
% end
% if el(1,2)==3
% gEdge=[me.elements(:,1) me.elements(:,2);
%     me.elements(:,2) me.elements(:,3);
%     me.elements(:,3) me.elements(:,1)];
% end
gEdge=[];
for i=1:(el(1,2)-1)
    gEdge=[gEdge;elements(:,i) elements(:,i+1)];
end
gEdge=[gEdge;elements(:,end) elements(:,1)];

adjMatrix=zeros(noNd,noNd);
for i=1:numel(gEdge(:,1))
%     d=norm(me.nodes(gEdge(i,2))-me.nodes(gEdge(i,1)));
    adjMatrix(gEdge(i,1),gEdge(i,2))=1;
    adjMatrix(gEdge(i,2),gEdge(i,1))=1;
end
% adjMatrix=adjMatrix+adjMatrix';
% adjMatrix(adjMatrix>0)=1;
G=graph(adjMatrix);
% edges=G.Edges.EndNodes;
end