function [W,edgeIndex,edges]=createEdgeOrientation(me)
graphMe = findEdges(me.elements,me.noNd);
edges=graphMe.Edges.EndNodes;
edgeIndex=zeros(me.noEl,3);
for itEl=1:me.noEl
    elConn=me.elements(itEl,:);
    tmp=[elConn(2:3) elConn(1)];
    elEdges=[elConn(:) tmp(:)];
    edgeIndex(itEl,:)=findedge(graphMe,elEdges(:,1),elEdges(:,2));
end
% 
% W=zeros(me.noNd,me.noNd);
% for itEl=1:me.noEl
%     elConn=me.elements(itEl,:);
%     tmp=[elConn(2:3) elConn(1)];
%     elOriEdges=[elConn(:) tmp(:)];
%     for i=1:3
%         id=elOriEdges(i,:);
%         o=W(id(1),id(2));
%         if abs(o)<0.001
%             W(id(1),id(2))=1;
%             W(id(2),id(1))=-1;
%         end
%     end
% end



W=zeros(me.noEl,3);
for itEl=1:me.noEl
    elConn=me.elements(itEl,:);
    tmp=[elConn(2:3) elConn(1)];
    elOriEdges=[elConn(:) tmp(:)];

    elEdgeIdx=edgeIndex(itEl,:);
    elEdges=edges(elEdgeIdx,:);
    o=zeros(1,3);
    for i=1:3
        if elOriEdges(i,1)==elEdges(i,1)
            o(i)=1;
        else
            o(i)=-1;
        end
    end
    W(itEl,:)=o;
end
