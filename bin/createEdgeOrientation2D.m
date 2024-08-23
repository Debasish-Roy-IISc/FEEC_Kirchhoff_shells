function [W,edgeIndex,edges]=createEdgeOrientation2D(me)
graphMe = findEdges(me.elements,me.noNd); %find edges from graph obtained by adjmatrix,

edges=graphMe.Edges.EndNodes; % edges from the graph created from adjmatrix
edgeIndex=zeros(me.noEl,3);
for i=1:me.noEl
    elConn1=me.elements(i,1);  % 1   %Element connectivity from given graph ie node numbers of element
    elConn2=me.elements(i,2);  % 2
    elConn3=me.elements(i,3);  % 3
    elConn=[elConn1  elConn1 elConn2]; %1 1 2 
    tmp=[elConn2 elConn3 elConn3 ];    % 2 3 3 
    elEdges=[elConn(:),tmp(:)]; 
    edgeIndex(i,:)= findedge(graphMe,elEdges(:,1),elEdges(:,2)); 
end

W=zeros(me.noEl,3);
for i=1:me.noEl
    elConn1=me.elements(i,1);  % 1   
    elConn2=me.elements(i,2);  % 2
    elConn3=me.elements(i,3);  % 3
    elConn=[elConn1  elConn1 elConn2]; %1 1 2 
    tmp=[elConn2 elConn3 elConn3 ];      % 2 3 3 
    elOriEdges=[elConn(:) tmp(:)];  
    
    elEdgeIdx=edgeIndex(i,:);
    elEdges=edges(elEdgeIdx,:);
    o=zeros(1,3);
    for j=1:3
        if elOriEdges(j,1)==elEdges(j,1)
            o(j)=1;
        else
            o(j)=-1;
        end
    end
    W(i,:)=o;
end