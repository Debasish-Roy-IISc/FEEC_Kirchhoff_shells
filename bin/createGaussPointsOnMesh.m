function meshGp=createGaussPointsOnMesh(me, gaussData)
meshGp=zeros(me.noEl*gaussData.noGp,2);
feL=basisLag2D(gaussData.pt);
itr=1;
for itEl=1:me.noEl
    elConn=me.elements(itEl,:);
    elNodes=me.nodes(elConn,:);
    for itGp=1:gaussData.noGp
        meshGp(itr,:)=feL.sp(:,itGp)'*elNodes;
        itr=itr+1;
    end
end