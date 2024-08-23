function shear=computeShearForce(me,ori,shearDof, edgeIndex, gpDat)
feL=basisLag2D(gpDat.pt);
feTO=basis2D1P1M(gpDat.pt);
shear=zeros(me.noEl*gpDat.noGp,2);
itr=1;
for itEl=1:me.noEl
    elConn=me.elements(itEl,:);
    elNodes=me.nodes(elConn,:);
    oriSt=createBasisCorrection2D(ori(itEl,:), 'P1M');
    for itGp=1:gpDat.noGp
        dsp=feL.dsp(:,:,itGp);
        T=dsp*elNodes;
        spD=oriSt*feTO.phi(:,:,itGp);
        spD=(T\spD');
        shear(itr,:)=spD*shearDof(edgeIndex(itEl,:),:);
        itr=itr+1;
    end
end