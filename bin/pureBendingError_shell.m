function error=pureBendingError_shell(l,deformation,me)
error=0;
gpDat=gaussData2D(2);
feL=basisLag2D(gpDat.pt);
for itEl=1:me.noEl
    elConn=me.elements(itEl,:);
    elNodes=me.nodes(elConn,:);
    elDef=deformation(elConn,:);
    elLen=[norm(elNodes(1,:)) norm(elNodes(2,:)) norm(elNodes(3,:))];
    for itGp=1:gpDat.noGp
        sp=feL.sp(:,itGp);
        dsp=feL.dsp(:,:,itGp);
        T=norm(dsp*elNodes(:,1:2));
        JxW=T*gpDat.wt(itGp);
        rComp=sp'*elDef;
        s=sp'*elLen(:);
        rAna=pureBendingAnalytical(s,l);
        e=rAna-rComp;
        error=error+JxW*dot(e,e)^0.5;
    end
end

