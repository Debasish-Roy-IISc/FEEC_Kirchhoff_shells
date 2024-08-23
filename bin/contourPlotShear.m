function contourPlotShear(elements, nodes, shear)
figure()
hold on;
noEl=size(elements,1);
for itEl=1:noEl
    elNodes=nodes(elements(itEl,:),:);
    fill(elNodes(:,1), elNodes(:,2), shear(elements(itEl,:),1),'LineStyle','none');
end
hold off;
end