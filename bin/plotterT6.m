function plotterT6(elements,nodes,figId,c)
figure (figId)
hold on;
noEl=size(elements,1);
for itEl=1:noEl
    elConn=elements(itEl,:);
    tmp=[elConn(1) elConn(2);
        elConn(2) elConn(3);
        elConn(3) elConn(1);];
    for itEg=1:3
        x=nodes(tmp(itEg,:),:);
        plot3(x(:,1),x(:,2),x(:,3),c,'linewidth',0.1);
    end
end
xlabel('X');
ylabel('Y');
zlabel('Z');

axis equal;
hold off;
end
