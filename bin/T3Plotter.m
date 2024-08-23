%Author: Bensingh Dhas
%Descripton: Plots the mesh with T3 element.
%Input: requires you to supply connectivity of the mesh and nodes
%associated with the connect.
%Remark: Crude code consider updating
%Created on: At the begining of time.

function T3Plotter(mesh,nodes,figId,c)
connect=[mesh.elements(:,1:end),mesh.elements(:,1)];
figure (figId)
hold on;
% for i=1:mesh.noEl
for i=1:mesh.noEl
    x=nodes(connect(i,:),:);
    plot(x(:,1),x(:,2),c,'linewidth',0.1);
    
end
% axis equal;
hold off;
end