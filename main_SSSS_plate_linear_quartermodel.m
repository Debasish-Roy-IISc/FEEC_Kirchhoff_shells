% Plate simply supported on all edges and subjected to central point load
close all;clear all; clc
addpath('./bin');
% %Material properties
mat.E=1.2e3; mat.h=0.1;
mat.nu=0.3;
mat.r=1;
mat.la=mat.E*mat.nu/((1+mat.nu)*(1-2*mat.nu));
mat.mu=mat.E/2*(1+mat.nu);
mat.C=(mat.E*mat.h^3)/(12*(1-mat.nu^2));

%Geometry and Mesh details
mat.l=5; mat.b=5;
l=mat.l; b=mat.b;
nx=7; ny=7;
me=mesh2DTriangle([0,l],[0,b],nx,ny);
noEl=me.noEl;     noNd=me.noNd;
me.nodes=[me.nodes,zeros(noNd,1)];

[W,edgeIndex,edges]=createEdgeOrientation2D(me);
[W1,edgeIndex1,edges1]=createEdgeOrientation2D(me);
edgeTract=zeros(size(edgeIndex1));
for i=1:noEl
    n=edgeIndex1(i,:);
   for j=1:3
       w=edges(n(j),:);
       x=me.nodes(w(1),:);
       y=me.nodes(w(2),:);
       if x(1)==y(1)
          edgeTract(i,j)=1;
       elseif x(2)==y(2)
            edgeTract(i,j)=2;
       else
           edgeTract(i,j)=3;
       end
   end
end
noLs=1;

noEg=size(edges,1);
gpDat=gaussDataTri(2);

% Boundary Nodes
ls=0.000001;
leftNodes=find(me.nodes(:,1)<ls);
rightNodes=find(me.nodes(:,1)>(l-ls));
botNodes=find(me.nodes(:,2)<(ls));
topNodes=find(me.nodes(:,2)>(b-ls));
loadNodes=find((me.nodes(:,1)<(0.5*l+ls)) & (me.nodes(:,1)>(0.5*l-ls)));

% Boundary Edges
edgeMidPt=zeros(noEg,2);
edgeMidPt(:,1)=0.5*(me.nodes(edges(:,1),1)+me.nodes(edges(:,2),1));
edgeMidPt(:,2)=0.5*(me.nodes(edges(:,1),2)+me.nodes(edges(:,2),2));
leftEdge=find(edgeMidPt(:,1)<0.000001);
rightEdge=find(edgeMidPt(:,1)>(l-0.000001));
loadEdge=find((edgeMidPt(:,1)<(0.5*l+ls)) & (edgeMidPt(:,1)>(0.5*l-ls)));
noDof=4*noEg+5*me.noNd;

t1Dof=(1:noEg);
t2Dof=noEg+(1:noEg);
rDof=4*noEg+(1:2*noNd);
r1Dof=4*noEg+(1:noNd);
r2Dof=4*noEg+noNd+(1:noNd);
u1Dof=4*noEg+2*noNd+(1:noNd);
u2Dof=4*noEg+2*noNd+noNd+(1:noNd);
u3Dof=4*noEg+2*noNd+2*noNd+(1:noNd);
uDof=4*noEg+2*noNd+(1:3*noNd);

deadnode=[leftNodes(:);botNodes(:)];
deadDispDoF=4*noEg+[2*noNd+[deadnode(:);rightNodes(:)];3*noNd+[deadnode(:); topNodes(:)];4*noNd+deadnode(:)];
deadRotDoF=4*noEg+[rightNodes(:);noNd+topNodes(:)];
loadNode=4*noEg+4*noNd+find(((me.nodes(:,1)<(l+ls)) & (me.nodes(:,1)>(l-ls))) & ((me.nodes(:,2)<(b+ls)) & (me.nodes(:,2)>(b-ls))));

deadDoF=[deadRotDoF(:); deadDispDoF(:)];
actDof=setdiff(1:noDof,deadDoF)';

theta1 =zeros(noEg,1);
theta2 =zeros(noEg,1);
stress1=zeros(noEg,1);
stress2=zeros(noEg,1);
Rotdof =zeros(noNd,2);
displacement=zeros(noNd,3);

fVec1=createBoundaryLoad_shear(me,edges(loadEdge,:),2);     %pulling
fVec2=createBoundaryLoad_moment(me,edges(rightEdge,:),1);    %moment

Shear=-20*mat.h^3/4;
Moment=60*mat.h^3*0;

elRot=cell(3,1);
rotation=cell(me.noNd,1);
for itNd=1:me.noNd
    rotation{itNd}=eye(3);
    rotation_ref{itNd}=eye(3);
end
q=[0;0;0];
for itLs=1:noLs
    itNr=1;
    disp('---------------------------')
    
    fGlob=zeros(noDof,1);
    fGlob(rDof)=(itLs/noLs)*Moment*fVec2;
    fGlob(loadNode)=(itLs/noLs)*Shear;
    
    kGlob=zeros(noDof);
    sol=zeros(noDof,1);
    
    for i=1:me.noEl
        elConn=me.elements(i,:);
        elEdIdx=edgeIndex(i,:);
        elEdIdxTr=edgeTract(i,:);
        elEdIdxTra=[elEdIdx(elEdIdxTr(1));elEdIdx(elEdIdxTr(2));elEdIdx(elEdIdxTr(3))]';
        elNodes=me.nodes(elConn,:);
        
        elTheta1=theta1(elEdIdx,:);
        elTheta2=theta2(elEdIdx,:);
        elStress1=stress1(elEdIdxTra,:);
        elStress2=stress2(elEdIdxTra,:);
        elRotdof=Rotdof(elConn,:);
        elDisp=displacement(elConn,:);
        elRot{1}=rotation{elConn(1)};elRot{2}=rotation{elConn(2)};elRot{3}=rotation{elConn(3)};
        
        [f,k]=force(elNodes, elTheta1, elTheta2,...
            elStress1,elStress2,elRot,rotation_ref, elDisp, W(i,:), gpDat, mat,q);
        
        t=[elEdIdx(:);noEg+elEdIdx(:)];
        p=2*noEg+[elEdIdxTra(:);noEg+elEdIdxTra(:)];
        r=(4*noEg)+[elConn(:);me.noNd+elConn(:)];
        u=(4*noEg)+(2*noNd)+[elConn(:);me.noNd+elConn(:);2*me.noNd+elConn(:)];
        
        dof=[t;p;r;u];
        kGlob(dof,dof)=kGlob(dof,dof)+k;
        fGlob(dof)=fGlob(dof)+f;
    end
    sol(actDof)=-kGlob(actDof,actDof)\fGlob(actDof);
    rotation=sol(rDof);
    Rotation1=rotation(noNd+1:end);
    displacement=displacement+reshape(sol(uDof),[me.noNd,3]);
    
    err=norm(fGlob(actDof));
    msg=['itLs:' num2str(itLs) ' itNr:' num2str(itNr) ' Error:' num2str(err)];
    disp(msg);
    itNr=itNr+1;
    
    nodes1=[me.nodes];
    deformation=nodes1+displacement;
    
end
plotterT6(me.elements,nodes1,2,'k');
plotterT6(me.elements,deformation,2,'b');
present_mid_point_displacement=-max(displacement(:,3))
   l=2*mat.l; b=2*mat.b;
   Shear=4*Shear;
%J N Reddy : Theory and analysis of elastic plates and shell  page 231
Analytical1=(0.0116*Shear*l^2)/(mat.C)        %for mu=0.3,  refer table 7.3.1

Analytical2=(1.1600*Shear*l^2)/(100*mat.C)     %for l/h=100  1.1600
