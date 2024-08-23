close all;clear all; clc
addpath('./bin');addpath('./mesh');

% %Material properties
mat.E=1.2e6; mat.h=0.1;           
mat.nu=0.0;
mat.r=1;
mat.la=mat.E*mat.nu/((1+mat.nu)*(1-2*mat.nu));
mat.mu=mat.E/2*(1+mat.nu);
mat.C=(mat.E*mat.h^3)/12;

%Geometry and Mesh details
mat.l=10; mat.b=1;
l=mat.l; b=mat.b;
nx=21; ny=3;
me=mesh2DTriangle([0,l],[0,b],nx,ny);
noEl=me.noEl;     noNd=me.noNd;
me.nodes=[me.nodes,zeros(noNd,1)];

[W,edgeIndex,edges]=createEdgeOrientation2D(me);
edgeTract=zeros(size(edgeIndex));

for i=1:noEl
    n=edgeIndex(i,:);
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

noEg=size(edges,1);
gpDat=gaussDataTri(2);

% Boundary Nodes
ls=0.000001;
leftNodes=find(me.nodes(:,1)<ls);
rightNodes=find(me.nodes(:,1)>(l-ls));
botNodes=find(me.nodes(:,2)<(ls));
topNodes=find(me.nodes(:,2)>(b-ls));

% Boundary Edges
edgeMidPt=zeros(noEg,2);
edgeMidPt(:,1)=0.5*(me.nodes(edges(:,1),1)+me.nodes(edges(:,2),1));
edgeMidPt(:,2)=0.5*(me.nodes(edges(:,1),2)+me.nodes(edges(:,2),2));
leftEdge=find(edgeMidPt(:,1)<0.000001);
rightEdge=find(edgeMidPt(:,1)>(l-0.000001));
noDof=4*noEg+5*me.noNd;

t1Dof=(1:noEg);
t2Dof=noEg+(1:noEg);
p1Dof=2*noEg+(1:noEg);
p2Dof=3*noEg+(1:noEg);
r1Dof=4*noEg+(1:noNd);
r2Dof=4*noEg+noNd+(1:noNd);
u1Dof=4*noEg+2*noNd+(1:noNd);
u2Dof=4*noEg+2*noNd+noNd+(1:noNd);
u3Dof=4*noEg+2*noNd+2*noNd+(1:noNd);

rDof=4*noEg+(1:2*noNd);
uDof=4*noEg+2*noNd+(1:3*noNd);

deadnode=leftNodes(:);
deadRotDoF=4*noEg+[deadnode,noNd+deadnode];

deadDispDoF=4*noEg+2*noNd+[deadnode,noNd+deadnode,2*noNd+deadnode];
deadDoF=[deadRotDoF, deadDispDoF];
actDof=setdiff(1:noDof,deadDoF)';

theta1 =zeros(noEg,1);
theta2 =zeros(noEg,1);
stress1=zeros(noEg,1);
stress2=zeros(noEg,1);
Rotdof =zeros(noNd,2);
displacement=zeros(noNd,3);

%Loading
fVec1=createBoundaryLoad_shear(me,edges(rightEdge,:),2);     %pulling
fVec2=createBoundaryLoad_moment(me,edges(rightEdge,:),0);    %moment

Shear=-mat.h^3*4000;
I=b*mat.h^3/12;
R=l/(2*pi);
Moment=(mat.E*I/R)*0;          
q=[0;0;0];

elRot=cell(3,1);
rotation=cell(me.noNd,1);
rotation_ref=cell(me.noNd,1);
for itNd=1:me.noNd
    rotation{itNd}=eye(3);
    rotation_ref{itNd}=eye(3);
end
noLs=4;
data=zeros(noLs,3);
for itLs=1:noLs
    err=100;
    itNr=1;
    disp('---------------------------')
    
    while err>1e-4
        
        fGlob=zeros(noDof,1);
        fGlob(rDof)=(itLs/noLs)*Moment*fVec2;
        fGlob(uDof)=(itLs/noLs)*Shear*fVec1;
        kGlob=zeros(noDof);
        sol=zeros(noDof,1);
        
        for i=1:me.noEl
            elConn=me.elements(i,:);
            elEdIdx=edgeIndex(i,:);
            elEdges=edges(elEdIdx,:);
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
            elRot_ref{1}=rotation_ref{elConn(1)};elRot_ref{2}=rotation_ref{elConn(2)};elRot_ref{3}=rotation_ref{elConn(3)};
            
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
        rotation_a=sol(rDof);
        Rotation1=rotation_a(noNd+1:end);
        displacement=displacement+reshape(sol(uDof),[me.noNd,3]);
        
        theta1=theta1+sol(t1Dof);
        theta2=theta2+sol(t2Dof);
        
        stress1=stress1+sol(p1Dof);
        stress2=stress2+sol(p2Dof);
        
        incRot=reshape(sol(rDof),[me.noNd,2]);
        
        rotation=updateRotationNodewise(me,rotation,incRot);
        err=norm(fGlob(actDof));
        msg=['itLs:' num2str(itLs) ' itNr:' num2str(itNr) ' Error:' num2str(err)];
        disp(msg);
        itNr=itNr+1;
    end
    
    deformation=me.nodes+displacement;
    data(itLs,:)=displacement(end-1,:);
end
plotterT6(me.elements,me.nodes,2,'k');
plotterT6(me.elements,deformation,2,'b');
xlim([-0.1*l,1.1*l]);
ylim([-0.1,b+0.2]);
zlim([-0.1*l,0.85*l]);
axis equal
hold off


