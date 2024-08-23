function [f,k]=force(elNodes, elTheta1, elTheta2,...
    stress1,stress2,elRot,elRot_ref, elDisp, ori, gpDat, mat,q)
k=zeros(27);
f=zeros(27,1);
E=mat.E;h=mat.h;nu=mat.nu ;

Cm=E*h/(4*(1-nu^2));
Cb=(E*h^3)/(12*(1-nu^2));

E11=[1;0;0];E22=[0;1;0]; E33=[0;0;1];  e1=[1;0];e2=[0;1];

feL=basisLag2D(gpDat.pt);
feTO=basis2D1P1M(gpDat.pt);

R1=elRot{1};

ta=(elRot{1}(:,3))'; tb=(elRot{2}(:,3))' ;tc=(elRot{3}(:,3))';   %Spatial Directors at the nodes
incR1=rotationBetweenUnitVec(ta,tb);
incR2=rotationBetweenUnitVec(ta,tc);
psi=zeros(3,1);
psi(:,2)=skew2axial(logMap(incR1));          %axial vector of incremental spatial rotation
psi(:,3)=skew2axial(logMap(incR2));

oriSt=createBasisCorrection2D(ori, 'P1M');

t1Dof=1:3; t2Dof=4:6;
p1Dof=7:9; p2Dof=10:12;
r1Dof=13:15; r2Dof=16:18;
u1Dof=19:21; u2Dof=22:24; u3Dof=25:27;

for itGp=1:gpDat.noGp
    sp=feL.sp(:,itGp);
    dsp=feL.dsp(:,:,itGp);
    T=dsp*elNodes(:,1:2);
    JxW=0.5*det(T)*gpDat.wt(itGp);
    spD=oriSt*feTO.phi(:,:,itGp);
    dsp=T\dsp;                 %Coordinate transform for shape functions
    spD=T\spD';

    db0=(dsp*elNodes);
    b1_0=db0(1,:);
    b2_0=db0(2,:);

    du=dsp*elDisp(:,1)+dsp*elNodes(:,1);
    dv=dsp*elDisp(:,2)+dsp*elNodes(:,2);
    dw=dsp*elDisp(:,3)+dsp*elNodes(:,3);

    t1=spD*elTheta1(:,1);
    t2=spD*elTheta2(:,1);

    p1=[E11,E22,E33]*stress1;  p2=[E11,E22,E33]*stress2;
    t1(1)=1+t1(1); t2(2)=1+t2(2);

    psigp=psi*sp;
    W=axial2Skew(psigp);
    t=(rodriguesFormula(W)*ta')';
    R=rodriguesFormula(W)*R1;             %Rotation matrix

    psin(:,1)=cross(psi(:,1),t');
    psin(:,2)=cross(psi(:,2),t');
    psin(:,3)=cross(psi(:,3),t');
    dt=psin*dsp';
    E1=R*E11;
    E2=R*E22;

    Me0=[b1_0*b1_0',b1_0*b2_0';
        b2_0*b1_0',b2_0*b2_0'];           %metric: reference config'n

    Hs=[Me0(1,1)*Me0(1,1),   nu*Me0(1,1)*Me0(2,2)+(1-nu)*Me0(1,2)*Me0(1,2),   2*Me0(1,1)*Me0(1,2);
        nu*Me0(1,1)*Me0(2,2)+(1-nu)*Me0(1,2)*Me0(1,2),  Me0(2,2)*Me0(2,2)  ,  2*Me0(2,2)*Me0(1,2);
        2*Me0(1,1)*Me0(1,2),   2*Me0(2,2)*Me0(1,2),   2*(1+nu)*(Me0(1,2)*Me0(1,2))+2*(1-nu)*Me0(1,1)*Me0(2,2)];

    spDdu=wedgeOneForms2D1(spD',du');
    spDdv=wedgeOneForms2D1(spD',dv');
    spDdw=wedgeOneForms2D1(spD',dw');
    spDt1=wedgeOneForms2D1(spD',t1');
    spDt2=wedgeOneForms2D1(spD',t2');
    dspt1=wedgeOneForms2D1(dsp',t1');
    dspt2=wedgeOneForms2D1(dsp',t2');

    spD_spD=(wedgeOneForms2D1(spD',spD'));
    spD_dsp=(wedgeOneForms2D1(spD',dsp'));

    t1_du=det([t1 du]);    t1_dv=det([t1 dv]) ;     t1_dw=det([t1 dw]);
    t2_du=det([t2 du]);    t2_dv=det([t2 dv]) ;     t2_dw=det([t2 dw]);
    t1_t2=det([t1 t2]);
    t1u=[t1_du, t1_dv ,t1_dw];
    t2u=[t2_du,t2_dv,t2_dw];
    spD_u=[spDdu, spDdv ,spDdw]';        

    %connection coefficients at current state  and derivatives of curvature

    O13=dt'*E1;                            %de3*e1
    O23=dt'*E2;                            %de3*e2
    k11=O13(1)*t1(1)+O23(1)*t2(1);
    k22=O13(2)*t1(2)+O23(2)*t2(2);         %     k11=b1*t_1'; k22=b2*t_2'; k12=b1*t_2'; k21=b2*t_1';
    k12=O13(1)*t1(2)+O23(1)*t2(2);
    k21=O13(2)*t1(1)+O23(2)*t2(1);
    K=[k11,k22,k12,k21];
    Ks=[k11,k22,0.5*(k12+k21)];                    %Curvature vector


    W1=axial2Skew(E11)   ;  W2=axial2Skew(E22);
    Dw1R=R*W2;    Dw2R=-R*W1;
    Dw1RtdR=W2;    Dw2RtdR=-W1;
    Dw1w1R=R*W2*W2;     Dw2w2R=R*W1*W1;
    Dw1w2R=-0.5*R*(W1*W2+W2*W1);

    Dw1w1RtdR_a=R'*(-W2*W2)*R+R'*(W2*W2)*R;       %for sp*dsp terms
    Dw2w2RtdR_a=R'*(-W1*W1)*R+R'*(W1*W1)*R;
    Dw1w2RtdR_a=0*(R'*(-W1*W2)*R+R'*(W2*W1)*R);   %Checked

    Dw1w1RtdR_b=R'*(-W2*W2)*R+R'*(W2*W2)*R;       %for dsp*sp terms
    Dw2w2RtdR_b=R'*(-W1*W1)*R+R'*(W1*W1)*R;
    Dw1w2RtdR_b=0*(R'*(-W2*W1)*R+R'*(W1*W2)*R);

    Dw1w1_k11=E11'*Dw1w1RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e1)')+E22'*Dw1w1RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e1)')...
        +E11'*Dw1w1RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e1)')'+E22'*Dw1w1RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e1)')';
    Dw1w1_k22=E11'*Dw1w1RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e2)')+E22'*Dw1w1RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e2)')+...
        E11'*Dw1w1RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e2)')'+E22'*Dw1w1RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e2)')';
    Dw1w1_k12=E11'*Dw1w1RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e1)')+E22'*Dw1w1RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e1)')+...
        E11'*Dw1w1RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e1)')'+E22'*Dw1w1RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e1)')';
    Dw1w1_k21=E11'*Dw1w1RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e2)')+E22'*Dw1w1RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e2)')+...
        E11'*Dw1w1RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e2)')'+E22'*Dw1w1RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e2)')';

    Dw2w2_k11=E11'*Dw2w2RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e1)')+E22'*Dw2w2RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e1)')+...
        E11'*Dw2w2RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e1)')'+E22'*Dw2w2RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e1)')';
    Dw2w2_k22=E11'*Dw2w2RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e2)')+E22'*Dw2w2RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e2)')+...
        E11'*Dw2w2RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e2)')'+E22'*Dw2w2RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e2)')';
    Dw2w2_k12=E11'*Dw2w2RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e1)')+E22'*Dw2w2RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e1)')+...
        E11'*Dw2w2RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e1)')'+E22'*Dw2w2RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e1)')';
    Dw2w2_k21=E11'*Dw2w2RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e2)')+E22'*Dw2w2RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e2)')+...
        E11'*Dw2w2RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e2)')'+E22'*Dw2w2RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e2)')';

    Dw1w2_k11=E11'*Dw1w2RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e1)')+E22'*Dw1w2RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e1)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e1)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e1)')';
    Dw1w2_k22=E11'*Dw1w2RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e2)')+E22'*Dw1w2RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e2)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e2)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e2)')';
    Dw1w2_k12=E11'*Dw1w2RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e1)')+E22'*Dw1w2RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e1)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e1)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e1)')';
    Dw1w2_k21=E11'*Dw1w2RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e2)')+E22'*Dw1w2RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e2)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e2)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e2)')';

    Dw2w1_k12=E11'*Dw1w2RtdR_a*E33*(t1'*e2)*(sp*(dsp'*e1)')+E22'*Dw1w2RtdR_a*E33*(t2'*e2)*(sp*(dsp'*e1)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e2)*(sp*(dsp'*e1)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e2)*(sp*(dsp'*e1)')';
    Dw2w1_k21=E11'*Dw1w2RtdR_a*E33*(t1'*e1)*(sp*(dsp'*e2)')+E22'*Dw1w2RtdR_a*E33*(t2'*e1)*(sp*(dsp'*e2)')+...
        E11'*Dw1w2RtdR_b*E33*(t1'*e1)*(sp*(dsp'*e2)')'+E22'*Dw1w2RtdR_b*E33*(t2'*e1)*(sp*(dsp'*e2)')';


    Dw1w1_k12s=0.5*(Dw1w1_k12+Dw1w1_k21);
    Dw2w2_k12s=0.5*(Dw2w2_k12+Dw2w2_k21);
    Dw1w2_k12s=0.5*(Dw1w2_k12+Dw1w2_k21);

    Dw1_k11=(E11'*Dw1RtdR*E33)*(t1'*e1)*(dsp'*e1)+(E22'*Dw1RtdR*E33)*(t2'*e1)*(dsp'*e1);
    Dw1_k22=(E11'*Dw1RtdR*E33)*(t1'*e2)*(dsp'*e2)+(E22'*Dw1RtdR*E33)*(t2'*e2)*(dsp'*e2);
    Dw1_k12=(E11'*Dw1RtdR*E33)*(t1'*e2)*(dsp'*e1)+(E22'*Dw1RtdR*E33)*(t2'*e2)*(dsp'*e1);
    Dw1_k21=(E11'*Dw1RtdR*E33)*(t1'*e1)*(dsp'*e2)+(E22'*Dw1RtdR*E33)*(t2'*e1)*(dsp'*e2);
    Dw1_Ks=[Dw1_k11,Dw1_k22,0.5*(Dw1_k12+Dw1_k21)];

    Dw2_k11=(E11'*Dw2RtdR*E33)*(t1'*e1)*(dsp'*e1)+(E22'*Dw2RtdR*E33)*(t2'*e1)*(dsp'*e1);
    Dw2_k22=(E11'*Dw2RtdR*E33)*(t1'*e2)*(dsp'*e2)+(E22'*Dw2RtdR*E33)*(t2'*e2)*(dsp'*e2);
    Dw2_k12=(E11'*Dw2RtdR*E33)*(t1'*e2)*(dsp'*e1)+(E22'*Dw2RtdR*E33)*(t2'*e2)*(dsp'*e1);
    Dw2_k21=(E11'*Dw2RtdR*E33)*(t1'*e1)*(dsp'*e2)+(E22'*Dw2RtdR*E33)*(t2'*e1)*(dsp'*e2);
    Dw2_Ks=[Dw2_k11,Dw2_k22,0.5*(Dw2_k12+Dw2_k21)];

    Dt1_k11=O13(1)*(spD'*e1);
    Dt1_k22=O13(2)*(spD'*e2);
    Dt1_k12=O13(1)*(spD'*e2);
    Dt1_k21=O13(2)*(spD'*e1);
    Dt1_Ks=[Dt1_k11,Dt1_k22,0.5*(Dt1_k12+Dt1_k21)];

    Dt2_k11=O23(1)*(spD'*e1);
    Dt2_k22=O23(2)*(spD'*e2);
    Dt2_k12=O23(1)*(spD'*e2);
    Dt2_k21=O23(2)*(spD'*e1);
    Dt2_Ks=[Dt2_k11,Dt2_k22,0.5*(Dt2_k12+Dt2_k21)];

    Dw1t1_k11=(E11'*Dw1RtdR*E33)*((dsp'*e1)*(spD'*e1)');
    Dw1t1_k22=(E11'*Dw1RtdR*E33)*((dsp'*e2)*(spD'*e2)');
    Dw1t1_k12=(E11'*Dw1RtdR*E33)*((dsp'*e1)*(spD'*e2)');
    Dw1t1_k21=(E11'*Dw1RtdR*E33)*((dsp'*e2)*(spD'*e1)');
    Dw1t1_k12s=0.5*(Dw1t1_k12+Dw1t1_k21);

    Dw1t2_k11=(E22'*Dw1RtdR*E33)*((dsp'*e1)*(spD'*e1)');
    Dw1t2_k22=(E22'*Dw1RtdR*E33)*((dsp'*e2)*(spD'*e2)');
    Dw1t2_k12=(E22'*Dw1RtdR*E33)*((dsp'*e1)*(spD'*e2)');
    Dw1t2_k21=(E22'*Dw1RtdR*E33)*((dsp'*e2)*(spD'*e1)');
    Dw1t2_k12s=0.5*(Dw1t2_k12+Dw1t2_k21);

    Dw2t1_k11=(E11'*Dw2RtdR*E33)*((dsp'*e1)*(spD'*e1)');
    Dw2t1_k22=(E11'*Dw2RtdR*E33)*((dsp'*e2)*(spD'*e2)');
    Dw2t1_k12=(E11'*Dw2RtdR*E33)*((dsp'*e1)*(spD'*e2)');
    Dw2t1_k21=(E11'*Dw2RtdR*E33)*((dsp'*e2)*(spD'*e1)');
    Dw2t1_k12s=0.5*(Dw2t1_k12+Dw2t1_k21);

    Dw2t2_k11=(E22'*Dw2RtdR*E33)*((dsp'*e1)*(spD'*e1)');
    Dw2t2_k22=(E22'*Dw2RtdR*E33)*((dsp'*e2)*(spD'*e2)');
    Dw2t2_k12=(E22'*Dw2RtdR*E33)*((dsp'*e1)*(spD'*e2)');
    Dw2t2_k21=(E22'*Dw2RtdR*E33)*((dsp'*e2)*(spD'*e1)');
    Dw2t2_k12s=0.5*(Dw2t2_k12+Dw2t2_k21);

    %Metric
    a11=(t1'*e1)*(t1'*e1)+(t2'*e1)*(t2'*e1);
    a22=(t1'*e2)*(t1'*e2)+(t2'*e2)*(t2'*e2);
    a12=(t1'*e1)*(t1'*e2)+(t2'*e1)*(t2'*e2);
    a21=(t2'*e1)*(t2'*e2)+(t1'*e1)*(t1'*e2);


    Ms=[a11-Me0(1,1),a22-Me0(2,2),(a12-Me0(1,2))];
    % derivatives of metric
    Dt1_a11=2*(t1'*e1)*(spD'*e1);
    Dt1_a22=2*(t1'*e2)*(spD'*e2);
    Dt1_a12=(t1'*e1)*(spD'*e2)+(t1'*e2)*(spD'*e1);
    Dt1_Ms=[Dt1_a11,Dt1_a22,Dt1_a12];

    Dt2_a11=2*(t2'*e1)*(spD'*e1);
    Dt2_a22=2*(t2'*e2)*(spD'*e2);
    Dt2_a12=(t2'*e1)*(spD'*e2)+(t2'*e2)*(spD'*e1);
    Dt2_Ms=[Dt2_a11,Dt2_a22,(Dt2_a12)];

    Dt1t1_a11=2*(spD'*e1)*(spD'*e1)';
    Dt2t2_a11=2*(spD'*e1)*(spD'*e1)';
    Dt1t2_a11=zeros(3,3);

    Dt1t1_a22=2*(spD'*e2)*(spD'*e2)';
    Dt2t2_a22=2*(spD'*e2)*(spD'*e2)';
    Dt1t2_a22=zeros(3,3);

    Dt1t1_a12=((spD'*e1)*(spD'*e2)'+(spD'*e2)*(spD'*e1)');
    Dt2t2_a12=((spD'*e1)*(spD'*e2)'+(spD'*e2)*(spD'*e1)');
    Dt1t2_a12=zeros(3,3);

    %***************************Residual Vector ***************************
    f(r1Dof)=f(r1Dof)+JxW*(Cb*(Ks*Hs*Dw1_Ks')'+(p1'*Dw1R'*t1u')*sp+(p2'*Dw1R'*t2u')*sp);
    f(r2Dof)=f(r2Dof)+JxW*(Cb*(Ks*Hs*Dw2_Ks')'+(p1'*Dw2R'*t1u')*sp+(p2'*Dw2R'*t2u')*sp);

    f(t1Dof)=f(t1Dof)+JxW*(Cm*(Ms*Hs*Dt1_Ms')'+Cb*(Ks*Hs*Dt1_Ks')'...
        -(p1(2)-p2(1))*spDt2+(p1'*R'*spD_u)');
    f(t2Dof)=f(t2Dof)+JxW*(Cm*(Ms*Hs*Dt2_Ms')'+Cb*(Ks*Hs*Dt2_Ks')'...
        +(p1(2)-p2(1))*spDt1+(p2'*R'*spD_u)');

    f(p1Dof(1))=f(p1Dof(1))+JxW*(E11'*R'*t1u');
    f(p1Dof(2))=f(p1Dof(2))+JxW*(-t1_t2+E22'*R'*t1u');
    f(p1Dof(3))=f(p1Dof(3))+JxW*(E33'*R'*t1u');

    f(p2Dof(1))=f(p2Dof(1))+JxW*(t1_t2+E11'*R'*t2u');
    f(p2Dof(2))=f(p2Dof(2))+JxW*(E22'*R'*t2u');
    f(p2Dof(3))=f(p2Dof(3))+JxW*(E33'*R'*t2u');

    f(u1Dof)=f(u1Dof)+JxW*(-(p1'*R'*E11)*dspt1-(p2'*R'*E11)*dspt2+sp*q(1));
    f(u2Dof)=f(u2Dof)+JxW*(-(p1'*R'*E22)*dspt1-(p2'*R'*E22)*dspt2+sp*q(2));
    f(u3Dof)=f(u3Dof)+JxW*(-(p1'*R'*E33)*dspt1-(p2'*R'*E33)*dspt2+sp*q(3));

    %***************************Tangent operator **************************
    %******************* (rotation, rotation) terms ***********************
    tmp=((Ks*Hs*E11)*Dw1w1_k11+Ks*Hs*E22*Dw1w1_k22+Ks*Hs*E33*Dw1w1_k12s);
    k(r1Dof,r1Dof)=k(r1Dof,r1Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dw1_Ks')+Cb*tmp+...
        (p1'*Dw1w1R'*t1u')*(sp*sp')+(p2'*Dw1w1R'*t2u')*(sp*sp'));

    tmp=(Ks*Hs*E11*Dw2w2_k11+Ks*Hs*E22*Dw2w2_k22+Ks*Hs*E33*Dw2w2_k12s);
    k(r2Dof,r2Dof)=k(r2Dof,r2Dof)+JxW*(Cb*(Dw2_Ks*Hs*Dw2_Ks')+Cb*tmp+...
        (p1'*Dw2w2R'*t1u')*(sp*sp')+(p2'*Dw2w2R'*t2u')*(sp*sp'));

    tmp=(Ks*Hs*E11*Dw1w2_k11+Ks*Hs*E22*Dw1w2_k22+Ks*Hs*E33*Dw1w2_k12s);
    k(r1Dof,r2Dof)=k(r1Dof,r2Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dw2_Ks')+Cb*tmp+...
        (p1'*Dw1w2R'*t1u')*(sp*sp')+(p2'*Dw1w2R'*t2u')*(sp*sp'));

    k(r2Dof,r1Dof)=k(r2Dof,r1Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dw2_Ks')+Cb*tmp+...
        (p1'*Dw1w2R'*t1u')*(sp*sp')+(p2'*Dw1w2R'*t2u')*(sp*sp'))';

    %************* (theta, theta) terms *****************

    tmp=Cb*(Dt1_Ks*Hs*Dt1_Ks');
    k(t1Dof,t1Dof)=k(t1Dof,t1Dof)+JxW*(tmp+Cm*(Dt1_Ms*Hs*Dt1_Ms'+...
        Ms*Hs*E11*Dt1t1_a11+Ms*Hs*E22*Dt1t1_a22+Ms*Hs*E33*Dt1t1_a12));

    tmp=Cb*(Dt2_Ks*Hs*Dt2_Ks');
    k(t2Dof,t2Dof)=k(t2Dof,t2Dof)+JxW*(tmp+Cm*(Dt2_Ms*Hs*Dt2_Ms'+...
        Ms*Hs*E11*Dt2t2_a11+Ms*Hs*E22*Dt2t2_a22+Ms*Hs*E33*Dt2t2_a12));

    tmp=Cb*(Dt1_Ks*Hs*Dt2_Ks')-(p1(2)-p2(1))*spD_spD;
    k(t1Dof,t2Dof)=k(t1Dof,t2Dof)+JxW*(tmp+Cm*(Dt1_Ms*Hs*Dt2_Ms'+...
        Ms*Hs*E11*Dt1t2_a11+Ms*Hs*E22*Dt1t2_a22+Ms*Hs*E33*Dt1t2_a12));

    k(t2Dof,t1Dof)=k(t2Dof,t1Dof)+JxW*(tmp+Cm*(Dt1_Ms*Hs*Dt2_Ms'+...
        Ms*Hs*E11*Dt1t2_a11+Ms*Hs*E22*Dt1t2_a22+Ms*Hs*E33*Dt1t2_a12))';

    %******************* (rotation, theta) terms ***********************
    tmp=((Ks*Hs*E11)*Dw1t1_k11+(Ks*Hs*E22)*Dw1t1_k22+(Ks*Hs*E33)*Dw1t1_k12s);
    k(r1Dof,t1Dof)=k(r1Dof,t1Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dt1_Ks'+tmp)+(sp*(p1'*Dw1R'*spD_u)));
    k(t1Dof,r1Dof)=k(t1Dof,r1Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dt1_Ks'+tmp)+(sp*(p1'*Dw1R'*spD_u)))';

    tmp=((Ks*Hs*E11)*Dw2t2_k11+(Ks*Hs*E22)*Dw2t2_k22+(Ks*Hs*E33)*Dw2t2_k12s);
    k(r2Dof,t2Dof)=k(r2Dof,t2Dof)+JxW*(Cb*(Dw2_Ks*Hs*Dt2_Ks'+tmp)+(sp*(p2'*Dw2R'*spD_u)));
    k(t2Dof,r2Dof)=k(t2Dof,r2Dof)+JxW*(Cb*(Dw2_Ks*Hs*Dt2_Ks'+tmp)+(sp*(p2'*Dw2R'*spD_u)))';

    tmp=((Ks*Hs*E11)*Dw1t2_k11+(Ks*Hs*E22)*Dw1t2_k22+(Ks*Hs*E33)*Dw1t2_k12s);
    k(r1Dof,t2Dof)=k(r1Dof,t2Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dt2_Ks'+tmp)+(sp*(p2'*Dw1R'*spD_u)));
    k(t2Dof,r1Dof)=k(t2Dof,r1Dof)+JxW*(Cb*(Dw1_Ks*Hs*Dt2_Ks'+tmp)+(sp*(p2'*Dw1R'*spD_u)))';

    tmp=((Ks*Hs*E11)*Dw2t1_k11+(Ks*Hs*E22)*Dw2t1_k22+(Ks*Hs*E33)*Dw2t1_k12s);
    k(r2Dof,t1Dof)=k(r2Dof,t1Dof)+JxW*(Cb*(Dw2_Ks*Hs*Dt1_Ks'+tmp)+(sp*(p1'*Dw2R'*spD_u)));
    k(t1Dof,r2Dof)=k(t1Dof,r2Dof)+JxW*(Cb*(Dw2_Ks*Hs*Dt1_Ks'+tmp)+(sp*(p1'*Dw2R'*spD_u)))';

    %******************* (rotation, displacement) terms ***********************

    k(r1Dof,u1Dof)=k(r1Dof,u1Dof)+JxW*-((p1'*Dw1R'*E11)*sp*dspt1'+(p2'*Dw1R'*E11)*sp*dspt2');
    k(u1Dof,r1Dof)=k(u1Dof,r1Dof)+JxW*-((p1'*Dw1R'*E11)*sp*dspt1'+(p2'*Dw1R'*E11)*sp*dspt2')';

    k(r1Dof,u2Dof)=k(r1Dof,u2Dof)+JxW*-((p1'*Dw1R'*E22)*sp*dspt1'+(p2'*Dw1R'*E22)*sp*dspt2');
    k(u2Dof,r1Dof)=k(u2Dof,r1Dof)+JxW*-((p1'*Dw1R'*E22)*sp*dspt1'+(p2'*Dw1R'*E22)*sp*dspt2')';

    k(r1Dof,u3Dof)=k(r1Dof,u3Dof)+JxW*-((p1'*Dw1R'*E33)*sp*dspt1'+(p2'*Dw1R'*E33)*sp*dspt2');
    k(u3Dof,r1Dof)=k(u3Dof,r1Dof)+JxW*-((p1'*Dw1R'*E33)*sp*dspt1'+(p2'*Dw1R'*E33)*sp*dspt2')';

    k(r2Dof,u1Dof)=k(r2Dof,u1Dof)+JxW*-((p1'*Dw2R'*E11)*sp*dspt1'+(p2'*Dw2R'*E11)*sp*dspt2');
    k(u1Dof,r2Dof)=k(u1Dof,r2Dof)+JxW*-((p1'*Dw2R'*E11)*sp*dspt1'+(p2'*Dw2R'*E11)*sp*dspt2')';

    k(r2Dof,u2Dof)=k(r2Dof,u2Dof)+JxW*-((p1'*Dw2R'*E22)*sp*dspt1'+(p2'*Dw2R'*E22)*sp*dspt2');
    k(u2Dof,r2Dof)=k(u2Dof,r2Dof)+JxW*-((p1'*Dw2R'*E22)*sp*dspt1'+(p2'*Dw2R'*E22)*sp*dspt2')';

    k(r2Dof,u3Dof)=k(r2Dof,u3Dof)+JxW*-((p1'*Dw2R'*E33)*sp*dspt1'+(p2'*Dw2R'*E33)*sp*dspt2');
    k(u3Dof,r2Dof)=k(u3Dof,r2Dof)+JxW*-((p1'*Dw2R'*E33)*sp*dspt1'+(p2'*Dw2R'*E33)*sp*dspt2')';

    %******************* (theta, displacement) terms ***********************
    k(t1Dof,u1Dof)=k(t1Dof,u1Dof)+JxW*(p1'*R'*E11)*spD_dsp;
    k(u1Dof,t1Dof)=k(u1Dof,t1Dof)+JxW*((p1'*R'*E11)*spD_dsp)';

    k(t1Dof,u2Dof)=k(t1Dof,u2Dof)+JxW*(p1'*R'*E22)*spD_dsp;
    k(u2Dof,t1Dof)=k(u2Dof,t1Dof)+JxW*((p1'*R'*E22)*spD_dsp)';

    k(t1Dof,u3Dof)=k(t1Dof,u3Dof)+JxW*(p1'*R'*E33)*spD_dsp;
    k(u3Dof,t1Dof)=k(u3Dof,t1Dof)+JxW*((p1'*R'*E33)*spD_dsp)';

    k(t2Dof,u1Dof)=k(t2Dof,u1Dof)+JxW*(p2'*R'*E11)*spD_dsp;
    k(u1Dof,t2Dof)=k(u1Dof,t2Dof)+JxW*((p2'*R'*E11)*spD_dsp)';

    k(t2Dof,u2Dof)=k(t2Dof,u2Dof)+JxW*(p2'*R'*E22)*spD_dsp;
    k(u2Dof,t2Dof)=k(u2Dof,t2Dof)+JxW*((p2'*R'*E22)*spD_dsp)';

    k(t2Dof,u3Dof)=k(t2Dof,u3Dof)+JxW*(p2'*R'*E33)*spD_dsp;
    k(u3Dof,t2Dof)=k(u3Dof,t2Dof)+JxW*((p2'*R'*E33)*spD_dsp)';

    %*************** (rotation, tractin)terms ***************************
    r1p1=zeros(3,3);
    r1p1(1:3,1)=((t1u*Dw1R*E11)*sp);                     %(r1Dof,p1Dof)
    r1p1(1:3,2)=((t1u*Dw1R*E22)*sp);
    r1p1(1:3,3)=((t1u*Dw1R*E33)*sp);
    k(r1Dof,p1Dof)=k(r1Dof,p1Dof)+JxW*r1p1;
    k(p1Dof,r1Dof)=k(p1Dof,r1Dof)+JxW*r1p1';

    r1p2=zeros(3,3);
    r1p2(1:3,1)=((t2u*Dw1R*E11)*sp);                     %(r1Dof,p2Dof)
    r1p2(1:3,2)=((t2u*Dw1R*E22)*sp);
    r1p2(1:3,3)=((t2u*Dw1R*E33)*sp);
    k(r1Dof,p2Dof)=k(r1Dof,p2Dof)+JxW*r1p2;
    k(p2Dof,r1Dof)=k(p2Dof,r1Dof)+JxW*r1p2';

    r2p1=zeros(3,3);
    r2p1(1:3,1)=((t1u*Dw2R*E11)*sp);                     %(r2Dof,p1Dof)
    r2p1(1:3,2)=((t1u*Dw2R*E22)*sp);
    r2p1(1:3,3)=((t1u*Dw2R*E33)*sp);
    k(r2Dof,p1Dof)=k(r2Dof,p1Dof)+JxW*r2p1;
    k(p1Dof,r2Dof)=k(p1Dof,r2Dof)+JxW*r2p1';

    r2p2=zeros(3,3);
    r2p2(1:3,1)=((t2u*Dw2R*E11)*sp);                     %(r2Dof,p2Dof)
    r2p2(1:3,2)=((t2u*Dw2R*E22)*sp);
    r2p2(1:3,3)=((t2u*Dw2R*E33)*sp);
    k(r2Dof,p2Dof)=k(r2Dof,p2Dof)+JxW*r2p2;
    k(p2Dof,r2Dof)=k(p2Dof,r2Dof)+JxW*r2p2';

    %***************(displacement, traction) terms ************************
    u1p1=zeros(3,3);
    u1p1(1:3,1)=-((E11'*R*E11)*dspt1) ;                  %(u1Dof,p1Dof)
    u1p1(1:3,2)=-((E11'*R*E22)*dspt1) ;
    u1p1(1:3,3)=-((E11'*R*E33)*dspt1) ;
    k(u1Dof,p1Dof)=k(u1Dof,p1Dof)+JxW*u1p1;
    k(p1Dof,u1Dof)=k(p1Dof,u1Dof)+JxW*u1p1';

    u1p2=zeros(3,3);
    u1p2(1:3,1)=-((E11'*R*E11)*dspt2);                   %(u1Dof,p2Dof)
    u1p2(1:3,2)=-((E11'*R*E22)*dspt2);
    u1p2(1:3,3)=-((E11'*R*E33)*dspt2);
    k(u1Dof,p2Dof)=k(u1Dof,p2Dof)+JxW*u1p2;
    k(p2Dof,u1Dof)=k(p2Dof,u1Dof)+JxW*u1p2';

    u2p1=zeros(3,3);                                     %(u2Dof,p1Dof)
    u2p1(1:3,1)=-((E22'*R*E11)*dspt1);
    u2p1(1:3,2)=-((E22'*R*E22)*dspt1) ;
    u2p1(1:3,3)=-((E22'*R*E33)*dspt1) ;
    k(u2Dof,p1Dof)=k(u2Dof,p1Dof)+JxW*u2p1;
    k(p1Dof,u2Dof)=k(p1Dof,u2Dof)+JxW*u2p1';

    u2p2=zeros(3,3);
    u2p2(1:3,1)=-((E22'*R*E11)*dspt2);                  %(u2Dof,p2Dof)
    u2p2(1:3,2)=-((E22'*R*E22)*dspt2);
    u2p2(1:3,3)=-((E22'*R*E33)*dspt2);
    k(u2Dof,p2Dof)=k(u2Dof,p2Dof)+JxW*u2p2;
    k(p2Dof,u2Dof)=k(p2Dof,u2Dof)+JxW*u2p2';

    u3p1=zeros(3,3);                                     %(u3Dof,p1Dof)
    u3p1(1:3,1)=-((E33'*R*E11)*dspt1);
    u3p1(1:3,2)=-((E33'*R*E22)*dspt1);
    u3p1(1:3,3)=-((E33'*R*E33)*dspt1);
    k(u3Dof,p1Dof)=k(u3Dof,p1Dof)+JxW*u3p1;
    k(p1Dof,u3Dof)=k(p1Dof,u3Dof)+JxW*u3p1';

    u3p2=zeros(3,3);                                     %(u3Dof,p2Dof)
    u3p2(1:3,1)=-((E33'*R*E11)*dspt2);
    u3p2(1:3,2)=-((E33'*R*E22)*dspt2);
    u3p2(1:3,3)=-((E33'*R*E33)*dspt2);
    k(u3Dof,p2Dof)=k(u3Dof,p2Dof)+JxW*u3p2;
    k(p2Dof,u3Dof)=k(p2Dof,u3Dof)+JxW*u3p2';

    %****************( traction, theta) terms *****************************

    p1t1=zeros(3,3);
    p1t1(1,1:3)=(E11'*R'*spD_u);                    %(p1Dof[0],t1Dof)
    p1t1(2,1:3)=-spDt2'+(E22'*R'*spD_u);             %(p1Dof[1],t1Dof)
    p1t1(3,1:3)=(E33'*R'*spD_u);                    %(p1Dof[2],t1Dof)
    k(p1Dof,t1Dof)=k(p1Dof,t1Dof)+JxW*p1t1;
    k(t1Dof,p1Dof)=k(t1Dof,p1Dof)+JxW*p1t1';

    p1t2=zeros(3,3);
    p1t2(2,1:3)=spDt1';                                %(p1Dof[1],t2Dof)
    k(p1Dof,t2Dof)=k(p1Dof,t2Dof)+JxW*p1t2;
    k(t2Dof,p1Dof)=k(t2Dof,p1Dof)+JxW*p1t2';

    p2t1=zeros(3,3);
    p2t1(1,1:3)=spDt2';                                %(p2Dof[0],t1Dof)
    k(p2Dof,t1Dof)=k(p2Dof,t1Dof)+JxW*p2t1;
    k(t1Dof,p2Dof)=k(t1Dof,p2Dof)+JxW*p2t1';

    p2t2=zeros(3,3);
    p2t2(1,1:3)=-spDt1'+(E11'*R'*spD_u) ;            %(p2Dof[0],t2Dof)
    p2t2(2,1:3)=(E22'*R'*spD_u);                    %(p2Dof[1],t2Dof)
    p2t2(3,1:3)=(E33'*R'*spD_u) ;                   %(p2Dof[2],t2Dof)
    k(p2Dof,t2Dof)=k(p2Dof,t2Dof)+JxW*p2t2;
    k(t2Dof,p2Dof)=k(t2Dof,p2Dof)+JxW*p2t2';

end
