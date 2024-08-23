function rotation=updateRotationElementwise(me,rotation,incRot)
for itEl=1:me.noEl
    %***************method 1 : Given by simo part 3 box 2, page 15
        elConn=me.elements(itEl,:);
        I=eye(3);
        %incremental reference rotations
        dT1=incRot(elConn(1),:); dT2=incRot(elConn(2),:); dT3=incRot(elConn(3),:);
        R1=rotation{elConn(1)} ; R2=rotation{elConn(2)} ; R3=rotation{elConn(3)} ;
    
        %incremental spatial rotations
        r1=R1(:,1:2)*dT1';  r2=R2(:,1:2)*dT2'; r3=R3(:,1:2)*dT3';
        t1=R1(:,3);         t2=R2(:,3);        t3=R3(:,3);
        dtheta1=cross(t1,r1);dtheta2=cross(t2,r2);dtheta3=cross(t3,r3);
    
        %incremental rotation matrix for each node
        a1=norm(r1);  a2=norm(r2);   a3=norm(r3);
        W1=axial2Skew(dtheta1) ; W2=axial2Skew(dtheta2) ; W3=axial2Skew(dtheta3) ;
        
        if a1==0
            dR1=I+W1;
        else
        dR1=cos(a1)*I+(sin(a1)/a1)*W1+((1-cos(a1))/(a1*a1))*(dtheta1*dtheta1');
        end
         if a2==0
            dR2=I+W2;
        else
        dR2=cos(a2)*I+(sin(a2)/a2)*W2+((1-cos(a2))/(a2*a2))*(dtheta2*dtheta2');
         end
         if a3==0
            dR3=I+W3;
        else
        dR3=cos(a3)*I+(sin(a3)/a3)*W3+((1-cos(a3))/(a3*a3))*(dtheta3*dtheta3');
         end

        rotation{elConn(1)}=dR1*R1;
        rotation{elConn(2)}=dR2*R2;
        rotation{elConn(3)}=dR3*R3;
    
    %****************************method 2******************************
%     elConn=me.elements(itEl,:);
%     dr1=[incRot(elConn(1),:),0]; dr2=[incRot(elConn(2),:),0]; dr3=[incRot(elConn(3),:),0];
%     R1=rotation{elConn(1)} ; R2=rotation{elConn(2)} ; R3=rotation{elConn(3)} ;
%     W1=axial2Skew(dr1) ; W2=axial2Skew(dr2) ; W3=axial2Skew(dr3) ;
%     dR1=expm(W1);dR2=expm(W2);dR3=expm(W3);
%     rotation{elConn(1)}=dR1*R1;
%     rotation{elConn(2)}=dR2*R2;
%     rotation{elConn(3)}=dR3*R3;
% %     rotation{elConn(1)}=R1*dR1;
% %     rotation{elConn(2)}=R2*dR2;
% %     rotation{elConn(3)}=R3*dR3;
end