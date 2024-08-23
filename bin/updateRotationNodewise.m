%Nodal rotation update is based on geometrically exact nodal director
%update procedure given in Simo part 3, Box 2.

function rotation=updateRotationNodewise(me,rotation,incRot)
for itEl=1:me.noNd
    dT=incRot(itEl,:);               %incremental reference rotations
    R=rotation{itEl} ;               %Rotation matrix at current confign
    
    dt=R(:,1:2)*dT';                  %spatial director increments
    t=R(:,3);                        %director at given state
    dtheta1=cross(t,dt);
    W=axial2Skew(dtheta1) ;
    
    %%Method 1
    rotation{itEl}=expm(W)*R;             %updated rotation
    %%Method 2
%     a1=norm(dt);
%     if a1==0
%         dR=eye(3)+W;
%     else
%         dR=cos(a1)*eye(3)+(sin(a1)/a1)*W+((1-cos(a1))/(a1*a1))*(dtheta1*dtheta1');
%         %          dR=expm(W);
%     end
%     rotation{itEl}=dR*R;               %updated rotation

     %%Checks
    
    % %     %check spatial director increment
    % %     %Checked , t1=t2 and t'*dt=0
    % %     t1=cos(a1)*t+(sin(a1)/a1)*dt;
    % %     t2=R(:,3);
    
    
end