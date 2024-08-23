%Nodal rotation update 

function rotation=updateRotationNodewise_test(me,rotation,incRot)
for itEl=1:me.noNd
    dT=incRot(itEl,:);               %incremental reference rotations
    R=rotation{itEl} ;               %Rotation matrix at current confign
    dtheta1=[dT,0]';
    W=axial2Skew(dtheta1) ;
    
    %%Method 1
    rotation{itEl}=expm(W)*R;             %updated rotation
  
end
