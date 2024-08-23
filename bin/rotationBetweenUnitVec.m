function R=rotationBetweenUnitVec(t1,t2)
% Need to check for corner cases
s=cross(t1,t2);
axS=axial2Skew(s);
c=t1*t2';
s=norm(s);
if c==1
    R=eye(3)+axS;
else
R=eye(3)+axS+((axS*axS)*((1-c)/s^2));
end

