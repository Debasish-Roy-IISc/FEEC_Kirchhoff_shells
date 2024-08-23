function R=rodriguesFormula(W)
theta=(0.5*trace(W'*W))^0.5;
R=eye(3);
if theta>1e-15
    R=eye(3)+(sin(theta)*W/theta)+((1-cos(theta))*W^2/theta^2);
end