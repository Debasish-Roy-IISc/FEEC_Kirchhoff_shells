function fe = basisP2Minus(p)

fe.noGp = size(p,1);
fe.phi  = zeros(8, 2, fe.noGp);
% M=eye(2);
M=-[0 1;-1 0];


for i=1:fe.noGp
    l1=p(i,1);
    l2=p(i,2);
    l3=1-l1-l2;
    dl1=[1 0];
    dl2=[0 1];
    dl3=-[1 1];
    W12=l1*dl2-l2*dl1;
    W23=l2*dl3-l3*dl2;
    W31=l3*dl1-l1*dl3;
    % calculate basis function values
    fe.phi(1,:,i) = l1*W12*M;
    fe.phi(2,:,i) = l2*W12*M;
    
    fe.phi(3,:,i) = l2*W23*M;
    fe.phi(4,:,i) = l3*W23*M;
    
    fe.phi(5,:,i) = l3*W31*M;
    fe.phi(6,:,i) = l1*W31*M;
    
    fe.phi(7,:,i) = l3*W12*M;
    fe.phi(8,:,i) = -l2*W31*M;
end