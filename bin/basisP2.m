function fe = basisP2(p)

fe.noGp = size(p,1);
fe.phi  = zeros(12, 2, fe.noGp);
for i=1:fe.noGp
    l1=p(i,1);
    l2=p(i,2);
    l3=1-l1-l2;
    dl1=[1 0];
    dl2=[0 1];
    dl3=-[1 1];
%     W12=l1*dl2-l2*dl1;
%     W23=l2*dl3-l3*dl2;
%     W31=l3*dl1-l1*dl3;
    % calculate basis function values
    fe.phi(1,:,i) = l1^2*dl2;
    fe.phi(2,:,i) = l2^2*dl1;
    fe.phi(3,:,i) = l1*l2*(dl2-dl1);
    
    fe.phi(4,:,i) = l2^2*dl3;
    fe.phi(5,:,i) = l3^2*dl2;
    fe.phi(6,:,i) = l2*l3*(dl3-dl2);
    
    fe.phi(7,:,i) = l3^2*dl1;
    fe.phi(8,:,i) = l1^2*dl3;
    fe.phi(9,:,i) = l1*l3*(dl1-dl3);
    
    fe.phi(10,:,i) = dl1*l2*l3;
    fe.phi(11,:,i) = dl2*l3*l1;
    fe.phi(12,:,i) = dl3*l1*l2;
    
%     fe.phi(10,:,i) = l1*l2*dl3;
%     fe.phi(11,:,i) = l2*l3*dl1;
%     fe.phi(12,:,i) = l3*l1*dl2;
end