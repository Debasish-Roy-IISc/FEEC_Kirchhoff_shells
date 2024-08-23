function fe = basis2D1P1(p)

fe.noGp = size(p,1);
fe.phi  = zeros(6, 2, fe.noGp);
for i=1:fe.noGp
    l2=p(i,1);
    l3=p(i,2);
    l1=1-l2-l3;
    
    dl1=-[1 1];
    dl2=[1 0];
    dl3=[0 1];

    % calculate basis function values
    fe.phi(1,:,i) = l1*dl2;      %lexicographic order
    fe.phi(4,:,i) = l2*dl1;
    
    fe.phi(2,:,i) = l1*dl3;
    fe.phi(5,:,i) = l3*dl1;
    
    fe.phi(3,:,i) = l2*dl3;
    fe.phi(6,:,i) = l3*dl2;

%     fe.phi(1,:,i) = l2*dl3;
%     fe.phi(4,:,i) = l3*dl2;
%     
%     fe.phi(2,:,i) = l1*dl3;
%     fe.phi(5,:,i) = l3*dl1;
%     
%     fe.phi(3,:,i) = l1*dl2;
%     fe.phi(6,:,i) = l2*dl1;
end