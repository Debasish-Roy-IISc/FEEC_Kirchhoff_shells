function fe = basisP1Minus(p)

fe.noGp = size(p,1);
fe.phi  = zeros(3, 2, fe.noGp);
fe.dphi = zeros(3, 1, fe.noGp);
% M=-[0 1;-1 0];
M=eye(2);
for i=1:fe.noGp
    l1=p(i,1);
    l2=p(i,2);
    l3=1-l1-l2;
    dl1=[1 0];
    dl2=[0 1];
    dl3=-[1 1];
    
    % calculate basis function values
    %     fe.phi(1,:,i) = [ -p(i,2) p(i,1)];
    %     fe.phi(2,:,i) = [ -p(i,2) p(i,1)-1];
    %     fe.phi(3,:,i) = [ -p(i,2)+1 p(i,1)];
    fe.phi(1,:,i) = (l1*dl2-l2*dl1)*M;
    fe.phi(2,:,i) = (l2*dl3-l3*dl2)*M;
    fe.phi(3,:,i) = (l3*dl1-l1*dl3)*M;
    
%     % calculate basis function values
%     fe.phi(1,:,i) = [ p(i,1)     p(i,2)   ];
%     fe.phi(2,:,i) = [ p(i,1)-1   p(i,2)   ];
%     fe.phi(3,:,i) = [ p(i,1)     p(i,2)-1 ];
%     
%     % calculate basis function divergence values
%     fe.dphi(:,:,i) = [sqrt(2); 2; 2];
end