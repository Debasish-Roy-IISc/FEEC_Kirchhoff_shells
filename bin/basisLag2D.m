function fe = basisLag2D(p)

fe.noGp = size(p,1);
fe.sp  = zeros(3,fe.noGp);
fe.dsp = zeros(2, 3, fe.noGp);
for i=1:fe.noGp

    % calculate basis function values
    fe.sp(:,i) = [1-p(i,1)-p(i,2); p(i,1); p(i,2)];
    % Derivative of shape function
    fe.dsp(:,:,i) = [-1 1 0 ;
                     -1 0 1 ];

%     % calculate basis function values
%     fe.sp(:,i) = [p(i,1); p(i,2); 1-p(i,1)-p(i,2)];
%     % Derivative of shape function
%     fe.dsp(:,:,i) = [1 0 -1;0 1 -1];
end