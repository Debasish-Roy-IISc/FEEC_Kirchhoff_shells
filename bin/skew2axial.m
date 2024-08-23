function w=skew2axial(W)
w=zeros(3,1);
w(1)=-W(2,3);
w(2)=W(1,3);
w(3)=-W(1,2);