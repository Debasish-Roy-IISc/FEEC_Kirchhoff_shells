function W=axial2Skew(v)
% Returns the axial vector associated with a skew symmetric matrix W
W=zeros(3);
W(1,2)=-v(3);
W(1,3)=v(2);
W(2,3)=-v(1);
W=W-W';
end
