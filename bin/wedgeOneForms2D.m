
%Input:Two one forms a and b
%Output: Wedge product of one forms
function we = wedgeOneForms2D(a,b)
% we=zeros(size(a,1),size(b,1),3);
we=zeros(size(a,1),1);
for j=1:size(a,1)
%     for k=1:size(b,1)
%         we(j,k,:)=cross(a(j,:),b(k,:));
        we(j,:)=det([(a(j,:))' b]);
%     end
end
