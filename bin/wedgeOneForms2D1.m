%author: Jamun
%Input:Two one forms a and b
%Output: Wedge product of one forms
function we = wedgeOneForms2D1(a,b)
we=zeros(size(a,1),size(b,1),1);
for j=1:size(a,1)
    for k=1:size(b,1)
        we(j,k,:)=det([a(j,:)' b(k,:)']);
    end
end
