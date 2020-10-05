function [ mx,mn ] = findmxmn( N,Map2DSave )
m1=zeros(N,1);
m2=zeros(N,1);
for i=1:N
    A=cell2mat(Map2DSave(i));
    m1(i)=max(max(A));
m2(i)=min(min(A));

end
mx=max(m1);
mn=min(m2);
end