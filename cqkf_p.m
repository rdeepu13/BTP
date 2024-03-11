%Points and weights generation in CQKF

%vpa(Weight), Here I have used vpa to approximate symbolic expression,
 %here "double" won't work; vpa(Weight)

 function [CQ_points,Weight] = cqkf_p(n,n1)
v=sym('v');
% num_dim=2;
% num_quadrature1=2;
m=1; %1 2 3
%m=1 Third degree of cubature rule, m=2 Fifth degree of cubature rule, 
%m=3 Seventh degree of cubature rule 

%n=num_dim;   %Dimension of the system
alph=(n/2)-1;  %alpha 
%n1=num_quadrature1;   %Number of roots
L(1,1)=1;
mult=-n1*(n1+alph);
L(1,2)=mult;
for i=1:n1-1
    mult=(-1)*mult*(n1-i)*(n1+alph-i);
    L(1,i+2)=mult/factorial(i+1);
end
o=L(1,:);  %Coefficient of chebyhev-Laguerre polynomial
r=roots(o);  %
pv=(v^(alph+n1))*exp(-v);
for i=1:n1
    L1=diff(pv, sym('v'));
    pv=L1;
end
L2=((-1)^n1)*(v^(-alph))*(exp(v))*pv;
L3=diff(L2);
points=subs(L3, 'v', r);  %Differentiated Chebyshev-Laguerre polynomial


y2=[eye(n) -eye(n)];
for i=1:n
    u(:,2*i)=y2(:,n+i);
    u(:,2*i-1)=y2(:,i);
end
U=u;

for i=1:n1
    cqkf(:,:,i)=sqrt(2*r(i,1))*U;
end
cqkf;

cqkf1=[];
for i=1:n1
    cqkf1=[cqkf1 cqkf(:,:,i)];
end
cqkf1;


for k=1:n1
    Weight1(k)=(1/(2*n*gamma(n/2)))*((factorial(n1)*gamma(alph+n1+1))/(r(k,1)*points(k,1)*points(k,1)));
end
for k=1:n1
for i=1:2*n
    Weight2(k,i)=Weight1(k);
end
end
Weight=[];
    for k=1:n1
        Weight=[Weight Weight2(k,:)];
    end
    Weight;
    CQ_points=cqkf1;
    
    Weight=double(Weight);
    %vpa(Weight)
 end
    