function X0 = S2SInMin( kmax )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    T=1;

xo(:,1)=[0;0;0.15433*sin(140*pi/180);0.15433*cos(140*pi/180)];
for k=1:12
    xo(1,k+1)=xo(1,k)+0.15433*sin(140*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(140*pi/180)*T;
    xo(3,k+1)=xo(3,k);
    xo(4,k+1)=xo(4,k);
end

for k=13
    xo(1,k+1)=xo(1,k)+0.15433*sin(125*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(125*pi/180)*T;
    xo(3,k+1)=0.15433*sin(125*pi/180);
    xo(4,k+1)=0.15433*cos(125*pi/180);
end 

for k=14
    xo(1,k+1)=xo(1,k)+0.15433*sin(96*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(96*pi/180)*T;
    xo(3,k+1)=0.15433*sin(96*pi/180);
    xo(4,k+1)=0.15433*cos(96*pi/180);
end

for k=15
    xo(1,k+1)=xo(1,k)+0.15433*sin(65*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(65*pi/180)*T;
    xo(3,k+1)=0.15433*sin(65*pi/180);
    xo(4,k+1)=0.15433*cos(65*pi/180);
end

for k=16
    xo(1,k+1)=xo(1,k)+0.15433*sin(34*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(34*pi/180)*T;
    xo(3,k+1)=0.15433*sin(34*pi/180);
    xo(4,k+1)=0.15433*cos(34*pi/180);
end

for k=17:kmax
    xo(1,k+1)=xo(1,k)+0.15433*sin(20*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(20*pi/180)*T;
    xo(3,k+1)=0.15433*sin(20*pi/180);
    xo(4,k+1)=0.15433*cos(20*pi/180);
end

    
    X0=xo;


end

