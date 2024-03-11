function X0 = S2SInMinHighlyNonlinear( kmax )
T=1;
xo(:,1)=[0 0 0.15433*sin(-80*pi/180) 0.15433*cos(-80*pi/180)]';
for k=1:14
    xo(1,k+1)=xo(1,k)+0.15433*sin(-80*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(-80*pi/180)*T;
    xo(3,k+1)=0.15433*sin(-80*pi/180);
    xo(4,k+1)=0.15433*cos(-80*pi/180);
end

for k=15:kmax
    xo(1,k+1)=xo(1,k)+0.15433*sin(146*pi/180)*T;
    xo(2,k+1)=xo(2,k)+0.15433*cos(146*pi/180)*T;
    xo(3,k+1)=0.15433*sin(146*pi/180);
    xo(4,k+1)=0.15433*cos(146*pi/180);
end
X0=xo;
end