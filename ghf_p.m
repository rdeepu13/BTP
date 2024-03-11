function [GH_points, GH_Weights] = ghf_p(alpha)
J = zeros(alpha);
% POINT GENERATION
for i = 1:alpha-1
    J(i,i+1) = sqrt(i/2);
    J(i+1,i) = sqrt(i/2);
end

[V,D]=eig(J);
points = zeros(1, alpha);
for i=1:alpha
    points(i)=sqrt(2)*D(i,i);
end
root = zeros(1, alpha); w = zeros(1, alpha);
for j=1:alpha
    s=0;
    for i=1:alpha
        s=s+(V(i,j)*V(i,j));
    end
    root(j)=sqrt(s);
    w1=V(1,j)/root(j);
    w(j)=w1*w1;
end

QP = []; % size of QP n*alpha^n
for i1=1:alpha
    for j1=1:alpha
        for k1=1:alpha
            for k11 = 1:alpha
                QP=[QP [points(i1);points(j1);points(k1);points(k11)]];
            end
        end
    end
end
GH_points = QP;
W = []; %size of W will be 1*alpha^n
for i1=1:alpha
    for j1=1:alpha
        for k1=1:alpha
            for k11 = 1:alpha
                W=[W w(i1)*w(j1)*w(k1)*w(k11)];
            end
        end
    end
end
GH_Weights = W;
end
