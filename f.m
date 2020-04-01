function dx = f(t,x,mu,options)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
mu = 4902.799;
r=x(1:3);
v=x(4:6);
ar=norm(r);
dx=[v;-mu*r/ar^3];
end

