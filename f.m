function dx = f(t,x,mu,options)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
mu = 4902.799;
r=x(1:3);
v=x(4:6);
ar=norm(r);
dx=[v;-mu*r/ar^3];
end

