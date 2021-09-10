function y = Lossalongtheway(d,L,Vm,lambda)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
global g
y=lambda*(L*Vm^2)/(d*2*g);
end

