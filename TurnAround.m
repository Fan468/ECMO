function [a,b] = TurnAround(Re1,Re2)
% 转捩区域模型
%   给定转捩开始雷诺数和结束雷诺数，
%   返回 $\lambda$ =a*Re^b的系数（双对数下是线性增加）
b=log10(64*Re2^(0.25)/(Re1*0.3164))/log10(Re1/Re2);
a=64*Re1^(-1-b);
end

