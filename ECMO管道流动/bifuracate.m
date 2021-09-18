function [outputArg1,outputArg2] =bifuracate(inputArg1,inputArg2)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
 [~,K] = JunctionLossCoefficient([100, -50, -50], [2, 2, 2], [pi, 0, pi/2]) 
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

