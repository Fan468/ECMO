% 计算ECMO的流量分离比
%% 主要的参数
%长度,unit=cm
cmTom=0.01;
cmTomm=10;
L0=15;
a1=15;
a2=30;
a3=2;
Ha=50;
La=18;
v1=15;
v2=30;%可更改
v3=2;
Lv=50;

%直径,unit=mm
inchtocm=2.54;
FrTomm=1/3;
Phi0=0.375*inchtocm*cmTomm;
Da=19*FrTomm;%可更改
Dv=21*FrTomm;%可更改


%流体参数
rho=1065;%单位kg/m^3;
mu=0.004;%单位kg/(m*s);
g=9.8;%单位m/s^2;
%边界条件
P0UT_a=101325;%单位Pa
P0UT_v=101325;%单位pa
UIN=0.468;%入口速度m/s 等价于入口流量2L/min 可更改

%% 统一单位及几何角度等


rad2deg(atan((Phi0-Dv)/(2*v3)));



