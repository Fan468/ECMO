function [C, varargout] = JunctionLossCoefficient(U, A, theta)
% This function estimates the pressure loss coefficients at a vascular or pipe junction
% Inputs:
%   U = velocities (+ supplier, - collector) 
%   A = vessel areas
%   theta = vessel angles relative to some arbitrary reference
% Outputs:
%   C = Loss coefficients, i.e. p_i = p_j + C_j(0.5*rho*u_i^2)  
%   K = Loss coefficients, i.e. p_i + 0.5*rho*u_i^2 = p_j + 0.5*rho*u_j^2 + K_j(0.5*rho*u_i^2)  
% Usage:
%   - Diverging flow in a T-junction 
%   C = JunctionLossCoefficient([100, -50, -50], [2, 2, 2], [pi, 0, pi/2])  
%   - Output K loss coefficients (only applicable for three-branch junctions
%   [~,K] = JunctionLossCoefficient([100, -50, -50], [2, 2, 2], [pi, 0, pi/2])  
%
% Note: This code accompanies the article titled "A unified method for estimating pressure losses at vascular junctions", 
%       International Journal of Numerical Methods in Biomedical Engineering (IJNMBE)
%       DOI:10.1002/cnm.2717
% Authors: 
% Jonathan P. Mynard (1,2,3), Kristian Valen-Sendstad (1)
% (1) Biomedical Simulation Laboratory, Department of Mechanical and Industrial Engineering, University of Toronto, Toronto, Canada
% (2) Heart Research, Clinical Sciences, Murdoch Childrens Research Institute, Parville VIC, Australia
% (3) Department of Paediatrics, University of Melbourne, Parkville VIC, Australia
% Correspondence to: jonathan.mynard@mcri.edu.au

% Revision 1 (20/8/2016) by Jonathan Mynard
% - Removed error (factor of 2) in denominator of Line 60 (calculation of TotPseudoArea). 
%   Note that the code used in the IJNMBE paper was correct.
% - Added comment regarding factor (1-exp(-FlowRatio/0.02)) on line 71
% - Updated help examples to use the correct function name (JunctionLossCoefficient rather than LossCoefficient)

U = U(:);
A = A(:);
theta = wrapToPi(theta(:));
theta = (theta(:));

Q = U.*A;                    % Volumetric Flow
Ci = Q < 0.0;                % Collector indexes
Si = Q >= 0.0;               % Supplier indexes 
Qtot = sum(Q(Si));           % Total flow through junction
FlowRatio = -Q(Ci)/Qtot;     % Flow ratio of collectors

%% Reorient all branch angles so that the average collector angle is 0 
PseudoColAngle = mean(theta(~Si));   
PseudoSupAngle = atan2(sum(sin(theta(Si)).*Q(Si)), sum((cos(theta(Si))).*Q(Si)));
if abs(PseudoSupAngle-PseudoColAngle) < pi/2  % Ensure that PseudoSupAngle will be in the second quadrant
   PseudoColAngle = PseudoColAngle + pi;
end
theta = wrapToPi(theta - PseudoColAngle); 

%% Calculate the pseudosupplier angle
pseudodirection = sign(mean(sin(theta(Si)).*Q(Si)));  % is the majority of supplier flow coming from positive or negative angles?
if pseudodirection < 0  % flip angles to ensure the majority of supplier flow is coming from positive angles
   theta = -theta;
end
PseudoSupAngle = atan2(sum(sin(abs(theta(Si))).*Q(Si)), sum(cos(abs(theta(Si))).*Q(Si)));  % always positive

%% Calculate effective pseudosupplier area
etransferfactor = (0.8*(pi-PseudoSupAngle)*sign(theta(Ci))-0.2).*(1-FlowRatio);
TotPseudoArea = Qtot ./ ((1 - etransferfactor) .* sum(U(Si).*Q(Si))/Qtot);  % Denominator is pseudosupplier velocity  

%% Calculate area ratios and relative angles
AreaRatio = TotPseudoArea ./ A(Ci);
phi = wrapTo2Pi(PseudoSupAngle - theta(Ci));

%% Calculate the C loss coefficients
C = zeros(length(U),1);
% Note: The factor (1-exp(-FlowRatio/0.02)) avoids infinite C when FlowRatio approaches zero.
% This factor is essentially equal to 1 when FlowRatio > 0.1. This factor was mistakenly missed
% from discussion in the IJNMBE paper.
C(Ci) = (1-exp(-FlowRatio/0.02)).*(1 - (1.0./(AreaRatio.*FlowRatio)) .* cos(0.75*(pi - phi)));

%% Calculate the K loss coefficients (if required)
if nargout > 1
    if length(U) <= 3  % cannot be applied in general to junctions with more than three branches
        if sum(Ci) == 1  % converging flow
            Ucom = U(Ci);
        else                % diverging flow
            Ucom = U(Si);
        end
            
        K = (U(Ci).^2/Ucom^2) .* (2*C(Ci) + U(Si).^2./U(Ci).^2 - 1);
        varargout{1} = K;
    else
        varargout{1} = [];
    end
end
