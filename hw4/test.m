%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exact solution to Sod Shock Tube problem
%
% Boundary conditions:
% pL = 1, rhoL = 1, uL = 0
% pR = 0.125, rhoR = 0.1, uR = 0
%
% Script by Sebastiaan ten Pas and Martin Goossens
% 29-1-2015 | University of Twente
%
% For more information: info@sebastiaantenpas.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, x2, x3, x4, rho, u, p] = exactSOD(rhoL, uL, pL, rhoR, uR, pR, gamma, x, x0, t)
   
aL = sqrt(gamma * pL / rhoL); % The speed of sound at region L
aR = sqrt(gamma * pR / rhoR); % The speed of sound at region R

f_p1R = @(p1R) pR / pL - ( 1 - aR/aL * (gamma - 1)/(2 * gamma) * (p1R - 1)/sqrt( 1 + (gamma + 1)/(2 * gamma) * (p1R - 1) ) )^(2 * gamma / (gamma - 1)) / p1R;
p1R = fzero(f_p1R, 3); % The ratio between p1 and pR

up = aR / gamma * (p1R - 1) / sqrt(1 + (gamma + 1) / (2 * gamma) * (p1R - 1)); % The speed at the interface
us = aR * ((gamma + 1) * up / (4 * aR) + sqrt( 1 + ( (gamma + 1) * up / (4 * aR) )^2)); % The speed of the shock

a2 = aL - (gamma - 1) / 2 * up; % The speed of sound at region 2

x1 = x0 - aL * t;
x2 = x0 + (up - a2) * t;
x3 = x0 + up * t;
x4 = x0 + us * t;
 
for i = 1:numel(x)
    if x(i) <= x1 %Region L
        rho(i) = rhoL;
        u(i) = uL;
        p(i) = pL;
    elseif (x(i) <= x2) %Region E
        p(i) = pL * (2 / (gamma + 1) + (gamma - 1) / (gamma + 1) * (x0 - x(i)) / (aL * t))^(2 * gamma / (gamma - 1));
        rho(i) = rhoL * (2 / (gamma + 1) + (gamma - 1) / (gamma + 1) * (x0 - x(i)) / (aL * t))^(2 / (gamma - 1));
        u(i) = 2 / (gamma + 1) * aL + 2 / (gamma + 1) * (x(i) - x0) / t;
    elseif (x(i) <= x3) %Region 2
        u(i) = up;
        p(i) = p1R * pR;
        rho(i) = rhoL * (1 - (gamma - 1) / 2 * up / aL)^(2 / (gamma - 1));
    elseif (x(i) <= x4) %Region 1
        u(i) = up;
        rho(i) = rhoR * (1 + (gamma + 1) / (2 * gamma) * (p1R - 1)) / (1 + (gamma - 1) / (2 * gamma) * (p1R - 1));
        p(i) = p1R * pR;
    else %Region R
        rho(i) = rhoR;
        u(i) = uR;
        p(i) = pR;
    end
end

end