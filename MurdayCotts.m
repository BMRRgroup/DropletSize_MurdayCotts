%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script is an implementation of the differential equation including bessel functions as shown in Eq.5 
% of the paper: 
% Self-Diffusion Coefficient of Liquid Lithium
% J. S. Murday, and R. M. Cotts
% The Journal of Chemical Physics 48, 4938 (1968); doi: 10.1063/1.1668160
% 
% It uses the chebfun toolbox:
% https://www.chebfun.org/examples/roots/BesselRoots.html
% 
% The signal attenuation due to restricted diffusion behavior assuming a perfectly spherical shape can be calculated with this formulation.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created:	April 2018
%
% Modified: -
%
% Authors: 
%
%   mingming.wu@tum.de
% 
% --------------------------------
% Dimitrios C. Karampinos
%
% phone: +49 (0) 89 4140 6972
% email: dimitrios.karampinos@tum.de
% 
% web: https://bmrrgroup.de
%
% Body Magnetic Resonance Research Group
% Department of Diagnostic and Interventional Radiology
% Technical University of Munich
% Klinikum rechts der Isar
% 22 Ismaninger St., 81675 Munich
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% R: radius in m
% D: diffusion constant of medium at certain temperature in m^2/s
% b: b-value in s/mm^2
% gamma: 2*pi*gamma_bar in (rad*Hz/Tesla)
% G: gradient strength in T/m
% delta: gradient pulse duration in s
% Delta: diffusion time in s







function snormalized=MurdayCotts(G,R,gamma,Delta,delta,D,r)
 
alpha = r/R;
cnt_max = size(r);

Reihe = 0;
for mcnt = 2:cnt_max
    part1 = (alpha(mcnt)^4)*(alpha(mcnt)^2*R^2-2);
    expression1 = 2*delta; 
    expression2 = 2+exp(-alpha(mcnt)^2*D*(Delta-delta))-2*exp(-alpha(mcnt)^2*D*delta) -2*exp(-alpha(mcnt)^2*D*Delta) + exp(-alpha(mcnt)^2*D*(Delta+delta));
    expression3 = (alpha(mcnt)^2*D);
    part2 =  expression1 - expression2/expression3;
    total_part = part2/part1;
    Reihe = Reihe + total_part;
end


snormalized = exp(((-2*gamma^2*G^2)/D)*Reihe);