%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script extracts the roots of a differential equation including bessel functions as shown in Eq.5 
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

% to be adjusted: maxval
% The roots are calculated in a range of 0 to 'maxval'
% The function becomes instable at a certain moment, probably not more than
% a hundred roots are necessary for a good approximation of the diffusion
% behavior anyways. 

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
% when typing in following two lines in the matlab command window, the equation is translated to the equivilant long expression with
% cos and sin functions:
%  syms x
%  x*diff(besselj(1.5,x))-0.5*besselj(1.5,x)

%% using chebfun
maxval = 10000;
J0 = chebfun(@(x) x*((2^(1/2)*(cos(x) - sin(x)/x))/(2*x^(3/2)*pi^(1/2)) + (2^(1/2)*(sin(x) + cos(x)/x - sin(x)/x^2))/(x^(1/2)*pi^(1/2))) + (2^(1/2)*(cos(x) - sin(x)/x))/(2*x^(1/2)*pi^(1/2)), [0 maxval]);
r = roots(J0);
figure, plot(J0,'linewidth',1.6), grid on
title('Bessel function J_0','fontsize',16)
hold on, plot(r, J0(r), '.r', 'markersize', 14)

save('roots_besselfunc_e4.mat', 'r', 'J0');




