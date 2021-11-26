function [t,z] = SolveKineticModel(z0, k, type, tspan)
if (nargin < 3)
    type = 'batch';
end
if nargin < 4
    % Choose time span
    tspan = [0 100]; % in minutes
end
% x'(t) = f(k,x(t))
% x - vector of component concentrations
% k - vector of parameters



options = odeset('RelTol',0.05, 'AbsTol', 0.05,'InitialStep',1,'NonNeg',1:6);

% Call ODE function
if (strcmp(type,'batch'))
    [t, z] = ode23t(@BatchKineticsModel, tspan, z0, options,  k);
elseif (strcmp(type,'membrane'))
    [t, z] = ode23t(@MembraneKineticsModel, tspan, z0, options,  k);
else
    error('SolveKineticModel:UnsupportedModel','Only batch and membrane models are supported');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Membrane model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = MembraneKineticsModel(t,z,k)
% Equations from the kinetic model of reaction

[MW_TG, MW_MEOH , rho_permeate, F_TGin, F_DGin, F_MGin, F_MEOHin, V_R] = getReactionConstans();

v_out = (F_MEOHin*MW_MEOH+F_TGin*MW_TG)/rho_permeate; % v_out is 0.6277 L/min, volumetric flow rate of the permeate, L/min. density of permeate phase is 0.815g/ml

%z(3)=TG; z(4)=DG; z(5)=MG; z(2)=GLY; z(1)=MEOH; z(6)=FAME

dz = zeros(6,1); 
dz(1) = (F_MEOHin-v_out*z(1))/V_R - k(1) * z(3) * z(1) + k(2) * z(4) * z(6) - k(3) * z(4) * z(1) + k(4) * z(5) * z(6) - k(5) * z(5) * z(1) + k(6) * z(2) * z(6);
dz(2) = -((z(2)*v_out)/V_R)+ k(5) * z(5) * z(1) - k(6) * z(2) * z(6);
dz(3) = (F_TGin-v_out*z(3))/V_R - k(1) * z(3) * z(1) + k(2) * z(4) * z(6);
dz(4) = (F_DGin-v_out*z(4))/V_R + k(1) * z(3) * z(1) - k(2) * z(4) * z(6) - k(3) * z(4) * z(1) + k(4) * z(5) * z(6);
dz(5) = (F_MGin-v_out*z(5))/V_R + k(3) * z(4) * z(1) - k(4) * z(5) * z(6) - k(5) * z(5) * z(1) + k(6) * z(2) * z(6);
dz(6) = -((z(6)*v_out)/V_R)+ k(1) * z(3) * z(1) - k(2) * z(4) * z(6) + k(3) * z(4) * z(1) - k(4) * z(5) * z(6) + k(5) * z(5) * z(1) - k(6) * z(2) * z(6);

end

function [MW_TG, MW_MEOH, rho_permeate, F_TGin, F_DGin, F_MGin, F_MEOHin, V_R] = getReactionConstans()
MW_TG = 885.45; %MW of TG = 885.45 g/mol
MW_MEOH = 32.04; %MW of MEOH = 32.04 g/mol

rho_permeate = 815; %Density of permeate phase = 815 g/L, extracted from literature "methanol recyclying ..."

zlh_TGin= 3 ;     %Volumetric inlet flow of TG = 3L/h
zlh_MEOHin= 3 ;   %Volumetric inlet flow of MEOH = 3L/h

zml_TGin = 0.482123;      %concentration of TG in the inlet flow (mole/l) from spreadsheet
zml_DGin = 0.036942;      %concentration of DG in the inlet flow (mole/l) from spreadsheet
zml_MGin = 0.001966;      %concentration of MG in the inlet flow (mole/l) from spreadsheet
zml_MEOHin = 12.18462;    %concentration of MEOH in the inlet flow (mole/l) from spreadsheet

F_TGin = zml_TGin*zlh_TGin/60;  % zml (mol/l); zlh (l/h)
F_DGin = zml_DGin*zlh_TGin/60;
F_MGin = zml_MGin*zlh_TGin/60;
F_MEOHin = zml_MEOHin*zlh_MEOHin/60; % mol/min
V_R = 6; % volume of the membrane reactor, L

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Batch model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dz = BatchKineticsModel(t,z,k)
% Equations from the kinetic model of reaction

dz=zeros(6,1);
dz(1) =  -k(1) * z(3) * z(1) + k(2) * z(4) * z(6) - k(3) * z(4) * z(1) + k(4) * z(5) * z(6) - k(5) * z(5) * z(1) + k(6) * z(2) * z(6);
dz(2) = k(5) * z(5) * z(1) - k(6) * z(2) * z(6);
dz(3) = - k(1) * z(3) * z(1) + k(2) * z(4) * z(6);
dz(4) = k(1) * z(3) * z(1) - k(2) * z(4) * z(6) - k(3) * z(4) * z(1) + k(4) * z(5) * z(6);
dz(5) = k(3) * z(4) * z(1) - k(4) * z(5) * z(6) - k(5) * z(5) * z(1) + k(6) * z(2) * z(6);
dz(6) = k(1) * z(3) * z(1) - k(2) * z(4) * z(6) + k(3) * z(4) * z(1) - k(4) * z(5) * z(6) + k(5) * z(5) * z(1) - k(6) * z(2) * z(6);
end
