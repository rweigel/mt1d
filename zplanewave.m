function Cs = zplanewave(s,h,f)
%ZPLANEWAVE Surface impedance for plane wave incident on infinite half space.
%
%   C = ZPLANEWAVE(s,h,f)
%   
%   Surface transfer function C in [m] at positive frequencies f [Hz]
%   above layered conductors with conductivities s [1/Ohm-m], thicknesses
%   h [m], and geometry
%
%   Apparent resistivity is rho_a = |C|^2 mu_0*2*pi*f [Ohm-m]
%
%   Conductivity    Thickness
%
%   0               Inf
%   --------------------------  Surface
%   s(1)            h(1)
%   --------------------------
%   s(2)            h(2)
%   --------------------------
%   .
%   .
%   .
%   --------------------------
%   s(end-1)        h(end-1)
%   --------------------------  Top of bottom layer
%   s(end)          h(end)
%
%   If h(end) is not Inf, an extra infinite bottom layer with s = s(end)
%   is assumed.
%   
%   Computed using Eqns. 2.33 and 2.34 of Simpson and Bahr, Practical
%   Magnetotellurics, 2005 (Eqn. 2.33 is recursion formula from Wait 1954).
%   All layers have vacuum permeability.
%
%   For an infinite half-layer with s = [s0], h = [Inf], and f = [f0],
%       C = zplanewave(s,h,f)
%   returns
%       C = 1/q, where
%       q = sqrt(mu_0*s0*2*pi*f0)*(1+i)/sqrt(2)
%   and
%   C = Ex/(i*2*pi*f*By) = -Ey/(i*2*pi*f*Bx)
% 
%   See also
%       ZPLANEWAVE_TEST - Verify calculations
%       ZPLANEWAVE_DEMO - Plot transfer fn, phase, and impulse response

% Code assumes f is vector with columns
if size(f,1) > size(f,2)
    f = f';
end

if size(f,1) > 1 && size(f,2) > 1
    error('frequency must be an array');
end

if ~isinf(h(end))
    % Set infinite layer to have same value as last layer.
    if (0)
    fprintf('zplanewave:\n  Warning: Last layer does not have infinite thickness\n');
    fprintf('               Inserting infinite bottom layer with conductivity\n');
    fprintf('               s = s(end)\n');
    s(end+1) = s(end);
    end
end

if any(f < 0)
  error('Frequencies must be positive');
end

nl    = length(s);      % Number of layers  
mu_0  = 4*pi*1e-7;      % Vacuum permeability

% Columns of C are frequencies.  Compute C at top of bottom layer.
    
% Eqn. 2.34 of Simpson and Bahr
q(nl,:) = sqrt(j*mu_0*s(nl)*2*pi*f);
C(nl,:) = 1./q(nl,:);

% Compute C(nl-1) given C(nl), C(nl-2) given C(nl-1), etc.
% (Eqn. 2.33 of Simpson and Bahr)
for l = nl-1:-1:1
    q(l,:) = sqrt(j*mu_0*s(l)*2*pi*f);  % Eqn. 2.28
    C(l,:) = (1./q(l,:)).*( q(l,:).*C(l+1,:) + tanh(q(l,:)*h(l)) ) ...
                        ./ (1 + q(l,:).*C(l+1,:).*tanh(q(l,:)*h(l)));
end

% Surface impedance.
Cs = C(1,:);

% Return Z having same shape as f.
if (size(f,1) > size(f,2))
    Cs = Cs';
end
