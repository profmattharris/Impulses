clear all; close all; clc

%---------------%
% PRELIMINARIES %
%---------------%
% Earth Data
mu = 398600;
R  = 6378;

% Impulse magnitude (in km/s)
Imag = 1;

% Initial Orbital Elements
sma = R + 1000;
ecc = 0.1;
inc = 10*pi/180;
W   = 20*pi/180;
w   = 30*pi/180;
f   = 40*pi/180;
h   = sqrt( mu*sma*(1-ecc^2) );
coe = [h ecc W inc w f sma];

% Convert the initial OEs to R,V...
[r,v] = sv_from_coe(coe,mu);



%-----------------------%
% PART 1 - V0 DIRECTION %
%-----------------------%
% Assume an impulse in the v0 direction...
u_r0 = r/norm( r );
v_r0 = v +  u_r0*Imag;
coe_r0 = coe_from_sv(r,v_r0,mu);


%--------------------------------------------%
% PART 2 - INTEGRATE THE R,V ATTACHED SYSTEM %
%--------------------------------------------%
ops   = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~,z] = ode45(@attached_rv,[0,1],[r';v'],ops,Imag,1);
v_att = z(end,4:6);
coe_att = coe_from_sv(r,v_att,mu);


%----------------------------------------------%
% PART 3 - INTEGRATE THE GAUSS ATTACHED SYSTEM %
%----------------------------------------------%
[~,s] = ode45(@attached_gauss,[0,1],[sma,ecc,inc,W,w,f],ops,Imag,1,mu);
coe_gatt = [sqrt(mu*s(end,1)*(1-s(end,2)^2)), s(end,2), s(end,4), s(end,3), s(end,5), s(end,6), s(end,1)];
[r_gatt,v_gatt] = sv_from_coe(coe_gatt,mu);


%--------------------------------%
% PART 4 - GROUP ALL THE RESULTS %
%--------------------------------%
v_all = [v_r0; v_att; v_gatt]
coe_all = [coe_r0; coe_att; coe_gatt]

%----------------------%
% SUPPORTING FUNCTIONS %
%----------------------%
function xdot = two_body(~,x,mu,Tmag)
r = x(1:3);
v = x(4:6);
thrust = Tmag*r / norm(r);
rdot = v;
vdot = -mu*r/norm(r)^3 + thrust;
xdot = [rdot; vdot];
end


function zdot = gauss_eqs(~,z,mu,Tmag)
a = z(1); e = z(2); inc = z(3);
W = z(4); w = z(5); f = z(6);
theta = w+f;
h = sqrt(mu*a*(1-e^2)); p = h^2/mu; r = p/(1+e*cos(f));
adot = 2*a^2/h*e*sin(f)*Tmag;
edot = p/h*sin(f)*Tmag;
idot = 0;
Wdot = 0;
wdot = -p/(h*e)*cos(f)*Tmag;
fdot = h/r^2 +p/(h*e)*cos(f)*Tmag;
zdot = [adot; edot; idot; Wdot; wdot; fdot];
end


function F = imp_dir(u,r0,v0,dv)
v1 = v0 + u*dv;
F(1) = norm(v0) / norm(v1) - 1;
F(2) = norm( cross(r0,v0) ) / norm( cross(r0,v1) ) - 1;
F(3) = norm(u) - 1;
end


function zdot = attached_rv(~,z,Imag,tf)
z1 = z(1:3);
z2 = z(4:6);
zdot = 1/tf*[zeros(3,1); Imag * z1/norm(z1)];
end

function zdot = attached_gauss(~,z,Imag,tf,mu)
a = z(1);
e = z(2);
inc = z(3);
W = z(4);
w = z(5);
f = z(6);
theta = w+f;
h = sqrt(mu*a*(1-e^2));
p = h^2/mu;
r = p/(1+e*cos(f));
adot = 2*a^2/h*e*sin(f)*Imag;
edot = p/h*sin(f)*Imag;
idot = 0;
Wdot = 0;
wdot = -p/(h*e)*cos(f)*Imag;
fdot = p/(h*e)*cos(f)*Imag;
zdot = 1/tf*[adot; edot; idot; Wdot; wdot; fdot];
end