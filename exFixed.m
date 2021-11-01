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



%----------------------------%
% PART 1 - R0 X V0 DIRECTION %
%----------------------------%
% Assume an impulse in the r0 x v0 direction...
Idir = cross(r,v) / norm( cross(r,v) );
u_h0 = Idir;
v_h0 = v +  u_h0*Imag;
coe_h0 = coe_from_sv(r,v_h0,mu);



%--------------------------------------------%
% PART 2 - INTEGRATE THE R,V ATTACHED SYSTEM %
%--------------------------------------------%
% Is the Frobenius condition satisfied? It is trivially satisfied because
% the number of columns in G is 1, i.e., the second column needed in the
% stated test is zero.
ops   = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~,z] = ode45(@attached_rv,[0,1],[r';v'],ops,Idir',Imag,1);
v_att = z(end,4:6);
coe_att = coe_from_sv(r,v_att,mu);


%----------------------------------------------%
% PART 3 - INTEGRATE THE GAUSS ATTACHED SYSTEM %
%----------------------------------------------%
[~,s] = ode45(@attached_gauss,[0,1],[sma,ecc,inc,W,w,f],ops,Idir',Imag,1,mu);
coe_gatt = [sqrt(mu*s(end,1)*(1-s(end,2)^2)), s(end,2), s(end,4), s(end,3), s(end,5), s(end,6), s(end,1)];
[r_gatt,v_gatt] = sv_from_coe(coe_gatt,mu);


%--------------------------------%
% PART 4 - GROUP ALL THE RESULTS %
%--------------------------------%
% The first row should be different than all subsequent rows, which should
% be equal.
v_all = [v_h0; v_att; v_gatt]
coe_all = [coe_h0; coe_att; coe_gatt]

%----------------------%
% SUPPORTING FUNCTIONS %
%----------------------%
function xdot = two_body(~,x,mu,Tdir,Tmag)
r = x(1:3);
v = x(4:6);
thrust = Tmag*Tdir;
rdot = v;
vdot = -mu*r/norm(r)^3 + thrust;
xdot = [rdot; vdot];
end


function zdot = gauss_eqs(~,z,mu,Tdir,Tmag)
a = z(1); e = z(2); inc = z(3);
W = z(4); w = z(5); f = z(6);
theta = w+f;
h = sqrt(mu*a*(1-e^2)); p = h^2/mu; r = p/(1+e*cos(f));

[R,V] = sv_from_coe([h e W inc w f a],mu);
er = R/norm(R); ar = dot(Tdir,er)*Tmag;
eh = cross(R,V)/norm( cross(R,V) ); ah = dot(Tdir,eh)*Tmag;
et = cross(eh,er); at = dot(Tdir,et)*Tmag;

adot = 2*a^2/h*e*sin(f)*ar + 2*a^2/h*p/r*at;
edot = p/h*sin(f)*ar + 1/h*(p*cos(f)+r*cos(f)+r*e)*at;
idot = r*cos(theta)/h*ah;
Wdot = r*sin(theta)/(h*sin(inc))*ah;
wdot = 1/(h*e)*(-p*cos(f)*ar+(p+r)*sin(f)*at)-r*sin(theta)*cos(inc)/(h*sin(inc))*ah;
fdot = h/r^2 + 1/(h*e)*(p*cos(f)*ar-(p+r)*sin(f)*at);
zdot = [adot; edot; idot; Wdot; wdot; fdot];
end


function F = imp_dir(u,r0,v0,dv)
v1 = v0 + u*dv;
F(1) = norm(v0) / norm(v1) - 1;
F(2) = norm( cross(r0,v0) ) / norm( cross(r0,v1) ) - 1;
F(3) = norm(u) - 1;
end


function zdot = attached_rv(~,z,Idir,Imag,tf)
z1 = z(1:3);
z2 = z(4:6);
zdot = 1/tf*[zeros(3,1); Imag * Idir];
end

function zdot = attached_gauss(~,z,Idir,Imag,tf,mu)
a = z(1); e = z(2); inc = z(3);
W = z(4); w = z(5); f = z(6);
theta = w+f;
h = sqrt(mu*a*(1-e^2)); p = h^2/mu; r = p/(1+e*cos(f));

[R,V] = sv_from_coe([h e W inc w f a],mu);
er = R/norm(R); ar = dot(Idir,er)*Imag;
eh = cross(R,V)/norm( cross(R,V) ); ah = dot(Idir,eh)*Imag;
et = cross(eh,er); at = dot(Idir,et)*Imag;

adot = 2*a^2/h*e*sin(f)*ar + 2*a^2/h*p/r*at;
edot = p/h*sin(f)*ar + 1/h*(p*cos(f)+r*cos(f)+r*e)*at;
idot = r*cos(theta)/h*ah;
Wdot = r*sin(theta)/(h*sin(inc))*ah;
wdot = 1/(h*e)*(-p*cos(f)*ar+(p+r)*sin(f)*at)-r*sin(theta)*cos(inc)/(h*sin(inc))*ah;
fdot = 1/(h*e)*(p*cos(f)*ar-(p+r)*sin(f)*at);
zdot = 1/tf*[adot; edot; idot; Wdot; wdot; fdot];
end