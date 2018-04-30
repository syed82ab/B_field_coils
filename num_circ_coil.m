function [B_ret,JacobianB_ret]=num_circ_coil(a,I,xx,yy,zz,x0,y0,z0)
% This function calculates the magammanetic field and its Jacobian at xx,yy,zz,
% for a circular coil at position x0,y0,z0 with its axis in +z.
% The radius of the coil is a0 with current I in Amperes.
% All distances are measured in metres. B is in Tesla, Jacobian in T/m.
%
% To rotate the coils in to another axis, swap the coordinates,
% e.gamma.: To place a coil on the axis pointing in -z use 
%       num_circ_coil(a0,-xx,yy,-zz,0,0,0).
%       To place the coil at z=10 pointing in -z 
%       num_circ_coil(a0,-xx,yy,-zz,0,0,-10).
%       To place the coil at x=10 pointing in x 
%       num_circ_coil(a0,-zz,yy,xx,10,0,).
%       etc.
%
% % Reference from 
% James Simpson, John Lane, Christopher Immer, and Robert Youngquist
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf
mu_0=4*pi*1e-7; %Tesla

xx=xx-x0;yy=yy-y0;zz=zz-z0;
rho_sq=xx.^2+yy.^2;
r_sq=xx.^2+yy.^2 +zz.^2;
gamma=xx.^2-yy.^2;
alp_c_sq=a.^2 + r_sq - 2*a*sqrt(rho_sq);
bet_c=sqrt(a.^2 + r_sq + 2*a*sqrt(rho_sq));
k_c_sq=1 - alp_c_sq./(bet_c.^2);

[K,E]=ellipke(k_c_sq);

const=mu_0*I/pi;
Bx=const*xx.*zz./(2*alp_c_sq.*bet_c.*rho_sq).*...
  ((a^2+r_sq).*E-alp_c_sq.*K);
By=const*yy.*zz./(2*alp_c_sq.*bet_c.*rho_sq).*...
  ((a^2+r_sq).*E-alp_c_sq.*K);
Bz=const./(2*alp_c_sq.*bet_c).*...
  ((a^2-r_sq).*E+alp_c_sq.*K);
B_ret(:,:,:,1)=Bx;
B_ret(:,:,:,2)=By;
B_ret(:,:,:,3)=Bz;

db_x_x= const*zz./(2*alp_c_sq.^2.*bet_c.^3.*rho_sq.^2).*...
  ((a^4*(-gamma.*(3*zz.^2 +a^2) + rho_sq.*(8*xx.^2-yy.^2))-...
   a^2*(rho_sq.^2.*(5*xx.^2 +yy.^2)-2*rho_sq.*zz.^2.*(2*xx.^2+yy.^2)+3*zz.^4.*gamma)-...
   r_sq.^2.*(2*xx.^4 +gamma.*(yy.^2 +zz.^2))).*E +...
   (a^2*(gamma.*(a^2 +2*zz.^2)-rho_sq.*(3*xx.^2-2*yy.^2))+...
   r_sq.*(2*xx.^4 +gamma.*(yy.^2 +zz.^2))).*alp_c_sq.*K);

db_x_y = const*xx.*yy.*zz./(2*alp_c_sq.^2.*bet_c.^3.*rho_sq.^2).*...
  ((3*a^4*(3*rho_sq-2*zz.^2) - r_sq.^2.*(2*r_sq+rho_sq)-...
   2*a^6-2*a^2*(2*rho_sq.^2-rho_sq.*zz.^2+3*zz.^4)).*E +...
  (r_sq.*(2*r_sq+rho_sq)-a^2*(5*rho_sq-4*zz.^2)+2*a^4).*alp_c_sq.*K);

db_x_z = const*xx./(2*alp_c_sq.^2.*bet_c.^3.*rho_sq).*...
  (((rho_sq-a^2).^2.*(rho_sq+a^2)+...
   2*zz.^2.*(a^4-6*a^2*rho_sq+rho_sq.^2)+zz.^4.*(a^2+rho_sq)).*E -...
   ((rho_sq-a^2).^2+zz.^2.*(rho_sq+a^2)).*alp_c_sq.*K);

db_y_x=db_x_y;

db_y_y = const*zz./(2*alp_c_sq.^2.*bet_c.^3.*rho_sq.^2).*...
  ((a^4*(gamma.*(3*zz.^2 +a^2) + rho_sq.*(8*yy.^2-xx.^2))-...
   a^2*(rho_sq.^2.*(5*yy.^2 +xx.^2)-2*rho_sq.*zz.^2.*(2*yy.^2+xx.^2)-3*zz.^4.*gamma)-...
   r_sq.^2.*(2*yy.^4 -gamma.*(xx.^2 +zz.^2))).*E +...
   (a^2*(-gamma.*(a^2 +2*zz.^2)-rho_sq.*(3*yy.^2-2*xx.^2))+...
   r_sq.*(2*yy.^4 -gamma.*(yy.^2 +zz.^2))).*alp_c_sq.*K);
 
db_y_z= const*yy./(2*alp_c_sq.^2.*bet_c.^3.*rho_sq).*...
  (((rho_sq-a^2).^2.*(rho_sq+a^2)+...
   2*zz.^2.*(a^4-6*a^2*rho_sq+rho_sq.^2)+zz.^4.*(a^2+rho_sq)).*E -...
   ((rho_sq-a^2).^2+zz.^2.*(rho_sq+a^2)).*alp_c_sq.*K);
 
db_z_x= db_x_z;
db_z_y= db_y_z;
db_z_z= const*zz./(2*alp_c_sq.^2.*bet_c.^3).*...
       ((6*a^2*(rho_sq-zz.^2)-7*a^4 +r_sq.^2).*E+...
       alp_c_sq.*(a^2-r_sq).*K);
     
 
JacobianB_ret(:,:,:,1,1)=db_x_x;
JacobianB_ret(:,:,:,1,2)=db_x_y;
JacobianB_ret(:,:,:,1,3)=db_x_z;

JacobianB_ret(:,:,:,2,1)=db_y_x;
JacobianB_ret(:,:,:,2,2)=db_y_y;
JacobianB_ret(:,:,:,2,3)=db_y_z;

JacobianB_ret(:,:,:,3,1)=db_z_x;
JacobianB_ret(:,:,:,3,2)=db_z_y;
JacobianB_ret(:,:,:,3,3)=db_z_z;


