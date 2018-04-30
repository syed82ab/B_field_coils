% Reference from 
% James Simpson, John Lane, Christopher Immer, and Robert Youngquist
% https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf

% Computes  the  a field  from  2  thin  loops  of  radius  a,
% located  at  z  = +|d2|  and  -|d1|,  
% current  I.
% and a third coil located at x=-|d3| with radius a3 and current I3

plot1d=0;
plot2d=3;
plot_fg=0;
save_plot=0;
%% Define coil parameters %%
LN = 48;             % Number  of  turns  in  the  loops
mu_0 = 4*pi*1e-7;   % Tesla
I = 0.0;      % Amperes in coils 1 and 2
%I = 60.0; 
a = 0.0502;         % Radius 1 and 2 (Metres)
% d=0.0596/2;         % coil position.
d=(0.067+2*4*0.00208)/2;
d1 = -d;            % Coil 1 at z= - d; Current loop point towards z.
d2 = d;             % Coil 2 at z= + d; Current loop point towards -z.
for I3= 0:5:20;              % Amperes in coil 3
% I3=0;             % Amperes in coil 3
% for LN3=32:5:64    
LN3=6*5;             % Number  of  turns  in  the  coil3
LN4=6*5;             % Number  of  turns  in  the  coil3
dI=(0.083+2*3*0.00208)/2;
aI=(0.033+0.033+2*5*cos(pi/6)*0.00208)/2;
a3 = aI;        % Radius 3
d3 = dI;         % Coil 3 at x= |d3|; Current loop point towards x.
a4 = aI;        % Radius 4
d4 = -dI;         % Coil 4 at x= -|d4|; Current loop point towards x.
axis_shift=0.000;

%% Define calculation parameters %%
% Define  the  plane  over  which  fields  are  computed
% N must  be  odd  to  include  the  point  (0,0)
% Use even N to avoid divide by zeros.
M=50;  % No.  of  points  along  the  x,y  axis
N=50;  % No.  of  points  along  the  z  axis
xr=0.02;
x_range=[-xr xr]; % in meters
y_range=[-xr xr];
z_range=[-xr xr];

%% Generate Axes for calculations %%
x_line=linspace(x_range(1),x_range(2),M);  
y_line=linspace(y_range(1),y_range(2),M);  
z_line=linspace(z_range(1),z_range(2),N);
[xx,yy,zz] = ndgrid(x_line,y_line,z_line);

%% Base calculations here with positions and axis included manually %%
% Simplifying substitutions. 
% For 1st coil at (z=-|d1|), 
% Cartesian: x,y,z == xx,yy,zz
rho_sq_1=xx.^2+yy.^2;
r_sq_1=xx.^2+yy.^2+(zz-d1).^2; %% Displace coil by z=d1 from origin.
gamma_1=xx.^2-yy.^2;

alp_c_sq_1=a^2 + r_sq_1 -2*a*sqrt(rho_sq_1);
bet_c_sq_1=a^2 + r_sq_1 +2*a*sqrt(rho_sq_1);
k_c_sq_1=1-alp_c_sq_1./bet_c_sq_1;

% For 2st coil at (z=+|d2|),
% Cartesian: x,y,z == -xx,yy,-zz
rho_sq_2=(-xx).^2+yy.^2; %%Displace z by d2, invert x,z-axes.
r_sq_2=(-xx).^2+yy.^2+(-(zz-d2)).^2; 
gamma_2=(-xx).^2-yy.^2;

alp_c_sq_2=a^2 + r_sq_2 -2*a*sqrt(rho_sq_2);
bet_c_sq_2=a^2 + r_sq_2 +2*a*sqrt(rho_sq_2);
k_c_sq_2=1-alp_c_sq_2./bet_c_sq_2;

% For 3rd coil at (x=-|d3|) with axis along x
% Cartesian:  x,y,z == -zz,yy,xx
rho_sq_3=(-zz).^2+(yy+axis_shift).^2; %% Rotate x,z-axis by pi/2 keeping y-axis fixed
r_sq_3=(-zz).^2+(yy+axis_shift).^2+(xx-d3).^2; %% Displace x-axis by d3.
gamma_3=(-zz).^2-(yy+axis_shift).^2;

alp_c_sq_3=a3^2 + r_sq_3 -2*a3*sqrt(rho_sq_3);
bet_c_sq_3=a3^2 + r_sq_3 +2*a3*sqrt(rho_sq_3);
k_c_sq_3=1-alp_c_sq_3./bet_c_sq_3;

% For 4th coil at (x=-|d4|) with axis along x
% Cartesian:  x,y,z == -zz,yy,xx
rho_sq_4=(-zz).^2+(yy+axis_shift).^2; %% Rotate x,z-axis by pi/2 keeping y-axis fixed
r_sq_4=(-zz).^2+(yy+axis_shift).^2+(xx-d4).^2; %% Displace x-axis by d4.
gamma_4=(-zz).^2-(yy+axis_shift).^2;

alp_c_sq_4=a4^2 + r_sq_4 -2*a4*sqrt(rho_sq_4);
bet_c_sq_4=a4^2 + r_sq_4 +2*a4*sqrt(rho_sq_4);
k_c_sq_4=1-alp_c_sq_4./bet_c_sq_4;

% Common constants
C=mu_0*LN*I/pi;
C3=mu_0*LN3*I3/pi;


%% Determine  E1liptic  integrals in cartesian %%
[K1,E1] = ellipke(k_c_sq_1);
[K2,E2] = ellipke(k_c_sq_2);
[K3,E3] = ellipke(k_c_sq_3);
[K4,E4] = ellipke(k_c_sq_4);

%% B-xyz fields anonymous function for coil with positive current at origin
%% alp2 == alpha.^2, rho2== rho.^2 , etc.
bx =@(C,a,x,y,z,alp2,bet,rho2,r2,E,K) C*x.*(z)./(2*alp2.*bet.*rho2).*...
     ((a^2+r2).*E - alp2.*K);
by =@(C,a,x,y,z,alp2,bet,rho2,r2,E,K) C*y.*(z)./(2*alp2.*bet.*rho2).*...
     ((a^2+r2).*E - alp2.*K);
bz =@(C,a,x,y,z,alp2,bet,rho2,r2,E,K) C./(2*alp2.*bet).*...
     ((a^2-r2).*E + alp2.*K);


%% Compute  B-xyz fields
% 1st coil at z=-|d1| and axis along z
bx1=bx(C,a,xx,yy,zz-d1... %% Displace by d1 in z-axis
       ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,E1,K1);
by1=by(C,a,xx,yy,zz-d1...
       ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,E1,K1);
bz1=bz(C,a,xx,yy,zz-d1...
       ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,E1,K1);
bx1(~isfinite(bx1))=0; % get rid of infinites if any
by1(~isfinite(by1))=0;

% 2nd coil at z=d2 and axis along -z
bx2=bx(C,a,-xx,yy,-(zz-d2)...
      ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,E2,K2);
by2=by(C,a,-xx,yy,-(zz-d2)...
      ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,E2,K2);
bz2=bz(C,a,-xx,yy,-(zz-d2)...
      ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,E2,K2);
bx2(~isfinite(bx2))=0;
by2(~isfinite(by2))=0;

% 3rd coil at x=-|d3| and axis along x
bx3=bx(C3,a3,-zz,yy+axis_shift,xx-d3...
      ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,E3,K3);
by3=by(C3,a3,-zz,yy+axis_shift,xx-d3...
      ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,E3,K3);
bz3=bz(C3,a3,-zz,yy+axis_shift,xx-d3...
      ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,E3,K3);
bx3(~isfinite(bx3))=0;
by3(~isfinite(by3))=0;

% 3rd coil at x=|d4| and axis along x
bx4=bx(C3,a4,-zz,yy+axis_shift,xx-d4...
      ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,E4,K4);
by4=by(C3,a4,-zz,yy+axis_shift,xx-d4...
      ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,E4,K4);
bz4=bz(C3,a4,-zz,yy+axis_shift,xx-d4...
      ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,E4,K4);
bx4(~isfinite(bx4))=0;
by4(~isfinite(by4))=0;

%% Add up the field components in the correct axes %%
% x1,y1,z1 == xx,yy,zz ; x2,y2,z2 == -xx,yy,-zz ;
% x3,y3,z3 == -zz,yy,xx ; x4,y4,z4 == -zz,yy,xx

bax = (bx1 - bx2 + bz3 + bz4); 
bay = (by1 + by2 + by3 + by4);
baz = (bz1 - bz2 - bx3 - bx4); 

bat = sqrt(bax.^2 + bay.^2 +baz.^2);
%% Plotting stuff %%

m_2_cm=100; Tes_2_gau=1e4; % metre to centimeter and Tesla to gauss

if plot1d==1
  
x_i=ceil(M/2);y_i=x_i;
x_j=x_i+1;y_j=y_i+1;
z_i=ceil(N/2);z_j=z_i+1;
figure
subplot(3,1,1)
plot(z_line*m_2_cm,reshape(baz(x_i,y_i,:),N,1)*Tes_2_gau...
    ,z_line*m_2_cm,reshape(baz(x_i,y_j,:),N,1)*Tes_2_gau...
    ,z_line*m_2_cm,reshape(baz(x_j,y_i,:),N,1)*Tes_2_gau...
    ,z_line*m_2_cm,reshape(baz(x_j,y_j,:),N,1)*Tes_2_gau...
    )
xlabel('z (cm)');
ylabel('B_z Field (g)');
title('AntiHelmholtz along z-axis \pm \Delta xy');
subplot(3,1,2)
plot(y_line*m_2_cm,reshape(bay(x_i,:,z_i),M,1)*Tes_2_gau...
    ,y_line*m_2_cm,reshape(bay(x_j,:,z_j),M,1)*Tes_2_gau...
    )
xlabel('y (cm)');
ylabel('B_y Field (g)');
title('AntiHelmholtz along y-axis \pm \Delta xz')
subplot(3,1,3)
plot(x_line*m_2_cm,reshape(bax(:,y_i,z_i),M,1)*Tes_2_gau...
    ,x_line*m_2_cm,reshape(bax(:,y_j,z_j),M,1)*Tes_2_gau...
     )
xlabel('x (cm)');
ylabel('B_x Field (g)');
title('AntiHelmholtz along x-axis \pm \Delta yz')

end
%% Field gradients anonymous function %%
db_x_x = @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K) C*z./(2*alp2.^2.*bet.^3.*rho2.^2).*...
  ((a^4*(-g.*(3*z.^2 +a^2) + rho2.*(8*x.^2-y.^2))-...
   a^2*(rho2.^2.*(5*x.^2 +y.^2)-2*rho2.*z.^2.*(2*x.^2+y.^2)+3*z.^4.*g)-...
   r2.^2.*(2*x.^4 +g.*(y.^2 +z.^2))).*E +...
   (a^2*(g.*(a^2 +2*z.^2)-rho2.*(3*x.^2-2*y.^2))+...
   r2.*(2*x.^4 +g.*(y.^2 +z.^2))).*alp2.*K);

db_x_y = @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K) C*x.*y.*z./(2*alp2.^2.*bet.^3.*rho2.^2).*...
  ((3*a^4*(3*rho2-2*z.^2) - r2.^2.*(2*r2+rho2)-...
   2*a^6-2*a^2*(2*rho2.^2-rho2.*z.^2+3*z.^4)).*E +...
  (r2.*(2*r2+rho2)-a^2*(5*rho2-4*z.^2)+2*a^4).*alp2.*K);

db_x_z = @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K) C*x./(2*alp2.^2.*bet.^3.*rho2).*...
  (((rho2-a^2).^2.*(rho2+a^2)+...
   2*z.^2.*(a^4-6*a^2*rho2+rho2.^2)+z.^4.*(a^2+rho2)).*E -...
   ((rho2-a^2).^2+z.^2.*(rho2+a^2)).*alp2.*K);

db_y_x=db_x_y;

db_y_y = @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K) C*z./(2*alp2.^2.*bet.^3.*rho2.^2).*...
  ((a^4*(g.*(3*z.^2 +a^2) + rho2.*(8*y.^2-x.^2))-...
   a^2*(rho2.^2.*(5*y.^2 +x.^2)-2*rho2.*z.^2.*(2*y.^2+x.^2)-3*z.^4.*g)-...
   r2.^2.*(2*y.^4 -g.*(x.^2 +z.^2))).*E +...
   (a^2*(-g.*(a^2 +2*z.^2)-rho2.*(3*y.^2-2*x.^2))+...
   r2.*(2*y.^4 -g.*(y.^2 +z.^2))).*alp2.*K);
 
db_y_z= @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K) C*y./(2*alp2.^2.*bet.^3.*rho2).*...
  (((rho2-a^2).^2.*(rho2+a^2)+...
   2*z.^2.*(a^4-6*a^2*rho2+rho2.^2)+z.^4.*(a^2+rho2)).*E -...
   ((rho2-a^2).^2+z.^2.*(rho2+a^2)).*alp2.*K);
 
db_z_x= db_x_z;
db_z_y= db_y_z;
db_z_z= @(C,a,x,y,z,alp2,bet,rho2,r2,g,E,K)C*z./(2*alp2.^2.*bet.^3).*...
       ((6*a^2*(rho2-z.^2)-7*a^4 +r2.^2).*E+...
       alp2.*(a^2-r2).*K);
     
     
%% coil1
dbxx1= db_x_x(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbxy1= db_x_y(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbxz1= db_x_z(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbyx1= db_y_x(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbyy1= db_y_y(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbyz1= db_y_z(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbzx1= db_z_x(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbzy1= db_z_y(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbzz1= db_z_z(C,a,xx,yy,zz-d1...
  ,alp_c_sq_1,sqrt(bet_c_sq_1),rho_sq_1,r_sq_1,gamma_1,E1,K1);
dbxx1(~isfinite(dbxx1))=0;
dbxy1(~isfinite(dbxy1))=0;
dbxz1(~isfinite(dbxz1))=0;
dbyx1(~isfinite(dbyx1))=0;
dbyy1(~isfinite(dbyy1))=0;
dbyz1(~isfinite(dbyz1))=0;
dbzx1(~isfinite(dbzx1))=0;
dbzy1(~isfinite(dbzy1))=0;

%% coil2
dbxx2= db_x_x(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbxy2= db_x_y(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbxz2= db_x_z(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbyx2= db_y_x(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbyy2= db_y_y(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbyz2= db_y_z(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbzx2= db_z_x(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbzy2= db_z_y(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbzz2= db_z_z(C,a,-xx,yy,-(zz-d2)...
  ,alp_c_sq_2,sqrt(bet_c_sq_2),rho_sq_2,r_sq_2,gamma_2,E2,K2);
dbxx2(~isfinite(dbxx2))=0;
dbxy2(~isfinite(dbxy2))=0;
dbxz2(~isfinite(dbxz2))=0;
dbyx2(~isfinite(dbyx2))=0;
dbyy2(~isfinite(dbyy2))=0;
dbyz2(~isfinite(dbyz2))=0;
dbzx2(~isfinite(dbzx2))=0;
dbzy2(~isfinite(dbzy2))=0;

%% coil3
dbxx3= db_x_x(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbxy3= db_x_y(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbxz3= db_x_z(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbyx3= db_y_x(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbyy3= db_y_y(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbyz3= db_y_z(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbzx3= db_z_x(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbzy3= db_z_y(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbzz3= db_z_z(C3,a3,-zz,yy+axis_shift,xx-d3...
  ,alp_c_sq_3,sqrt(bet_c_sq_3),rho_sq_3,r_sq_3,gamma_3,E3,K3);
dbxy3(~isfinite(dbxy3))=0;
dbxz3(~isfinite(dbxz3))=0;
dbyx3(~isfinite(dbyx3))=0;
dbyy3(~isfinite(dbyy3))=0;
dbyz3(~isfinite(dbyz3))=0;
dbzx3(~isfinite(dbzx3))=0;
dbzy3(~isfinite(dbzy3))=0;
dbzz3(~isfinite(dbzz3))=0;

%% coil4
dbxx4= db_x_x(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbxy4= db_x_y(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbxz4= db_x_z(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbyx4= db_y_x(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbyy4= db_y_y(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbyz4= db_y_z(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbzx4= db_z_x(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbzy4= db_z_y(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbzz4= db_z_z(C3,a4,-zz,yy+axis_shift,xx-d4...
  ,alp_c_sq_4,sqrt(bet_c_sq_4),rho_sq_4,r_sq_4,gamma_4,E4,K4);
dbxy4(~isfinite(dbxy4))=0;
dbxz4(~isfinite(dbxz4))=0;
dbyx4(~isfinite(dbyx4))=0;
dbyy4(~isfinite(dbyy4))=0;
dbyz4(~isfinite(dbyz4))=0;
dbzx4(~isfinite(dbzx4))=0;
dbzy4(~isfinite(dbzy4))=0;
dbzz4(~isfinite(dbzz4))=0;
%% Add up the field gradient components in the correct axes %%
% x1,y1,z1 == xx,yy,zz ; x2,y2,z2 == -xx,yy,-zz ;
% x3,y3,z3 == -zz,yy,xx ; x4,y4,z4 == -zz,yy,xx

% Diagonal terms first since they are easier.
dbxx=dbxx1+dbxx2+dbzz3+dbzz4;
dbyy=dbyy1+dbyy2+dbyy3+dbyy4;
dbzz=dbzz1+dbzz2+dbxx3+dbxx4;
% Other terms.
dbxy=dbxy1-dbxy2+dbzy3+dbzy4;
dbxz=dbxz1+dbxz2-dbzx3-dbzx4;
dbyx=dbyx1-dbyx2+dbyz3+dbyz4;
dbyz=dbyz1-dbyz2-dbyx3-dbyx4;
dbzx=dbzx1+dbzx2-dbxz3-dbxz4;
dbzy=dbzy1-dbzy2-dbxy3-dbxy4;

%% Plotting stuff %%

Tes_p_m_2_gau_p_cm=1e4/1e2; % Tesla per metre to gauss per cm

if plot1d==1
  
figure
plot(z_line*m_2_cm,reshape(dbzz(x_i,y_i,:),N,1)*Tes_p_m_2_gau_p_cm)
xlabel('z (cm)');
ylabel('Field gradient dB_z/dz (G/cm)');
title(['Field gradient along x=' num2str(xx(x_i,y_i,1)) ',y='...
       num2str(yy(x_i,y_i,1))]);
     
figure
plot(x_line*m_2_cm,reshape(dbzz(:,y_i,z_i),N,1)*Tes_p_m_2_gau_p_cm)
xlabel('x (cm)');
ylabel('Field gradient dB_z/dz (G/cm)');
title(['Field gradient along z=' num2str(zz(1,y_i,z_i)) ',y='...
       num2str(yy(x_i,y_i,1))]);
end

if plot2d~=0
  %% Plot in xz plane
  z_plot_range=[-0.015 0.015];
  y_plot_range=[-0.015 0.015];
  x_plot_range=[-0.015 0.015];
  z_pl_idx=find(z_line>=z_plot_range(1)&z_line<=z_plot_range(2));
  y_pl_idx=find(y_line>=y_plot_range(1)&y_line<=y_plot_range(2));
  x_pl_idx=find(x_line>=x_plot_range(1)&x_line<=x_plot_range(2));
  [xz,zx]=ndgrid(x_line(x_pl_idx)*m_2_cm,z_line(z_pl_idx)*m_2_cm);
  [yz,zy]=ndgrid(y_line(y_pl_idx)*m_2_cm,z_line(z_pl_idx)*m_2_cm);
  x_i=ceil(M/2);x_j=x_i+1;
  y_i=ceil(M/2);y_j=y_i+1;
  z_i=ceil(N/2);
  bz_2d=reshape((baz(x_pl_idx,y_i,z_pl_idx)+...
                 baz(x_pl_idx,y_j,z_pl_idx))/2 ...
               ,size(xz))*Tes_2_gau;
  bx_2d=reshape((bax(x_pl_idx,y_i,z_pl_idx)+...
                 bax(x_pl_idx,y_j,z_pl_idx))/2 ...
               ,size(xz))*Tes_2_gau;
  by_2d=reshape((bay(x_i,y_pl_idx,z_pl_idx)+...
                 bay(x_j,y_pl_idx,z_pl_idx))/2 ...
               ,size(yz))*Tes_2_gau;
  bt_2d=reshape((bat(x_pl_idx,y_i,z_pl_idx)+...
                 bat(x_pl_idx,y_j,z_pl_idx))/2 ...
               ,size(xz))*Tes_2_gau;
  bzz_2d=reshape((dbzz(x_pl_idx,y_i,z_pl_idx)+...
                 dbzz(x_pl_idx,y_j,z_pl_idx))/2 ...
                 ,size(xz))*Tes_p_m_2_gau_p_cm;
%% Plot in xy plane
%   z_plot_range=[-0.015 0.015];
%   y_plot_range=[-0.015 0.015];
%   x_plot_range=[-0.015 0.015];
  z_plot_range=[-0.15 0.15];
  y_plot_range=[-0.15 0.15];
  x_plot_range=[-0.15 0.15];

  z_pl_idx=find(z_line>=z_plot_range(1)&z_line<=z_plot_range(2));
  y_pl_idx=find(y_line>=y_plot_range(1)&y_line<=y_plot_range(2));
  x_pl_idx=find(x_line>=x_plot_range(1)&x_line<=x_plot_range(2));
  [xy,yx]=ndgrid(x_line(x_pl_idx)*m_2_cm,z_line(y_pl_idx)*m_2_cm);
%   [yz,zy]=ndgrid(y_line(y_pl_idx)*m_2_cm,z_line(z_pl_idx)*m_2_cm);
  x_i=ceil(M/2);x_j=x_i+1;
  y_i=ceil(M/2);y_j=y_i+1;
  z_i=ceil(N/2);z_j=z_i+1;
  bz2_2d=reshape((baz(x_pl_idx,y_pl_idx,z_i)+...
                  baz(x_pl_idx,y_pl_idx,z_j))/2 ...
               ,size(xy))*Tes_2_gau;
  bx2_2d=reshape((bax(x_pl_idx,y_pl_idx,z_i)+...
                  bax(x_pl_idx,y_pl_idx,z_j))/2 ...
               ,size(xy))*Tes_2_gau;
  by2_2d=reshape((bay(x_pl_idx,y_pl_idx,z_i)+...
                  bay(x_pl_idx,y_pl_idx,z_j))/2 ...
               ,size(xy))*Tes_2_gau;
  bt2_2d=reshape((bat(x_pl_idx,y_pl_idx,z_i)+...
                 bat(x_pl_idx,y_pl_idx,z_j))/2 ...
               ,size(xy))*Tes_2_gau;
  bxx2_2d=reshape((dbxx(x_pl_idx,y_pl_idx,z_i)+...
                 dbxx(x_pl_idx,y_pl_idx,z_j))/2 ...
                 ,size(xy))*Tes_p_m_2_gau_p_cm;               
  if plot2d==1
  figure
  subplot(3,1,1)
  contour(xz,zx,bz_2d,50)
%   colormap('Copper');
  xlabel('x (cm)');
  ylabel('z (cm)');
  title('Bz at y=0');
  subplot(3,1,2)
  contour(xz,zx,bx_2d,50)
%   colormap('Copper');
  xlabel('x (cm)');
  ylabel('z (cm)');
  title('Bx at y=0');
  subplot(3,1,3)
  contour(yz,zy,by_2d,50)
%   colormap('Copper');
  xlabel('y (cm)');
  ylabel('z (cm)');
  title('By at x=0');
  elseif plot2d==2
  figure('units','pixel','Position',[680 100 560 420*2]);
  subplot(2,1,1)
  contour(xz,zx,bt_2d,50)
   xlabel('x (cm)');
  ylabel('z (cm)');
  title(['B_{total} at y=0 at I_3=' num2str(I3) ....
         'A and I_0=' num2str(I) 'A']);
  subplot(2,1,2)
  plot(xz(:,z_i),bzz_2d(:,z_i))
  xlabel('x (cm)');
  ylabel('dBz/dz (G/cm)');
  title('dBz/dz at y=0,z\approx 0');
  elseif plot2d==3
  figure('units','pixel','Position',[680 100 560 420*2]);
  subplot(2,1,1)
  contour(xy,yx,bt2_2d,50)
   xlabel('x (cm)');
  ylabel('y (cm)');
  title(['B_{total} at z=0 at I_3=' num2str(I3) ....
         'A and I_0=' num2str(I) 'A']);
  subplot(2,1,2)
  plot(xy(:,y_i),bxx2_2d(:,y_i))
  xlabel('x (cm)');
  ylabel('dBx/dx (G/cm)');
  title('dBx/dx at z=0,y\approx 0');
  if save_plot==1
  saveLooseFigure(gcf,['Figure\Ioffe_field_I3_' num2str(I3) 'A.fig']);
  saveLooseFigure(gcf,['Figure\Ioffe_field_I3_' num2str(I3) 'A.pdf']);
  end
  elseif plot3d==3
    figure('units','pixel','Position',[680 100 560 420*2]);
    quiver(xz,zx,bx_2d,bz_2d);
     xlabel('x (cm)');
  ylabel('z (cm)');
  
  end
  if plot_fg==1
    
    figure;
    contour(xz,zx,bzz_2d,100)
    xlabel('x (cm)');
    ylabel('z (cm)');
    title('dBz/dx at y=0');
    figure
    surf(xz,zx,bzz_2d)
             
  end
  
end
  
end