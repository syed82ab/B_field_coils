function [B_ret,JacobianB_ret]=num_rect_coil(a,b,I,xx0,yy0,zz0,x0,y0,z0)
% This function calculates the magnetic field and its Jacobian at xx,yy,zz,
% for a rectangular coil at position x0,y0,z0 with its axis in +z.
% The length of the side along x direction (perpendicular to y) is 2b0 and 
% along y direction (perp to x) is 2a0.
%
%function [B,JacobianB]=num_rect_coil(a0,b0,xx,yy,zz,x0,y0,z0)
%size(B)= [size(xx) ,3 ]
%size(JB)= [size(xx, 3, 3] ;

% Reference from Appendix of
% Dejana Herceg, Anamarija Juhas, and Miodrag Milutinov
% http://facta.junis.ni.ac.rs/eae/fu2k93/2herceg.pdf

mu_0=4*pi*1e-7; %Tesla
% syms a b 'positive'
% syms x y z 'real'
% r=sym('r',[1 4]);
% r=sym(r,'positive');
Bx=zeros(size(xx0));
By=zeros(size(yy0));
Bz=zeros(size(zz0));
xx=repmat(xx0-x0,[1,1,1,4]);
yy=repmat(yy0-y0,[1,1,1,4]);
zz=repmat(zz0-z0,[1,1,1,4]);
c_factor=[a -a -a a];
d_factor=[b b -b -b];
c=zeros(size(xx));
d=zeros(size(yy));
for i=1:4
c(:,:,:,i)=xx(:,:,:,1)+ c_factor(i);
d(:,:,:,i)=yy(:,:,:,1)+ d_factor(i);
r=sqrt(c.^2+d.^2 +zz.^2);
end
exp_fac=(1:4)+1;

const=mu_0*I/4/pi;

for i=1:4
Bx=(-1).^exp_fac(i).*zz(:,:,:,i)./r(:,:,:,i)./(r(:,:,:,i)+d(:,:,:,i)) + Bx;
By=(-1).^exp_fac(i).*zz(:,:,:,i)./r(:,:,:,i)./(r(:,:,:,i)+c(:,:,:,i)) + By;
Bz=(-1).^(exp_fac(i)-1).*(c(:,:,:,i)./r(:,:,:,i)./(r(:,:,:,i)+d(:,:,:,i))...
                    +d(:,:,:,i)./r(:,:,:,i)./(r(:,:,:,i)+c(:,:,:,i))) + Bz;

end
% B=[Bx ; By ;Bz]*const;


B_ret(:,:,:,1)=Bx*const;     %%Bx component
B_ret(:,:,:,2)=By*const;   %%By component
B_ret(:,:,:,3)=Bz*const; %%Bz component

s_r=r;
A1=-(c./s_r.*(r+d) + r.*c./s_r);
B1=-(d./s_r.*(r+d) + r.*(d./s_r +1));
G1= r.*(r+d);
G2= G1.*G1;
C1=(-zz.*((r+d)./r + 1));

D1=-(c./s_r.*(r+c) + r.*(c./s_r +1));
E1=-(d./s_r.*(r+c) + r.*d./s_r);
H1= r.*(r+c);
H2= H1.*H1;
F1=(-zz.*((r+c)./r + 1));

JBxX=zeros(size(xx0));
JBxY=zeros(size(yy0));
JBxZ=zeros(size(zz0));
JByX=zeros(size(xx0));
JByY=zeros(size(yy0));
JByZ=zeros(size(zz0));
JBzX=zeros(size(xx0));
JBzY=zeros(size(yy0));
JBzZ=zeros(size(zz0));
% JBzZ1=zeros(size(zz0));

for i=1:4
JBxX=(-1).^exp_fac(i).*A1(:,:,:,i).*zz(:,:,:,i)./G2(:,:,:,i) + JBxX;
JBxY=(-1).^exp_fac(i).*B1(:,:,:,i).*zz(:,:,:,i)./G2(:,:,:,i) + JBxY;
JBxZ=(-1).^exp_fac(i).*(G1(:,:,:,i)+C1(:,:,:,i).*zz(:,:,:,i))./G2(:,:,:,i) + JBxZ;
JByX=(-1).^exp_fac(i).*D1(:,:,:,i).*zz(:,:,:,i)./H2(:,:,:,i) + JByX;
JByY=(-1).^exp_fac(i).*E1(:,:,:,i).*zz(:,:,:,i)./H2(:,:,:,i) + JByY;
JByZ=(-1).^exp_fac(i).*(H1(:,:,:,i)+F1(:,:,:,i).*zz(:,:,:,i))./H2(:,:,:,i) + JByZ;
JBzX=(-1).^(exp_fac(i)-1).*((G1(:,:,:,i)+A1(:,:,:,i).*c(:,:,:,i))./G2(:,:,:,i) + ...
                            D1(:,:,:,i).*d(:,:,:,i)./H2(:,:,:,i)) + JBzX;
JBzY=(-1).^(exp_fac(i)-1).*(B1(:,:,:,i).*c(:,:,:,i)./G2(:,:,:,i) + ...
                           (H1(:,:,:,i)+E1(:,:,:,i).*d(:,:,:,i))./H2(:,:,:,i)) + JBzY;
JBzZ=(-1).^(exp_fac(i)-1).*(C1(:,:,:,i).*c(:,:,:,i)./G2(:,:,:,i) + ...
                            F1(:,:,:,i).*d(:,:,:,i)./H2(:,:,:,i)) + JBzZ;
end


JacobianB_ret(:,:,:,1,1)=JBxX*const;
JacobianB_ret(:,:,:,1,2)=JBxY*const;
JacobianB_ret(:,:,:,1,3)=JBxZ*const;

JacobianB_ret(:,:,:,2,1)=JByX*const;
JacobianB_ret(:,:,:,2,2)=JByY*const;
JacobianB_ret(:,:,:,2,3)=JByZ*const;

JacobianB_ret(:,:,:,3,1)=JBzX*const;
JacobianB_ret(:,:,:,3,2)=JBzY*const;
JacobianB_ret(:,:,:,3,3)=JBzZ*const;
