% from http://dtic.mil/dtic/tr/fulltext/u2/a286081.pdf

% Like  loops.m  but  symetric  plane  through  z  axis
% Computes  the  a field  from  2  thin  loops  of  radius  a,
% located  at  z  = +d  and  -d,  (d-a/2,  felaboltz  Coils),  and
% current  I.

LN = 48;             % Number  of  turns  in  the  loops
mu_0 = 4*pi*1e-7;   % Telsa
I = 0.019052; %1.0;              % Amperes
a = 0.6096;          % Radius (Metres)
% d = 0.60;           % Coils at z= ± d;
d=a/2;

dend=d-0.1;
pend=a-0.05;

% Define  the  plane  over  which  fields  are  computed
% N must  be  odd  to  include  the  point  (0,0)
M=26;  % No.  of  points  along  the  rho  axis
N=101;  % No.  of  points  along  the  z  axis
P1=linspace(0,pend,M);  % rho  always  positive
P =[  P1(M:-1:2) P1];
z=linspace(-dend,dend,N);
[p1,z1]=meshgrid(P,z);
[p,z] = meshgrid(P1,z);  % rho=O,  p(:,1);  z= -d,  z(1,:)

% Determine  modulus  and  E1liptic  integrals
k1=(4*a*p)./( (a+p).^2 + (z-d).^2 );
[K1,E1] = ellipke(k1);
k2=(4*a*p)./( (a+p).^2 + (z+d).^2 );
[K2,E2] = ellipke(k2);

% Compute  B-rho fields
bpl = ((z-d)./(p.*sqrt( (a+p).^2 + (z-d).^2)))  .* ...
( ( (a^2  +  p.^2  +  (z-d).^2)./((a-p).^2 + (z-d).^2 ) ).*E1  - K1 );
bp2 = ((z+d)./(p.*sqrt( (a+p).^2 + (z+d).^2)))  .* ...
( ( (a^2  +  p.^2  +  (z+d).^2)./((a-p).^2 + (z+d).^2 ) ).*E2  - K2  );
bp  = ( (LN*mu_0*I)/(2*pi) )*(bpl + bp2);
bp(1:N) = zeros(N,1);  % Remove  NaN's  from  rho=0 axis
bp = [ bp(:,M:-1:2) bp ];  % Make  it  mymetric

bap  = ( (LN*mu_0*I)/(2*pi) )*(bpl - bp2);
bap(1:N) = zeros(N,1);  % Remove  NaN's  from  rho=0 axis
bap = [ bap(:,M:-1:2) bap ];  % Make  it  mymetric

% Compute  B-z fields
bzl = ( 1 ./sqrt( (a+p).^2 + (z-d).^2 ) ) .* ...
( ((a^2 - p.^2 - (z-d).^2)./((a-p).^2 + (z-d).^2)).*E1 + K1 );
bz2 = ( 1 ./sqrt( (a+p).^2 + (z+d).^2 ) ) .* ...
( ((a^2 - p.^2 - (z+d).^2)./((a-p).^2 + (z+d).^2)).*E2 + K2 );
bz = ( (LN*mu_0*I)/(2*pi)  )*(bzl  +  bz2);
bz(~isfinite(bz))=0;
bz = [ bz(:,M:-1:2) bz ];  % Make  it  symmetric

baz = ( (LN*mu_0*I)/(2*pi)  )*(bzl  -  bz2);
baz = [ baz(:,M:-1:2) baz ];  % Make  it  symmetric


%  Compute  the  total  B  field
bt = sqrt(bp.^2  +  bz.^2);

bat = sqrt(bap.^2  +  baz.^2);
% test  along  the  z  axis  where  we  know  the  solution
bzz = 1  ./(sqrt(a^2  +  (z(:,ceil(M/2))-d).^2)).^3  +  1 ./(sqrt(a^2  +  (z(:,ceil(M/2))+d).^2)).^3;
bzz = 0.5*LN*mu_0*I*a^2*bzz;

bazz= 1  ./(sqrt(a^2  +  (z(:,ceil(M/2))-d).^2)).^3  -  1 ./(sqrt(a^2  +  (z(:,ceil(M/2))+d).^2)).^3;
bazz=0.5*LN*mu_0*I*a^2*bazz;
figure
subplot(2,1,1)
plot(z(:,ceil(M/2)),bzz)
xlabel('z');
ylabel('B Field (T)');
title('Helmholtz 1D equation along axis');
subplot(2,1,2)
plot(z(:,ceil(size(bz,2)/2)),bz(:,ceil(size(bz,2)/2)))
xlabel('z');
ylabel('B Field (T)');
title('Helmholtz Complete equation, \rho=0')


