% This script computes the field from thin loops of coils.
% Magnetic fields for Square and Circular coils are computed using
% functions sym_rect_coil num_circ_coil respectively.
% Computes  the  a field  from   thin  loops  of  radius  a,
% located  at  z  = +|d2|  and  -|d1|, in anti-helmholtz configuration
% current  I.
% and a pair helmholtz coil located at x=+|d3| and -|d3|
% with radius a3 and current I3
clear all
% close all

plot1d=0;
plot2d=2;
plot_fg=0;
save_plot=0;
%% Define coil parameters %%
layer1 = 6;             % Number  of  turns  in  the  loops
player1 = 8;
I = 20;      % Amperes in coils 1 and 2
% a = 0.0502;         % Radius 1 and 2 (Metres) 
% d=(0.067+2*4*0.00208)/2;% coil position.

a=0.045:2.08e-3*cos(pi/6):54.01e-3;
d=0.067/2:2.08e-3:0.0482;
d1 = -d;            % Coil 1 at z= - d; Current loop point towards z.
d2 = d;             % Coil 2 at z= + d; Current loop point towards -z.

I3=5;             % Amperes in coil 3
% dI=((0.083+2*3*0.00208)/2)-15e-3;
% aI=(0.033+0.033+2*5*cos(pi/6)*0.00208)/2;

% dI=(0.083/2:2.08e-3:0.052); %% Coils from 41.5 mm to 51.9mm. ~93cm separation
dI_sep = [83e-3+1.04e-3,repmat(2.08e-3*2,1,5)];
dI = cumsum(dI_sep);
% dI=83e-3+1.04e-3,83e-3+1.04e-3+2.08e-3*6,6); %% Coils from 41.5 mm to 51.9mm. ~93cm separation
aI=linspace(0.043-2.08e-3*2,0.043+2.08e-3*2,5);   %% Coils with sides 33.56mm to 48mm.
a3 = aI;        % Sides 3
d3 = dI;         % Coil 3 at x= |d3|; Current loop point towards x.
a4 = aI;        % Sides 4
d4 = -dI;         % Coil 4 at x= -|d4|; Current loop point towards x.


%% Define calculation parameters %%
% Define  the  plane  over  which  fields  are  computed
% N must  be  odd  to  include  the  point  (0,0)
% Use even N to avoid divide by zeros.
M=40;  % No.  of  points  along  the  x,y  axis
N=40;  % No.  of  points  along  the  z  axis
xr=0.02;
x_range=[-xr xr]; % in meters
y_range=[-xr xr];
z_range=[-xr xr];

%% Generate Axes for calculations %%
x_line=linspace(x_range(1),x_range(2),M);  
y_line=linspace(y_range(1),y_range(2),M);  
z_line=linspace(z_range(1),z_range(2),N);
[xx,yy,zz] = ndgrid(x_line,y_line,z_line);

%%
%%%%%%%%%%%%%%%%%% Generate loops position and size %%%%%%%%%%%%%%%%%%%%%%%
[aa1, dd1]=ndgrid(a,d1);
[aa2, dd2]=ndgrid(a,d2);
[aa3, dd3]=ndgrid(a3,d3);
[aa4, dd4]=ndgrid(a3,d4);

%%%%%%%%%%%%%%%%%% Calculate Here %%%%%%%%%%%%%%%%%%%%%%%
tic
[B1,JB1]=sum_many_loops('num_circ_coil',aa1,I,xx,yy,zz,0,0,dd1);
[B2,JB2]=sum_many_loops('num_circ_coil',aa2,I,-xx,yy,-zz,0,0,-dd2);
[B3,JB3]=sum_many_loops('num_rect_coil',aa3/2,I3,-zz,yy,xx,0,0,dd3);
[B4,JB4]=sum_many_loops('num_rect_coil',aa4/2,I3,-zz,yy,xx,0,0,dd4);
toc
% [B1,JB1]=num_circ_coil(a,I,xx,yy,zz,0,0,d1);
% [B2,JB2]=num_circ_coil(a,I,-xx,yy,-zz,0,0,-d2);
% [B3,JB3]=num_circ_coil(a3,I3,-zz,yy,xx,0,0,d3);
% [B4,JB4]=num_circ_coil(a3,I3,-zz,yy,xx,0,0,d4);

%%% Add up the field components in the correct axes %%
% x1,y1,z1 == xx,yy,zz ; x2,y2,z2 == -xx,yy,-zz ;
% x3,y3,z3 == -zz,yy,xx ; x4,y4,z4 == -zz,yy,xx

Bsum(:,:,:,1)=B1(:,:,:,1)-B2(:,:,:,1)+B3(:,:,:,3)+B4(:,:,:,3);
Bsum(:,:,:,2)=B1(:,:,:,2)+B2(:,:,:,2)+B3(:,:,:,2)+B4(:,:,:,2);
Bsum(:,:,:,3)=B1(:,:,:,3)-B2(:,:,:,3)-B3(:,:,:,1)-B4(:,:,:,1);
baz=Bsum(:,:,:,3);bay=Bsum(:,:,:,2);bax=Bsum(:,:,:,1);
bat(:,:,:)=sqrt(Bsum(:,:,:,1).^2 +Bsum(:,:,:,2).^2+Bsum(:,:,:,3).^2);

%%% Add up the field gradient components in the correct axes %%
% x1,y1,z1 == xx,yy,zz ; x2,y2,z2 == -xx,yy,-zz ;
% x3,y3,z3 == -zz,yy,xx ; x4,y4,z4 == -zz,yy,xx

% Diagonal terms first since they are easier.
dbxx=JB1(:,:,:,1,1)+JB2(:,:,:,1,1)+JB3(:,:,:,3,3)+JB4(:,:,:,3,3);
dbyy=JB1(:,:,:,2,2)+JB2(:,:,:,2,2)+JB3(:,:,:,2,2)+JB4(:,:,:,2,2);
dbzz=JB1(:,:,:,3,3)+JB2(:,:,:,3,3)+JB3(:,:,:,1,1)+JB4(:,:,:,1,1);
% Other terms.
dbxy=JB1(:,:,:,1,2)-JB2(:,:,:,1,2)+JB3(:,:,:,3,2)+JB4(:,:,:,3,2);
dbxz=JB1(:,:,:,1,3)+JB2(:,:,:,1,3)-JB3(:,:,:,3,1)-JB4(:,:,:,3,1);
dbyx=JB1(:,:,:,2,1)-JB2(:,:,:,2,1)+JB3(:,:,:,2,3)+JB4(:,:,:,2,3);
dbyz=JB1(:,:,:,2,3)-JB2(:,:,:,2,3)-JB3(:,:,:,2,1)-JB4(:,:,:,2,1);
dbzx=JB1(:,:,:,3,1)+JB2(:,:,:,3,1)-JB3(:,:,:,1,3)-JB4(:,:,:,1,3);
dbzy=JB1(:,:,:,3,2)-JB2(:,:,:,3,2)-JB3(:,:,:,1,2)-JB4(:,:,:,1,2);

%% Plotting stuff %%

m_2_cm=100; Tes_2_gau=1e4; % metre to centimeter and Tesla to gauss
Tes_p_m_2_gau_p_cm=1e4/1e2; % Tesla per metre to gauss per cm

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
  plot(zx(x_i,:),bzz_2d(x_i,:))
  xlabel('z (cm)');
  ylabel('dBz/dz (G/cm)');
  title('dBz/dz at y=0,x\approx 0');
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
  
% end

  