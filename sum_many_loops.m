function [B,JB]= sum_many_loops(fun,a,I,xx,yy,zz,x0,y0,z0)
% This program calculates the sum of many loops of coils pointing in the
% same direction. fun can either be num_circ_coil or sym_rect_coil(square
% with a==b). a and one axis separation, either x0,y0 or z0, should be an
% array with same size created by ndgrid
no_a=numel(a);
x_f=(numel(x0)==no_a);
y_f=(numel(y0)==no_a);
z_f=(numel(z0)==no_a);
if ~x_f
  x0=x0*ones(size(a));
end
if ~y_f
  y0=y0*ones(size(a));
end
if ~z_f
  z0=z0*ones(size(a));
end
B0=zeros([size(xx) 3 no_a]);
JB0=zeros([size(xx) 3 3 no_a]);
if strcmp(fun,'sym_rect_coil')
    b=a;
    for i=1:no_a
      [B0(:,:,:,:,i),JB0(:,:,:,:,:,i)]=sym_rect_coil(a(i),b(i),I,xx,yy,zz,x0(i),y0(i),z0(i));
    end
elseif strcmp(fun,'num_rect_coil')
    b=a;
    for i=1:no_a
      [B0(:,:,:,:,i),JB0(:,:,:,:,:,i)]=num_rect_coil(a(i),b(i),I,xx,yy,zz,x0(i),y0(i),z0(i));
    end
elseif strcmp(fun,'num_circ_coil')
  for i=1:no_a
      [B0(:,:,:,:,i),JB0(:,:,:,:,:,i)]=num_circ_coil(a(i),I,xx,yy,zz,x0(i),y0(i),z0(i));
      
  end
end
B=sum(B0,5);
JB=sum(JB0,6);

end


