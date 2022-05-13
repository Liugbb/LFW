clc
clear
%% % This is a matlab script that generates the input data
% Dimensions of grid
a=5;%%km
nx=20*10;      ny=29*14;          nz=100;
Ho=2000;  % ocean depth in meters

% Vertical resolution (m)
dz=zeros(nz,1);    zc=zeros(nz,1);
dz(1) = 20;        zc(1) = -10;
for i=2:nz
    dz(i)=20;         zc(i)=zc(i-1)-20;
end
fid=fopen('delZvar','w','b'); fwrite(fid,dz,'real*8'); fclose(fid);

% Flat bottom at z=-Ho
h=-Ho*ones(nx,ny);

% Walls (surrounding domain) - generate bathymetry file
h(:,[1 end])=0;
fid=fopen('bathy.bin','w','b'); fwrite(fid,h,'real*8'); fclose(fid);


%% wind_force
dx=[];xc=[];
dy=[];yc=[];
dx=zeros(nx,1);  xc=zeros(nx,1);
dx(1) =    1000*a;  xc(1) = 500*a;
for i=2:nx/2
    dx(i) = 1000*a;      xc(i) = xc(i-1)+1000*a;
end
dy=zeros(ny,1);  yc1=zeros(ny,1);
dy(1) =    1000*a;  yc(1) = 500*a;
for i=2:ny
    dy(i) = 1000*a;      yc(i) = yc(i-1)+1000*a;
end
%%%Lower Frequency Wind-stress
taumax=1*0.1;%wind stress maximum
tau_s=zeros(size(h));
[X,Y]=ndgrid(xc,yc);
tau_s(:,:)=taumax*sin(pi*Y/(ny*dy(1))).*sin(pi*Y/(ny*dy(1)));%% wind stress profile
% % % tau_s(:,47:86)=tau_1;
fid=fopen('windx_siny.bin','w','b'); fwrite(fid,tau_s,'real*8'); fclose(fid);

% % % %% restoring temperature
% % % % 2-D SST field for relaxation
% % % Tmax=10.0; Tmin= 2.0; y=0.1:ny;
% % % % % T_surf=repmat(Tmin+(Tmax-Tmin)*y/(ny),[nx 1]);
% % % T_surf=repmat(Tmin+(Tmax-Tmin)*y(ny/2)/(ny),[nx ny]);
% % % fid=fopen('SST_relax.bin','w','b');fwrite(fid,T_surf,'real*8');fclose(fid);
% % % %% Temeture file
% % % load Tempeture
% % % Tref=Tempeture;
% % % t=0.0*rand([nx,ny,nz]);
% % % for k=1:nz
% % % t(:,:,k) = t(:,:,k) + Tref(k);
% % % end
% % % fid=fopen('T.init.bin','w','b'); fwrite(fid,t,'real*8'); fclose(fid);
% % % fid=fopen('tRef','w','b'); fwrite(fid,Tref,'real*8'); fclose(fid);

%% restoring temperature
% 2-D SST field for relaxation
Tmax=10.0; Tmin= -2.0; y=0.05:0.1:ny/10;
T_surf=repmat(Tmin+(Tmax-Tmin)*y/(ny/10),[nx 1]);
fid=fopen('SST_relax.bin','w','b');fwrite(fid,T_surf,'real*8');fclose(fid);
%% temprature
% 3-D Temperature field for initial conditions and RBCS relaxation
hh=500; %e-folding scale for temperature decrease with depth
T_3D=zeros(nx,ny,nz);
for k=1:length(dz)
    T_3D(:,:,k)=(T_surf - Tmin)*(exp(zc(k)/Ho) - exp(-Ho/hh))/(1-exp(-Ho/hh)) + Tmin;
end
fid=fopen('hydrogThetaFile.bin','w','b');fwrite(fid,T_3D,'real*8');fclose(fid);
% % TT=T_3D(:,:,1);