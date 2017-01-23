function [tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,vw,vs,M,depth,isplot)
% MOOR 2016 CUMCM Problem A - Optimal Design of Mooring System
% Reference: http://canuck.seos.uvic.ca/rkd/mooring/moordyn.php
% 
% USAGE: [tiltdrum,elevanch,xsbed,xbuoy,f] = moor(Lc,chain,v,M,isplot)
%
%      tilt = tilt angle of the drum
%      elev = angle between the seabed and the welded chain at the anchor
%     xsbed = length of the welded chain that reclined on the seabed
%     xbuoy = swimming range of the buoy
%         f = immersion ratio of the of the buoy
%        Lc = length ofthe welded chain
%         vw = speeds of wind 
%         vs = speeds of sea water at sea surface
%         M = mass of the heavy ball
%    isplot = is plot?
%
% Zhou Lvwen: zhou.lv.wen@gmail.com
% September 12, 2016

if nargin==0;Lc=22.05;chain=2;vw=12;vs=0;M=1200;depth=18;isplot=1; end
if nargin==6;isplot=0; end

g = 9.81;                 % acceleration of gravity [m/s^2]
rho = 1.025e3;            % density of seawater [kg/m^3]
rhoFe = 7.9e3;            % density of ferrum [kg/m^3]
cdwin = 0.625;            % drag coefficient in wind

[lc, mc, dc] = chainpara(chain);
nc = round(Lc/lc);        % number of chain units

%   [buoy, 4 steel tubes, cylindrical drum, a series of welded chain]
h = [   2,      ones(1,4),   1, lc*ones(1,nc)];  % m
d = [   2, 5e-2*ones(1,4), 0.3, dc*ones(1,nc)];  % m
m = [1000,   10*ones(1,4), 100, mc*ones(1,nc)];  % kg

phi = zeros(1,length(h));

Fb = pi*(d/2).^2.*h*rho*g - m*g;         % Turn masses/buoyancies into forces

fmin = 0; fmax = 1;
while fmax-fmin>1e-10
    
    f = (fmax+fmin)/2;                   % determined by f dichotomy 
    
    Fb(1) = f*pi*(d(1)/2).^2.*h(1) * rho * g - m(1)*g;
    
    Fw = cdwin * (1-f)*h(1).*d(1) * vw.^2;  % wind load
    
    Fs = waterload(vs, h, d, phi, depth, f);% water load
    
    % solve equilibrium equations to determine tilt angle
    phi = solvequileq(Fb, Fw, Fs, M, f);
        
    x = h.*sin(phi); z = h.*cos(phi);    % projected length
    
    zw = sum(z(2:end)) + h(1)*f;         % waterline
    
    if zw>depth; fmax = f; else; fmin = f; end
end

% integrate from the bottom to the top
x = cumsum([0 fliplr(x)]);  z = cumsum([0 fliplr(z)]);

tilt = phi(6)*180/pi;                    % tilt angle of the drum
elev = 90- phi(end)*180/pi;              % elevation angle of chain at anchor
xsbed = max(x(z<1e-10));                 % length of the chain on the seabed
xbuoy = x(end-1);                        % swimming range of the buoy

if isplot;plotmoor(x,z,h,phi,Lc,chain,vw,vs,M,depth,tilt,elev,xbuoy,f);end

% ------------------------------------------------------------------------

function Fs = waterload(vs, h, d, phi, depth, f)
cd = 374;                 % drag coefficient in sea water 
z = h.*cos(phi);
zi = fliplr(cumsum(fliplr(z))) - h.*cos(phi)/2; 
zi(1) = depth - f*h(1)/2;
vsi = vs./sqrt(depth)*sqrt(zi);          % water volecity proflile
Fs = cd * h.*d.*cos(phi) .* vsi.^2;      % water load
Fs(1) = Fs(1)*f; 

% ------------------------------------------------------------------------

function phi = solvequileq(Fb, Fw, Fs, M, f)
g = 9.81;                 % acceleration of gravity [m/s^2]
rho = 1.025e3;            % density of seawater [kg/m^3]
rhoFe = 7.9e3;            % density of ferrum [kg/m^3]
N = length(Fb);
[theta, phi, Ft] = deal(zeros(1,N));
for i = 1:N-1
    fx = Ft(i)*sin(theta(i)) + Fs(i);
    if i==1; fx = fx + Fw; end
    
    fz = Fb(i) + Ft(i)*cos(theta(i));
    % add load from heavy ball = weight - buoyancy
    if i==6; fz = fz - M*g + rho*(M/rhoFe)*g; end
    
    Ft(i+1) = sqrt(fx^2+fz^2);
    
    theta(i+1) = acos(fz/Ft(i+1));   % the tilt of Ft
    
    if theta(i+1)>pi/2; theta(i+1) = pi/2; end
end

phi =atan2( Ft.*sin(theta)+Fs/2, Ft.*cos(theta)+Fb/2);
phi(phi>pi/2) = pi/2;                % the tilt of each element
phi(1) = atan2( Fs(1)*f/2+Fw(1)*(f+(1-f)/2), Fb(1)*f/2 );

% ------------------------------------------------------------------------

function [lc, mc, dc] = chainpara(typeid)
% The types and parameters of the welded chains

rhoFe = 7.9e3;  % density of ferrum [kg/m^3];
rho= [3.2 7.0 12.5 19.5 28.12];
% length(lc) and mass(mc) of one unit of the welded chain.
lc = [ 78 105  120  150   180]*1e-3; % m    
mc = rho.*lc;                        % kg

lc = lc(typeid);
mc = mc(typeid);
% one unit of the welded chain was assumed as a cylinder.
dc = 2*sqrt(rho(typeid)/rhoFe/pi);   % diameter of the welded chain
% -------------------------------------------------------------------------


function plotmoor(x,z,h,phi,Lc,chain,vw,vs,M,depth,tilt,elev,xbuoy,f)
figure; N = length(x);
% draw buoy
x1 = x(N)  -cos(phi(1));   x2 = x(N)  +cos(phi(1)); 
x3 = x(N-1)+cos(phi(1));   x4 = x(N-1)-cos(phi(1));
z1 = z(N)  +sin(phi(1));   z2 = z(N)  -sin(phi(1));   
z3 = z(N-1)-sin(phi(1));   z4 = z(N-1)+sin(phi(1));
fill([x1 x2 x3 x4], [z1 z2 z3 z4], 'r');hold on
text(x(N-1)+1.3, z(N-1), sprintf('%2.1f m',f*h(1)));

% draw steel tubes
plot(x(N-5:N-1),z(N-5:N-1),'o-','linewidth',1, 'markersize',2);

% draw the drum
x1 = x(N-5)-0.15*cos(phi(6));   x2 = x(N-5)+0.15*cos(phi(6)); 
x3 = x(N-6)+0.15*cos(phi(6));   x4 = x(N-6)-0.15*cos(phi(6));
z1 = z(N-5)+0.15*sin(phi(6));   z2 = z(N-5)-0.15*sin(phi(6));   
z3 = z(N-6)-0.15*sin(phi(6));   z4 = z(N-6)+0.15*sin(phi(6));
fill([x1 x2 x3 x4], [z1 z2 z3 z4], 'r');hold on
text(x(N-6)-2,z(N-6)+0.5,sprintf('%3.1f^\\circ',tilt));

% draw the welded chain
plot(x(1:N-6),z(1:N-6),'-m')
text(0.5,1.5,sprintf('%3.1f^\\circ',elev));
types = {'I','II','III','IV','V'};
text(x(ceil(N/2))+1, z(ceil(N/2)), sprintf('%s:%5.1f m',types{chain}, Lc));

% draw the heavy ball 
plot([x(N-6) x(N-6)],[z(N-6) z(N-6)-4],'-r');
plot(x(N-6),z(N-6)-4,'.r','markersize',50)
text(x(N-6)-1,z(N-6)-5,sprintf('%d kg',M));

axismax = max(depth+2,xbuoy+4);
% draw waterline
plot([0,axismax],[depth,depth],'--');

% wind velocity
quiver(x(N)-8,depth+1,2,0,'MaxHeadSize',0.5,'color','r');
text(x(N)-5.5,depth+1,sprintf('%3.1f m/s', vw),'color','r');

% water velocity
quiver(x(N)-8,depth-1,2,0,'MaxHeadSize',0.5,'color','b');
text(x(N)-5.5,depth-1,sprintf('%3.1f m/s', vs),'color','b');

box on; axis image; 
set(gca,'xtick',[max(x(z<1e-10)),round(xbuoy*10)/10],'xgrid','on')
axis([0,axismax,0,axismax]); 
xlabel('x'); ylabel('z')
