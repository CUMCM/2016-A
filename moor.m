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
cdsea = 374;              % drag coefficient in sea water 

[lc, mc, dc] = chainpara(chain);
nc = round(Lc/lc);        % number of chain units

%   [buoy, 4 steel tubes, cylindrical drum, a series of welded chain]
h = [   2,      ones(1,4),   1, lc*ones(1,nc)];  % m
m = [1000,   10*ones(1,4), 100, mc*ones(1,nc)];  % kg
d = [   2, 5e-2*ones(1,4), 0.3, dc*ones(1,nc)];  % m

N = length(h);
[Ft, Fw, theta, phi] = deal(zeros(1,N));

Fb = pi*(d/2).^2.*h*rho*g - m*g;         % Turn masses/buoyancies into forces
vsi = @(zi)(vs./sqrt(18)*sqrt(zi));      % water volecity proflile

zi = fliplr(cumsum(fliplr(h.*cos(phi)))); zi = zi/max(zi)*depth;

fmin = 0; fmax = 1;
while fmax-fmin>1e-10
    
    f = (fmax+fmin)/2;                   % determined by f dichotomy 
    
    Fb(1) = pi*(d(1)/2).^2.*h(1)*rho*g*f-m(1)*g;
    
    S = h.*d;                            % exposed area of cylinder/wire
    Fw(1) = cdwin*(1-f)*S(1)*vw.^2;      % wind load
    
    zi(1) = depth-f*h(1)/2;
    Fs = cdsea.*S.*cos(phi).*vsi(zi).^2; % water load
    
    Ft(1) = 0;
    theta(1) = 0;   
    
    for i = 1:N-1
        fx = Ft(i)*sin(theta(i)) + Fw(i) + Fs(i);
        
        fz = Fb(i) + Ft(i)*cos(theta(i));
        % add load from heavy ball = weight - buoyancy
        if i==6; fz = fz - (M-M/7.9)*g; end
        
        Ft(i+1) = sqrt(fx^2+fz^2);
        
        theta(i+1) = acos(fz/Ft(i+1));   % the tilt of Ft
        if theta(i+1)>=pi/2; theta(i+1) = pi/2; end
    end
    
    phi =atan2( Ft.*sin(theta)+Fs(i)/2, Ft.*cos(theta)-Fb/2 );
    phi(theta==pi/2)=pi/2; phi(1)=0;     % the tilt of each element
    
    x = h.*sin(phi); z = h.*cos(phi);    % projected length
    
    zi = fliplr(cumsum(fliplr(z))) - h.*cos(phi)/2; 
    zw = sum(z)-h(1)*(1-f);              % waterline
    
    if zw>depth; fmax = f; else; fmin=f; end
    
end

% integrate from the bottom to the top
x = cumsum([0 fliplr(x)]);
z = cumsum([0 fliplr(z)]);

tilt = phi(6)*180/pi;                    % tilt angle of the drum
elev = 90- phi(end)*180/pi;              % elevation angle of chain at anchor
xsbed = max(x(z<1e-10));                 % length of the chain on the seabed
xbuoy = max(x);                          % swimming range of the buoy

if isplot; plotmoor(x,z,Lc,chain,vw,M,depth,tilt,elev); end
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


function plotmoor(x,z,Lc,chain,v,M,depth,tilt,elev)
figure('name',sprintf('Lc=%4.2f m, chain=%d, v=%3.1f m/s, M=%4.1f kg',...,
                       Lc,         chain,    v,           M));
x1 = x(end)-1; x2 = x1+2; x3 = x(end-1)+1; x4 = x3-2;
z1 = z(end);   z2 = z1;   z3 = z(end-1);   z4 = z3;
fill([x1 x2 x3 x4], [z1 z2 z3 z4], 'r')
hold on
plot(x(end-6:end-1),z(end-6:end-1),'o-','linewidth',2, 'markersize',3);
plot(x(1:end-6),z(1:end-6),'-m','linewidth',1)

plot([x(end-6) x(end-6)],[z(end-6) z(end-6)-2],'-r','linewidth',1);
plot(x(end-6),z(end-6)-2,'.r','linewidth',1,'markersize',50)

plot([0,depth+2],[depth,depth],'--')
text(1,2,sprintf('%3.1f degree',elev));
text(x(end-6)-5,z(end-6)+0.5,sprintf('%3.1f degree',tilt));
box on; axis image; grid on
set(gca,'xtick',[max(x(z<1e-10)),round(x(end)*10)/10])
axis([0,depth+2,0,depth+2]); xlabel('x'); ylabel('z')
title(sprintf('Lc=%4.2f m, chain=%d, v=%3.1f m/s, M=%4.1f kg',...,
                       Lc,         chain,    v,           M));
