% 2016 CUMCM problem A - Optimal Design of Mooring System
%
% Zhou Lvwen: zhou.lv.wen@gmail.com

Lc = 22.05;               % m
chain = 2;
M = 1200;                 % kg;
depth = 18;               % m

% -------------------------------------------------------------------------
% Question 1
%
[tilt,elev,xsbed,xbuoy,f] = moor(Lc, chain, 12, 0, M, depth, 1);
[tilt,elev,xsbed,xbuoy,f] = moor(Lc, chain, 24, 0, M, depth, 1);

% -------------------------------------------------------------------------
% Question 2
%
i = 0;
Mi = [];
tilti = [];
elevi = [];
while f<0.7
    i = i + 1;
    mi = 1200+10*i;
    [tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,36,0,mi, depth);
    Mi = [Mi, mi];
    tilti = [tilti, tilt];
    elevi = [elevi, elev];
end

figure('name','problem 2'); 
plot(Mi, tilti,'r', Mi, elevi, 'b');

xi = [min(Mi), max(Mi)];
hold on; plot(xi,[5,5],'r--',xi,[16,16],'b--');  
set(gca,'ytick',[0,5,10,16,20]); grid on
legend('Tilt angle of the drum','Elevation angle of the chain at the anchor')
xlabel('Mass of the heavy ball (kg)'); ylabel('Angle (degree)')

% -------------------------------------------------------------------------
% Question 3
Lc = 22.05;               % m
chain = 2;
vw = 36;                  % m/s
vs = 1.5;                 % m/s
M = 3000;                 % kg
depth = 20;               % m
[tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,vw,vs,M,20,1);