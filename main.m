% 2016 CUMCM problem A - Optimal Design of Mooring System
%
% Zhou Lvwen: zhou.lv.wen@gmail.com

% -------------------------------------------------------------------------
% Question 1
%
Lc = 22.05; chain = 2; M = 1200; depth = 18;
[tilt,elev,xsbed,xbuoy,f] = moor(Lc, chain, 12, 0, M, depth, 1);
[tilt,elev,xsbed,xbuoy,f] = moor(Lc, chain, 24, 0, M, depth, 1);

% -------------------------------------------------------------------------
% Question 2
%
Mi = 1200:10:4000; tilti = []; elevi = [];
for mi = Mi
    [tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,36,0,mi, depth);
    tilti = [tilti, tilt]; elevi = [elevi, elev];
end

figure('name','problem 2'); 
plot(Mi, tilti,'r', Mi, elevi, 'b');

xi = [min(Mi), max(Mi)];
hold on; plot(xi,[5,5],'r--',xi,[16,16],'b--');  
set(gca,'ytick',[0,5,10,16,20]); grid on
legend('Tilt angle of the drum','Elevation angle of the chain at the anchor')
xlabel('Mass of the heavy ball (kg)'); ylabel('Angle (degree)')

idx = find( tilti<=5 & elevi<=16 );
[tilt,elev,xbuoy,f] = moor(Lc, chain, 36, 0, min(Mi(idx)), depth, 1);

% -------------------------------------------------------------------------
% Question 3
Lc = 15.84; chain = 5; vw = 36; vs = 1.5; M = 4020; depth = 16;
[tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,vw,vs,M,depth,1);

Lc = 20.88; chain = 5; vw = 36; vs = 1.5; M = 4000; depth = 20;
[tilt,elev,xsbed,xbuoy,f] = moor(Lc,chain,vw,vs,M,depth,1);
