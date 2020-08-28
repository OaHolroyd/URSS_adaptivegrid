%% RUN CUSTOM
% does a single run with a custom setup.

test_setup;

dx = 0.25;
dt = dx;
dxsup = 0.2;
dxinf = 0.1;
K = 1.1;
d = 0;

Tout = [0,Tmax]';
interpmethod = 'pchip';

cl = clock;
filename = ['test_',name,'_uniform'];
save('bparams_custom','g','f','L','Tmax','b','R','C');
save('mparams_custom','dx','dxsup','dxinf','K','d','dt','Tout','interpmethod','filename');
disp(['starting ',filename, ' (at ',num2str(cl(4),'%02.0f'),':',num2str(cl(5),'%02.0f'),')'])
record_all;
disp(['    finished (',num2str(elapsed_time),'secs)'])

delete bparams_custom.mat
delete mparams_custom.mat
