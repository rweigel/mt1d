% Executes three tests for ZPLANEWAVE
%
% First two tests compare direct calculation with zplanewave calculation.
% Third test creates plots for comparison with Figure 2.4(c) of
% Simpson and Bahr.

clear;

saveimg = true;

% Output files will be named ./figures/base_...
base = 'zplanewave_test'; 

logfile = [mfilename('fullpath'),'.log'];
delete(logfile);
diary(logfile);

mu_0 = 4*pi*1e-7; % Vacuum permeability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infinite half-space
f = 0.1;
h = Inf;
s = 2;
C = zplanewave(s,h,f);

Ctest = sqrt(2)/( (1j+1)*sqrt(mu_0*s*2*pi*f) );

% Should get conductivty back.
stest = 1/(C*conj(C)*mu_0*2*pi*f);

fprintf('Test 1: Infinite half-space, sigma = %d, f = 0.1.\n',s);

rho_a = C.*conj(C)*mu_0*2*pi*f;
phi_a = (180/pi)*atan2(imag(C),real(C));
fprintf('Test 1     error/eps : %6.2f\n',(s-stest)/eps);
fprintf('Test 1         rho_a : %6.2f [Ohm m]\n',rho_a);
fprintf('Test 1         phi_a : %6.2f [degrees]\n',phi_a);
fprintf('Test 1 Re(error)/eps : %6.2f\n',real(C-Ctest(1))/eps);
fprintf('Test 1 Im(error)/eps : %6.2f\n',imag(C-Ctest(1))/eps);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two layers
f = 0.1;
s = [100,1]';
p = sqrt(2./(mu_0*s(1)*2*pi*f));
h = [0.1*p(1),Inf];

C = zplanewave(s,h,f);

% Manual calculation of C.
q(2) = sqrt(1j*mu_0*s(2)*2*pi*f);
Ctest(2) = 1/q(2); 

q(1) = sqrt(1j*mu_0*s(1)*2*pi*f);
Ctest(1) = (1/q(1)) * ( q(1)*Ctest(2) + tanh(q(1)*h(1)) ) / (1 + q(1)*Ctest(2)*tanh(q(1)*h(1)));

fprintf('Test 2: Two layers\n');
fprintf('Test 2: Two layers,\n');
fprintf('        sigma = [100,1],');
fprintf('        h     = (penetration depth of first layer)/10,\n');
fprintf('        f     = 0.1.\n');
rho_a = C.*conj(C)*mu_0*2*pi*f;
sigma_a = 1/rho_a;

phi_a = (180/pi)*atan2(imag(C),real(C));

fprintf('Test 2         rho_a   : %6.2f [Ohm m]\n',rho_a);
fprintf('Test 2         sigma_a : %6.2f [1/(Ohm m)]\n',sigma_a);
fprintf('Test 2         phi_a   : %6.2f [degrees]\n',phi_a);
fprintf('Test 2 Re(error)/eps   : %6.2f\n',real(C-Ctest(1))/eps);
fprintf('Test 2 Im(error)/eps   : %6.2f\n',imag(C-Ctest(1))/eps);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison with Figure 2.4(c) of Simpson and Baer
fprintf('Test 3: Visual comparison of Figure 2.4(c) of Simpson and Baer\n');

T = logspace(-3,3,30);
f = 1./T;

s = [1/1000,1/10,1/100,1/5]';
h = 1e3*[10,20,400,Inf];

C    = zplanewave(s,h,f);
Z    = 1j*2*pi*f.*C;
Zmag = sqrt(Z.*conj(Z));

rho_a = C.*conj(C)*mu_0*2*pi.*f;
phi_C = (180/pi)*atan2(imag(C),real(C));
phi_Z = (180/pi)*atan2(imag(Z),real(Z));

figure(1);clf;
    d = cumsum(h);
    d = [10^3,d(1:end-1),10^6];
    s(end+1) = s(end);
    for i = 1:length(s)-1
        loglog([1/s(i),1/s(i)],[d(i),d(i+1)]/1e3,'k','LineWidth',5);
        hold on;
        loglog([1/s(i),1/s(i+1)],[d(i+1),d(i+1)]/1e3,'k','LineWidth',1);
    end
    grid on;
    xlabel('Resistivity [\Omega\cdotm]');
    ylabel('Depth [km]');
    axis([.9 1.1*10^3 1 10^3]);

    set(gca,'YDir','reverse');
    % Fix for https://www.mathworks.com/matlabcentral/answers/431242-reversing-y-axis-direction-causes-xticks-to-invert-and-overlap-with-labels
    set(gca,'XAxisLocation','Top');
    xtl = get(gca,'XTickLabel');
    for i = 1:length(xtl)
        xtl{i} = [xtl{i},'_{ }'];
    end
    set(gca,'XTickLabel',xtl);    

    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');
    if saveimg
        figurep(base, 1, 'phi_geometry')
    end

figure(2);clf;
    loglog(1./f,Zmag*1e-3,'k','LineWidth',3,'Marker','.','MarkerSize',20);
    hold on;
    loglog(1./f,rho_a,'b','LineWidth',3,'Marker','.','MarkerSize',20);
    grid on;
    xlabel('Period [s]');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    lh = legend(' $|\widetilde{Z}| = |\widetilde{E}_x/\widetilde{B}_y| = \omega|\widetilde{C}|\;\left[\frac{\mbox{mV/km}}{\mbox{nT}}\right]$',...
                ' $\rho_a = \omega\mu_0|\widetilde{C}|^2\quad[\Omega\cdot\mbox{m}]$');
    set(lh,'Interpreter','Latex','Location','SouthWest','FontSize',14);
    set(gca,'FontName','Times');
    if saveimg
        figurep(base, 1, 'rho_a')
    end
    
figure(3);clf;
    semilogx(1./f,phi_Z,'k','LineWidth',3,'Marker','.','MarkerSize',20);
    hold on;
    semilogx(1./f,phi_C,'b','LineWidth',3,'Marker','.','MarkerSize',20);
    grid on;
    set(gca,'YLim',[-90 90]);
    set(gca,'YTick',[-90:15:90]);
    xlabel('Period [s]');
    lh = legend(' $\phi_{\widetilde{Z}}\quad\mbox{[deg]}$',...
                ' $\phi_{\widetilde{C}}\quad\mbox{[deg]}$',...
                'Location','NorthWest');
    set(lh,'Interpreter','Latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');
    if saveimg
        figurep(base, 1, 'phi_a')
    end

diary off