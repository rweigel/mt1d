% Compute and plot transfer function, phase, and impulse response
% for various profiles.

clear;

addpath([fileparts(mfilename('fullpath')),pathsep(),'m']);
addpath([fileparts(mfilename('fullpath')),pathsep(),'models']);

saveimg = true;

% Output files will be named ./figures/base_...
base = 'zplanewave_demo'; 
mu_0 = 4*pi*1e-7; % Vacuum permeability

ne = 5;
no = 1; 
N = no*10^ne; 
f = [1:N/2]/N;

Info = modelinfo();
titlestr = '';

for profile = 1:7
%for profile = 1:1

profilestr = sprintf('S%d',profile);
s = 1./Info.(profilestr).('rho');
h = Info.(profilestr).('thickness');

C    = zplanewave(s,h,f);   % = (1/(2*pi*f*j))*(Ex(f)/By(f)) = Ex(f)/B'y(f)
Cmag = sqrt(C.*conj(C));
Z    = 1j*2*pi*f.*C;         % Ex(f)/By(f)
Zmag = sqrt(Z.*conj(Z));

% Full array - zplanewave only returns f > 0.
Zf  = [0,Z,fliplr(conj(Z))]; 

% Full array for Ex(f)/B'y(f), B'y = dBy/dt
dZf = [0,C,fliplr(conj(C))];

% Apparent resistivity.  For infinite half-space will equal actual.
rho_a = C.*conj(C)*mu_0*2*pi.*f;
phi_C = (180/pi)*atan2(imag(C),real(C));

% Ex ~ exp(j*2*pi*f*t + phi) if By ~ exp(j*2*pi*f*t).
phi_Z = (180/pi)*atan2(imag(Z),real(Z));

tft  = [-N/2:1:N/2];
dhft = fftshift(ifft(dZf));
hft  = fftshift(ifft(Zf));

p = sqrt(rho_a./(2*pi*mu_0*f));
m = (4/pi)*(phi_Z*pi/180) - 1;
rho_nb = rho_a.*(1-m)./(1+m);

[~,idx] = min(abs(p - 100e3));
%p(idx)
if (isinf(h(1)))
    s = [s;s;s];
    h = 1e3*[1;100];
end
if length(h) > 1 && isinf(h(end))
    h(end) = 100*h(end-1);
end

figurex(1);clf;
    d = cumsum(h); % Depth
    %d = [10^(round(log10(d(1)))-1),d,10^(round(log10(d(end)))+1)];
    d = [10^2;d;10^6];
    s(end+1) = s(end);
    loglog(rho_nb,p/1e3,'.');
    hold on;
    for i = 1:length(s)-1
        loglog([1/s(i),1/s(i)],[d(i),d(i+1)]/1e3,'k','LineWidth',5);
        loglog([1/s(i),1/s(i+1)],[d(i+1),d(i+1)]/1e3,'k','LineWidth',1);
    end
    grid on;

    % Add some space to left and right.
    if length(unique(s)) == 1
        xl = get(gca,'XLim');
        set(gca,'XLim',[s(1)/10,s(1)*10]);
    else
        xl = get(gca,'XLim');
        set(gca,'XLim',[xl(1)*0.9,xl(2)*1.1]);
    end
    
    set(gca,'YDir','reverse');
    % Fix for https://www.mathworks.com/matlabcentral/answers/431242-reversing-y-axis-direction-causes-xticks-to-invert-and-overlap-with-labels
    xtl = get(gca,'XTickLabel');
    for i = 1:length(xtl)
        xtl{i} = [xtl{i},'_{ }'];
    end
    set(gca,'XTickLabel',xtl);    
    set(gca,'XAxisLocation','Top');

    xlabel('$\mbox{Resistivity}\;[\Omega\cdot \mbox{m}]$','Interpreter','Latex');
    lh = legend('$\rho_{NB}=\rho_a\,\frac{(\pi/2-\phi_{\tilde{C}})}{(\pi/2+\phi_{\tilde{C}})}\quad$',...
        '$\rho$','Location','SouthEast','Orientation','Horizontal');
    set(lh,'Interpreter','Latex');
    ylabel('\mbox{Depth}\,\mbox{[km]}','Interpreter','Latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');
    if saveimg
        figurep(base, profile, 'geometry')
    end

figurex(2);clf;
    loglog(f,rho_a,'b','LineWidth',3,'Marker','.','MarkerSize',10);
    hold on;
    loglog(f,rho_nb,'m','LineWidth',3,'Marker','.','MarkerSize',10);    
    loglog(f,Zmag*1e-3,'k','LineWidth',3,'Marker','.','MarkerSize',10);
    loglog(f,Cmag/1e3,'g','LineWidth',3,'Marker','.','MarkerSize',10);
    %loglog(f,p,'r','LineWidth',3,'Marker','.','MarkerSize',10);
    th = title(titlestr); 
    set(th,'Interpreter','Latex');    
    grid on;
    xlabel('$f \mbox{ [Hz]}$','Interpreter','Latex');
    set(gca,'XLim',[1/N-eps,0.51]);
    
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    lh = legend(' $\rho_a = \omega\mu_0|\widetilde{C}|^2\;[\Omega\cdot\mbox{m}]$',...
                ' $\rho_{NB}\; [\Omega\cdot\mbox{m}]$',...
                ' $\widetilde{Z}=|\widetilde{E}_x/\widetilde{B}_y| = \omega|\widetilde{C}|\;\left[\frac{\mbox{mV/km}}{\mbox{nT}}\right]$',...
                ' $|\widetilde{C}|=p=|\widetilde{E}_x/\widetilde{B}''_y|\;\mbox{[km]}$',...
                'Location','SouthEast','FontSize',14);

    set(lh,'Interpreter','Latex');
    set(gca,'FontName','Times');
        
    lh.BoxFace.ColorType = 'truecoloralpha';
    lh.BoxFace.ColorData = uint8(255*[1 1 1 0.5]');

    if saveimg
        figurep(base, profile, 'transferfn')
    end

figurex(3);clf;
    semilogx(f,phi_C,'g','LineWidth',3,'Marker','.','MarkerSize',10);
    hold on;
    semilogx(f,phi_Z,'k','LineWidth',3,'Marker','.','MarkerSize',10);
    grid on;
    set(gca,'YLim',[-90 90]);
    set(gca,'YTick',[-90:15:90]);
    th = title(titlestr);
    set(th,'Interpreter','Latex');  
    ylabel('$\mbox{[degrees]}$','Interpreter','Latex');
    xlabel('$f \mbox{ [Hz]}$','Interpreter','Latex');
    set(gca,'XLim',[1/N-eps,0.51]);    
    lh = legend(' $\phi_{\widetilde{C}}\quad$',...
                ' $\phi_{\widetilde{Z}}$',...
                'Location','SouthWest',...
                'Orientation','Horizontal');
    set(lh,'Interpreter','Latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');

    if saveimg
        figurep(base, profile, 'phase')
    end
 
figurex(4);clf;
    plot(tft,hft/1000,'k','LineWidth',2);
    hold on;grid on;
    plot(tft,dhft/1000,'g','LineWidth',2);
    %plot(dhft(N/2+2)./sqrt([1:60]),'b','LineWidth',2);
    set(gca,'XLim',[-30,60])
    th = title(titlestr);
    set(th,'Interpreter','Latex');      
    xlabel('$t\,\mbox{[s]}$','Interpreter','Latex');
    ylabel('$E_x(t)\,\mbox{[mV/km]}$','Interpreter','Latex');
    lh = legend('$B_y(t) = \delta(t)\,\mbox{[nT]}$',...
                '$B''_y(t) = \delta(t)\,\mbox{[nT/s]}$');
    set(lh,'Interpreter','Latex');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');
    if saveimg
        figurep(base, profile, 'irf')
    end
end % profile
