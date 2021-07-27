clear;

saveimg = true;

addpath([fileparts(mfilename('fullpath')),pathsep(),'m']);
addpath([fileparts(mfilename('fullpath')),pathsep(),'models']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models to compare
models     = {'Q1','CO1','BR1'};
fnamebase  = models{1};
for i = 2:length(models)
    fnamebase = [fnamebase,'_',models{i}];
end
fname{1} = sprintf('z_planewave_model_compare_depth_vs_rho_%s',fnamebase);
fname{2} = sprintf('z_planewave_model_compare_Z_vs_T_%s',fnamebase);
fname{3} = sprintf('z_planewave_model_compare_phi_vs_T_%s',fnamebase);
fname{4} = sprintf('z_planewave_model_compare_irf_%s',fnamebase);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ne = 5;
no = 1; 
N = no*10^ne; 
f = [1:N/2]/N;
T = 1./f;
mu_0 = 4*pi*1e-7; % Vacuum permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Info = modelinfo();

tmp1 = ' $|\widetilde{Z}| = |\widetilde{E_x}/\widetilde{B_y}|';
tmp2 = '= \omega\cdot|\widetilde{C}|\quad\mbox{[V/m/T]}$';
legstr{1} = [tmp1,tmp2];

for i = 1:length(fname)
    figurex(i);clf;
end

for i = 1:length(models)
    thickness = Info.(models{i}).('thickness');
    d = cumsum(thickness(1:end-1));
    if length(thickness) > 1
        dmin(i) = 10^(round(log10(d(1)))-1);
        dmax(i) = 10^(round(log10(d(end)))+1);
    else
        dmin(i) = nan;
        dmax(i) = nan;
    end
end
dmin = nanmin(dmin);
dmax = nanmin(dmax);

for i = 1:length(models)
    legendstrings{i} = Info.(models{i}).('longname');
    rho = Info.(models{i}).('rho');
    thickness = Info.(models{i}).('thickness');
    sigma = 1./rho;
    
    C{i}    = zplanewave(1./rho',thickness',f);
    Z{i}    = 1j*2*pi*f.*C{i};
    Zmag{i} = sqrt(Z{i}.*conj(Z{i}));

    rho_a{i} = C{i}.*conj(C{i})*mu_0*2*pi.*f;
    phi_C{i} = (180/pi)*atan2(imag(C{i}),real(C{i}));
    phi_Z{i} = (180/pi)*atan2(imag(Z{i}),real(Z{i}));

    % Full array - zplanewave only returns f > 0.
    Zf   = [0,Z{i},fliplr(conj(Z{i}))]; 
    % Full array for Ex(f)/B'y(f), B'y = dBy/dt
    dZf  = [0,C{i},fliplr(conj(C{i}))];
    dhft = fftshift(ifft(dZf));

    hft{i} = fftshift(ifft(Zf));
    tft{i} = [-N/2:1:N/2];

    figurex(1);
        c = get(gca,'ColorOrder');
        %axis([.7 1.1*10^4 1 1.1*10^4]);

        if (i == 1)
            for z = 1:3
                loglog([1,1],[NaN,NaN],'Color',c(z,:),'LineWidth',2);
                hold on;grid on;
            end
            ax = gca;ax.ColorOrderIndex = 1;
        end

        if length(sigma) > 1
            d = [dmin;d;dmax];
            sigma(end+1) = sigma(end);
            for l = 1:length(sigma)-1
                loglog([1/sigma(l),1/sigma(l)],[d(l),d(l+1)]/1e3,...
                       'Color',c(i,:),'LineWidth',2);
            if (i==1),hold on;grid on;end		
                loglog([1/sigma(l),1/sigma(l+1)],[d(l+1),d(l+1)]/1e3,...
                       'Color',c(i,:),'LineWidth',1);
            end
            set(gca,'YDir','reverse');
            xlabel('Resistivity [\Omega\cdot m]');
            ylabel('Depth [km]');
        else
            loglog([1/sigma(1),1/sigma(1)],[dmin,dmax]/1e3,...
                    'Color',c(i,:),'LineWidth',2);
        end
        
        % Fix for https://www.mathworks.com/matlabcentral/answers/431242-reversing-y-axis-direction-causes-xticks-to-invert-and-overlap-with-labels
        xtl = get(gca,'XTickLabel');
        for l = 1:length(xtl)
            xtl{l} = [xtl{l},'_{ }'];
        end
        set(gca,'XTickLabel',xtl);    
        set(gca,'XAxisLocation','Top');

    figurex(2);
        loglog(1./f,Zmag{i},'-','LineWidth',2);
        hold on;grid on;
        xlabel('Period [s]');
        if 0
            yl = get(gca,'YLim');
            xl = get(gca,'XLim');
            text(xl(1),yl(2),...
             ' $|\widetilde{Z}| \;\left[\frac{\mbox{mV/km}}{\mbox{nT}}\right]$','Interpreter','Latex',...
             'VerticalAlignment','bottom');
        end
        ylabel('$|\widetilde{Z}|$ [V/m/T]','Interpreter','Latex');
        set(gca,'XLim',[min(1./f),max(1./f)])

    figurex(3);
        semilogx(1./f,phi_Z{i},'LineWidth',2);
        hold on;grid on;
        ylabel(' $\phi_{|\widetilde{Z}|}$ [deg]','Interpreter','Latex');
        set(gca,'YLim',[-90,90]);
        set(gca,'YTick',[-90:15:90]);
        xlabel('Period [s]');
        set(gca,'XLim',[min(1./f),max(1./f)])	
		
    figurex(4);
        plot(tft{i},hft{i}/1000,'LineWidth',2);
        hold on;grid on;
        set(gca,'XLim',[-20,30])
        xlabel('$t\,\mbox{[s]}$','Interpreter','Latex');
        ylabel('$E_x(t)\,\mbox{[mV/m]}$ for $B_y(t) = \delta(t)\,\mbox{[nT]}$','Interpreter','Latex');
end

for i = 1:6
    figurex(i);
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    set(gca,'FontName','Times');    
end

figurex(1);
    legend(legendstrings,'Location','NorthWest');
figurex(2);
    legend(legendstrings,'Location','NorthEast');
figurex(3);
    legend(legendstrings,'Location','SouthWest');
figurex(4);
    legend(legendstrings,'Location','SouthWest');

if saveimg
    for i = 1:length(fname)
        figurex(i);
        figurep(fname{i});
    end
end
