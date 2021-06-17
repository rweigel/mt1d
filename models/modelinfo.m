function M = modelinfo()

filename = mfilename('fullpath');
path = fileparts(filename); % Directory of this script.

% Read file from disk if found.
matfile = [path,filesep(),'modelinfo.mat'];
if exist(matfile,'file')
    load(matfile,'M');
    return
end

M = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quebec Models
modelsQ     = {'Q1','Q2','Q3'};
modelnamesQ = {'Quebec Model 1','Quebec Model 2','Quebec Model 3'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USGS Models
modelsU = {...
            'AK1','AK2','AP1','AP2','BR1','CL1','CO1',...
            'CP1','CP2','CS1','IP1','IP2',...
            'IP3','IP4','NE1','PB1','PB2','PT1','SL1','SU1'...
          };

modelnamesU = {...
    'Adirondack Mountains',...
    'Appalachan Plateaus',...
    'Northern Appalachan Plateaus',...
    'Northwest Basin and Range',...
    'Colorado Plateau',...
    'Columbia Plateau',...
    'Coastal Plain (South Carolina)',...
    'Coastal Plain (Georgia)',...
    'Cascade-Sierra Mountains',...
    'Florida Peninsula',...
    'Interior Plains (North Dakota)',...
    'Interior Plains',...
    'Interior Plains (Michigan)',...
    'Interior Plains (Great Plains)',...
    'New England',...
    'Pacific Border (Willamette Valley)',...
    'Pacific Border (Puget Lowlands)',...
    'Piedmont',...
    'St. Lawrence Lowlands',...
    'Superior Upland'...
};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Models
modelsS     = {'S1','S2','S3','S4','S5','S6','S7'};
modelnamesS = {...
    'Constant',...
    'Two Layer - High/Low Resistivity',...
    'Two Layer - Low/High Resistivity',...
    'Decreasing Resistivity',...
    'Increasing Resistivity',...
    'High Resistivity Sheet',...
    'Low Resistivity Sheet'...
   };
rhoS = {};
tS = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(modelnamesS)
    rho_t = load(sprintf('%s/Simple/%s.txt',path,modelsS{i}));

    M.(modelsS{i}) = struct();
    M.(modelsS{i}).('shortname') = modelsS{i};
    M.(modelsS{i}).('longname') = modelnamesS{i};
    M.(modelsS{i}).('rho') = rho_t(:,1);
    M.(modelsS{i}).('rho_units') = 'ohm-m';
    M.(modelsS{i}).('thickness') = rho_t(:,2);
    M.(modelsS{i}).('thickness_units') = 'm';
end

for i = 1:length(modelnamesQ)
    rho_t = load(sprintf('%s/Quebec/%s.txt',path,modelsQ{i}));

    M.(modelsQ{i}) = struct();
    M.(modelsQ{i}).('shortname') = modelsQ{i};
    M.(modelsQ{i}).('longname') = modelnamesQ{i};
    M.(modelsQ{i}).('rho') = rho_t(:,1);
    M.(modelsQ{i}).('rho_units') = 'ohm-m';
    M.(modelsQ{i}).('thickness') = rho_t(:,2);
    M.(modelsQ{i}).('thickness_units') = 'm';
end

for i = 1:length(modelnamesU)
    M.(modelsU{i}) = struct();
    M.(modelsU{i}).('shortname') = modelsU{i};
    M.(modelsU{i}).('longname') = modelnamesU{i};
    rho_t = load(sprintf('%s/USGS/%s_GroundModel.txt',path,modelsU{i}));
    M.(modelsU{i}).('rho') = rho_t(:,1);
    M.(modelsU{i}).('rho_units') = 'ohm-m';
    M.(modelsU{i}).('thickness') = rho_t(:,2);
    M.(modelsU{i}).('thickness_units') = 'm';
end

save(matfile,'M');
