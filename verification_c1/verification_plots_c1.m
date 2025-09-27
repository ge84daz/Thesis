clear all

%--------------------- Load uzPz_xyf_hscyl_fem_midline_4node --------------
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
matPath = fullfile(thisDir, ['variables_verification/uzPz_xyf_hscyl_fem_midline_4node.mat']);
tmp = load(matPath);
reqPath = {'uzPz_xyf_hscyl_fem_midline_4node'};

for k = 1:numel(reqPath)
        tmp = tmp.(reqPath{k});
end
uzPz_xyf_hscyl_fem_midline_4node = tmp;


%--------------------- Load uzPz_xyf_hscyl_fem_midline_3node --------------
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
matPath = fullfile(thisDir, ['variables_verification/uzPz_xyf_hscyl_fem_midline_3node.mat']);
tmp = load(matPath);
reqPath = {'uzPz_xyf_hscyl_fem_midline_3node'};

for k = 1:numel(reqPath)
        tmp = tmp.(reqPath{k});
end
uzPz_xyf_hscyl_fem_midline_3node = tmp;


%--------------------- Load uzPz_xyf_hscyl_fem_midline_3node --------------
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
matPath = fullfile(thisDir, ['variables_verification/nod_midline.mat']);
tmp = load(matPath);
reqPath = {'nod_midline'};

for k = 1:numel(reqPath)
        tmp = tmp.(reqPath{k});
end
nod_midline = tmp;

%--------------------- Input Variables (might change)----------------------
f_tot = [-30,0,30];
fev = 30;
x = -64:2:62;
%plot settings
color_cyan = [0,0.447000000000000,0.741000000000000];
MarkerSize = 5;
LineWidth = 1.3;

%-------------------------- PLOT THAT SHIT --------------------------------

plot(nod_midline(:,2),abs(uzPz_xyf_hscyl_fem_midline_4node(size(x,2)/2+1,:,f_tot==-fev)),'color',color_cyan,'MarkerSize',MarkerSize,'LineWidth',LineWidth)



