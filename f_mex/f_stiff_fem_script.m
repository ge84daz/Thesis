% F_STIFF_FEM_SCRIPT_N   Generate MEX-function f_stiff_fem_mex from f_stiff_fem.
% 
% Script generated from project 'f_stiff_fem.prj' on 06-Dec-2019.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.EnableMexProfiling = true;
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'f_stiff_fem'.
ARGS = cell(1,1);
ARGS{1} = cell(16,1);
ARGS{1}{1} = coder.typeof('X',[1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[Inf  13],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{5} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{6} = coder.typeof(0,[Inf  13],[1 0]);
ARGS{1}{7} = coder.typeof(0);
ARGS{1}{8} = coder.typeof(0);
ARGS{1}{9} = coder.typeof(1i,[6 6]);
ARGS{1}{10} = coder.typeof(0);
ARGS{1}{11} = coder.typeof(1i,[6 6]);
ARGS{1}{12} = coder.typeof(0);
ARGS{1}{13} = coder.typeof(0);
ARGS{1}{14} = coder.typeof(0);
ARGS{1}{15} = coder.typeof(0);
ARGS{1}{16} = coder.typeof(0);

%% Invoke MATLAB Coder.
%cd(path.functions);

% codegen options function -args {func_inputs}
% -config  = use configuration in cfg. ...
% -o       = output_file_name
% -args    = func_inputs in Cell arry ARGS{1}

codegen -config cfg f_stiff_fem -args ARGS{1}

cd(path.current)