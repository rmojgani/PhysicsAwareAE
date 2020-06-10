function addpath_host(par)
% Edit 2017-04-24
%%
if nargin == 0; par = 0; end % par == 0 is Serial Code
if par == 0; fprintf('The code is initiated on Serial \n'); end
%% Machine dependent paths
[~,hostname]= system('hostname');
if ~isempty(strfind(hostname,'bridges'));
    disp(['ON SSH connection: ',hostname])
    % read pwd and split
    my_pwd=pwd;
    my_pwd=strsplit(my_pwd,'/');
    
    % my home directory
    my_home=strcat('/',my_pwd(2),'/',my_pwd(3));
    
    % ********************************************
    
    % add path to library
    path_to_library=strcat(my_home,'/library/matlab');
    addpath(path_to_library{1});  
elseif strncmpi(hostname(1,:),'rmojgani',7);
    disp(['ON Local: ',hostname])
    addpath('/media/rmojgani/Data/PhD Work/Paper/for rambod/matlab')
    addpath('/media/rmojgani/Data/PhD Work/Theory/Matlab Utilities/Computational Tools')
    addpath('/media/rmojgani/Data/PhD Work/Library/Mesh/Matlab Implementaion/Matmol_1.4/Tutorial_Examples/Dynamic_regridding/Burgers_movgrid')
elseif ~isempty(strfind(hostname,'golubh'));
    disp(['ON SSH connection: ',hostname])
    addpath('/home/mojgani2/library/matlab');       
end
%% Prarallel
% Note: use ALLPROFILES = parallel.clusterProfiles to know list of the
% profiles on your machine
if (par>0);
    fprintf('The code is initiated on Parallel \n');
    delete(gcp('nocreate'));
    if strfind(hostname,'rmojgani') == 1;           % Rambod on local machine
        par = min (par, 8);
        parpool('local8',par,'IdleTimeout', Inf)
    elseif or( ~isempty(strfind(hostname,'bridges')) ,...
               ~isempty(strfind(hostname,'golubh'))  );
        par = min (par, 72);
        parpool('local72',par,'IdleTimeout', Inf)   % Rambod on Bridges
    end
end