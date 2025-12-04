function install()
    baseDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(baseDir, 'src'));
    addpath(genpath(fullfile(baseDir, 'src')));
    addpath(genpath(fullfile(baseDir, 'data')));
    addpath(genpath(fullfile(baseDir, 'tests')));
    addpath(genpath(fullfile(baseDir, 'tutorials')));
    addpath(fullfile(baseDir, 'bin')); 
    savepath; clc;
    fprintf('MagPackage installed successfully.\n');
end