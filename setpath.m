%% Set path

addpath(genpath(pwd),'-begin');
rmpath(genpath('./.git'));
rmpath(genpath('./docs'));
rmpath(genpath('./generate'));
savepath;
clear;clc;