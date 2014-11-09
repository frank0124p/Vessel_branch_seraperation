%Vessel_Separation
clc;close all;clear
%Add big branch radius and pair to the node
path=pwd;
path=strcat(pwd,'\Data');
addpath(path);

load center_pts.mat % parameter-result 
load tubular_response_bigBranch.mat %parameter --smooth_result with big branch
load BranchNode.mat %key point detection

big_branch=smooth_result;
