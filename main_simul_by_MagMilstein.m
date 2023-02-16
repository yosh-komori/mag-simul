% main_simul_by_MagMilstein.m
%   Ver. 0.1 (15-Feb-2023)
%      Comments were updated for GitHub.
%   Ver. 0 (30-May-2022)
%      This program manages simulations, and it saves
%      the results in files.
%%%
%%% input
mm=512; % this is for the base step size: 2^{-mm}.
i_step=256; % this is for a step size: i_step*(base step size).
traject=1000; % number of trajectories.
lam=-1.0/4; % a parameter in an SDE.
sig1=1.0/2; % a parameter in an SDE.
sig2=2.0/5; % a parameter in an SDE.
%%% output %%%
% yVec: a numerical solution for an SDE.
%%%%%%%%%%%%%%
%
simul_MagMilstein;
%
formatSpec=...
    'yVec_MagMilstein_step_%ddiv%d_tr%d_lm_%4.3f_s1_%3.2f_s2_%3.2f_t_%2.1f';
tmpData=[yVec(1,1:traject)' yVec(2,1:traject)'];
tmpStr=sprintf(formatSpec,i_step,mm,traject,lam,sig1,sig2,Tend);
save(tmpStr, 'tmpData','-ASCII');
