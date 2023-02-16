% main_simul_by_MagMilstein.m
%   Ver. 1.1 (16-Feb-2023)
%      Comments were updated for GitHub.
%   Ver. 1 (30-Aug-2022)
%      1) Revised to set seed=5489 and rng(seed,'twister') in this file.
%      2) Revised to output a file for each batch.
%%%
%   Ver. 0 (30-May-2022)
%      This program manages simulations, and it saves
%      the results in files.
%%%
seed=5489;
rng(seed,'twister'); % Setting a seed. */
%%% input
mm=2; % this is for the base step size: 2^{-mm}.
i_step=1; % this is for a step size: i_step*(base step size).
traject=100000; % number of trajectories.
batchMax=10; % number of batches.
lam=-0.2; % a parameter in an SDE.
sig1=1.0; % a parameter in an SDE.
sig2=1.0; % a parameter in an SDE.
%%% output %%%
% yVec: a numerical solution for an SDE.
%%%%%%%%%%%%%%
%
for ib=1:batchMax
    simul_MagMilstein;
    %
    formatSpec=...
        'yVec_MagMilstein_step_%ddiv%d_tr%d_b%d_lm_%4.3f_s1_%3.2f_s2_%3.2f_t_%2.1f';
    tmpData=[yVec(1,1:traject)' yVec(2,1:traject)'];
    tmpStr=sprintf(formatSpec,i_step,mm,traject,ib,lam,sig1,sig2,Tend);
    save(tmpStr, 'tmpData','-ASCII');
end
