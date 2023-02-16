% output_example.m
%   Ver. 0 (15-Feb-2023)
%%%%%
% If you have carried out 
%    main_simul_by_Milstein.m, main_simul_by_MagEuler.m,
%    and main_simul_by_MagMilstein,
% then this program shows the errors of Magnus-type Euler
% and Magnus-type Milstein methods for h=1/2 in a noncommutative
% test SDE.
%%%%%
%
RefSol=load('yVec_Milstein_step_1div512_tr1000_lm_-0.25_s1_0.50_s2_0.40');
MESol=load(...
    'yVec_MagEuler_step_256div512_tr1000_lm_-0.25_s1_0.50_s2_0.40_t_1.0');
MMSol=load(...
 'yVec_MagMilstein_step_256div512_tr1000_lm_-0.250_s1_0.50_s2_0.40_t_1.0');
%
tmp=MESol-RefSol;
err1 = sqrt(sum(tmp(:,1).^2+tmp(:,2).^2))/length(tmp(:,1));
fprintf("Error of MagEuler for h=1/2: %f\n",log2(err1));
%
tmp=MMSol-RefSol;
err2 = sqrt(sum(tmp(:,1).^2+tmp(:,2).^2))/length(tmp(:,1));
fprintf("Error of MagMilstein for h=1/2: %f\n",log2(err2));
