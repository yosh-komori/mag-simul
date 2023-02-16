% output_example_for_stab.m
%   Ver. 0 (16-Feb-2023)
%%%%%
% If you carry out 
%    main_simul_by_Milstein.m,
% then you will have 10 files such as 
%"yVec_MagMilstein_step_1div2_tr100000_b1_lm_-0.200_s1_1.00_s2_1.00_t_5.0".
%
% After that, this program shows mean and standard deviation
% related to ||y vector||^2 for h=1/2;
%%%%%
%
str1='yVec_MagMilstein_step_1div2_tr100000_';
str3='_lm_-0.200_s1_1.00_s2_1.00_t_5.0';
%
ibmax=10;
tmpVList=zeros(1,ibmax);
%
warning('off','MATLAB:namelengthmaxexceeded');
for ii=1:ibmax
    str2=append('b',num2str(ii));
    tmpfname=append(str1,str2,str3);
    tmpSol=load(tmpfname);
    %
    tmpVList(ii) = sum(tmpSol(:,1).^2+tmpSol(:,2).^2)/length(tmpSol(:,1));
end
warning('on','MATLAB:namelengthmaxexceeded');
%
meanV=mean(tmpVList);
stdV=sqrt(sum((tmpVList-meanV).^2)/length(tmpVList));
fprintf("(Mean,Std) for h=1/2: (%f,%f)\n",meanV,stdV);
