% simul_Milstein.m
%   Ver. 0.1 (15-Feb-2023)
%      Comments were updated for GitHub.
%   Ver. 0 (30-May-2022)
%%%%%
% This was copied from simul_MagMilstein.m (Ver. 0), and 
% was revised for the Milstein method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CAUTION: tranij means I(j,i)/h, not I(i,j)/h 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% output %%%
% yVec: a numerical solution for an SDE.
%%%%%%%%%%%%%%
%
TRAJ=150000;
wDimMax=5;
MaxK=128;
yDimMax=5;
%%%
randI=zeros(1,wDimMax*TRAJ);
rani=zeros(1,wDimMax*TRAJ);
tranij=zeros(1,wDimMax*(wDimMax-1)*TRAJ);
ran_diag=zeros(1,wDimMax*TRAJ);
%
yVec=zeros(yDimMax,TRAJ);
%
seed=5489;
%
%%%%% For ran_gene_using_genrand_real1_by_MATLAB
pi2_6=1.6449340668482264; % pi*pi/6 */
pi2x24=236.87050562614461; % 24*pi*pi */
sq2=1.4142135623730950; % sqrt(2)  */
sq2pi=4.4428829381583662; % sqrt(2)*pi */
%
tilX=zeros(1,wDimMax);
tilZ=zeros(wDimMax,wDimMax-1);
X_dv_k=zeros(wDimMax,MaxK);
Y=zeros(wDimMax,MaxK);
G=zeros(1,(wDimMax-1)*wDimMax/2);
xi_dv_sqstep=zeros(1,wDimMax);
loc_rani=zeros(1,wDimMax);
%%%%%
%
ydim = 2; % dimension of y
wdim = 2; % dimension of Wiener process
wdimM1=wdim-1;
if TRAJ<traject
    fprintf("traject is too large!");
    return;
end
Tend = 1.0;
transrijDim=wdim*wdimM1;
%
matF0=zeros(2,2);
matF1=zeros(2,2);
matF2=zeros(2,2);
matF0(1,1)=lam;
matF0(2,2)=lam;
matF1(1,1)=sig1;
matF1(2,2)=-sig1;
matF2(1,2)=sig2;
matF2(2,1)=sig2;
%
base_step=1.0/mm;  
errC=1.0/2.0; % 1.0/2.0
%
% Initialization */
ibase_index=0;
ibase=ibase_index;
for i=ibase:traject-1
    for j=1:wdim % starts from 1
        randI(ibase+j)=0;
    end
    ibase=ibase+wdim;
end
%
rng(seed,'twister'); % Setting a seed. */
tpoint=0;
yinit1=1; % 1st component of initial vector
yinit2=1; % 2nd component of initial vector
step = i_step*base_step;
eps=0.1*step;
sqstep=sqrt(step);
rann_empty=1;
%%%
loc_rann_empty=rann_empty;
%
for i=1:traject
    yVec(1,i)=yinit1;
    yVec(2,i)=yinit2;
end
%
while (eps<abs(Tend-tpoint))
    %
    ran_gene_using_genrand_real1_by_MATLAB;
    %%%
    ibase_index=0;
    ibase=ibase_index;
    for i=1:traject % starts from 1
        for j=1:wdim % starts from 1
            randI(ibase+j)=randI(ibase+j)+rani(ibase+j);
        end
         ibase=ibase+wdim;
    end
    %%%
    ibase_index=0;
    ibase=ibase_index;
    for i=1:traject % starts from 1
        % Calculation part for SDEs (begin)
        tmpVec=[yVec(1,i) yVec(2,i)]'+...
            matF0*[yVec(1,i) yVec(2,i)]'*step;
        tmpVec=tmpVec+...
            matF1*[yVec(1,i) yVec(2,i)]'*rani(ibase+1)*sqstep+...
            matF2*[yVec(1,i) yVec(2,i)]'*rani(ibase+2)*sqstep;
        j=1;
        tmp_ibase_index=ibase_index+(j-1)*wdimM1; % note j-1
        l=1; % note l=j:wdimM1
        tmpVec=tmpVec+...
            matF1*matF2*[yVec(1,i) yVec(2,i)]'*...
            tranij(tmp_ibase_index+l)*step;
        %
        j=2;
        tmp_ibase_index=ibase_index+(j-1)*wdimM1; % note j-1
        l=1; % note l=1:j-1
        tmpVec=tmpVec+...
            matF2*matF1*[yVec(1,i) yVec(2,i)]'*...
            tranij(tmp_ibase_index+l)*step;
        %%%
        j=1;
        tmpVec=tmpVec+...
            matF1*matF1*[yVec(1,i) yVec(2,i)]'*...
            ran_diag(ibase+j)*step;
        %
        j=2;
        tmpVec=tmpVec+...
            matF2*matF2*[yVec(1,i) yVec(2,i)]'*...
            ran_diag(ibase+j)*step;
        %%%%%
        yVec(1,i)=tmpVec(1);
        yVec(2,i)=tmpVec(2);
        % Calculation part for SDEs (end)
        ibase=ibase+wdim;
        ibase_index=ibase_index+transrijDim;
    end
    tpoint= tpoint+step;
end % end of while related to eps





