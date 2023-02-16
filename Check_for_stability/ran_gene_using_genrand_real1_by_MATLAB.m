% ran_gene_using_genrand_real1_by_MATLAB
%   Ver. 0.1 (15-Feb-2023)
%      Comments were updated for GitHub.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For an approximation to a double stochastic
%   integral, this program uses the polar method
%   [Knuth:1981, p. 118] and the approximation
%   method by [Wiktorsson:2001].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Ver. 0 (17-May-2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CAUTION: tranij means I(j,i)/h, not I(i,j)/h 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
rate_step=step/base_step;
sq_rate_step=sqrt(rate_step);
eps=0.1*base_step;
%
% Preparation */
trajectM1=traject-1;
wdimM1=wdim-1;
%
if 0~=wdimM1
    tranijDim=wdim*wdimM1;
    gdim=tranijDim/2;
    gdimM1=gdim-1;
    maxDiscardNum=MaxK*wdim+gdim;
    CBaseStep24PI2=errC*base_step*pi2x24;
end
%
% Initializing rani, ran_diag and tmpStep */
ibase=0;
for ii=ibase:trajectM1
    for jj=1:wdim % starts from 1
        rani(ibase+jj)=0;
        ran_diag(ibase+jj)=0;
    end
    ibase=ibase+wdim;
end
tmpStep=0.0;
%
if 0~=wdimM1
    % Initializing tranij */
    ibase_index=0;
    for ii=ibase_index:trajectM1
      for ll=1:wdim % starts from 1
          tmp_ibase_index=ibase_index+(ll-1)*wdimM1; % note ll-1
          for jj=1:wdimM1 % starts from 1
              tranij(tmp_ibase_index+jj)=0;
          end
      end
      ibase_index=ibase_index+tranijDim;
    end
end
%
while (eps<abs(step-tmpStep))
    ibase_index=0;
    ibase=ibase_index;
    for ii=ibase:trajectM1
        %
        % Generating delW/sqrt(h) and Iii/h */
        for jj=1:2:wdimM1 % starts from 1
            while (1)
                ranu1 = 2*rand(1,1)-1;
                ranu2 = 2*rand(1,1)-1;
                caps = ranu1*ranu1 + ranu2*ranu2;
                if ((1>caps) && (0<caps))
                    break;
                end
            end
            caps = sqrt(-2.0*log(caps)/caps);
            if (1==loc_rann_empty)
                loc_rani(jj) = ranu1*caps;
                loc_rani(jj+1) = ranu2*caps;
            else
                loc_rani(jj) = loc_rann_rest;
                loc_rani(jj+1) = ranu1*caps;
                loc_rann_rest = ranu2*caps;
                loc_rann_empty=0;
            end
            rani(ibase+jj)=rani(ibase+jj) + loc_rani(jj);
            rani(ibase+jj+1)=rani(ibase+jj+1) + loc_rani(jj+1);
        end
        %
        if (wdim~=(jj+2)-1) % note (jj+2)-1
            if (0==loc_rann_empty)
                loc_rani(wdim) = loc_rann_rest; % note wdim
                loc_rann_empty = 1;
            else
                while (1)
                    ranu1 = 2*rand(1,1)-1;
                    ranu2 = 2*rand(1,1)-1;
                    caps = ranu1*ranu1 + ranu2*ranu2;
                    if ((1>caps) && (0<caps))
                        break;
                    end
                end
                caps = sqrt(-2.0*log(caps)/caps);
                loc_rani(wdim) = ranu1*caps; % note wdim
                loc_rann_rest = ranu2*caps;
                loc_rann_empty = 0;
            end % note wdim
            rani(ibase+wdim)=rani(ibase+wdim) + loc_rani(wdim);
        end
        %
        for jj=1:wdim % starts from 1
            ran_diag(ibase+jj)=ran_diag(ibase+jj)+...
                (randI(ibase+jj)*sq_rate_step+...
                rani(ibase+jj)-loc_rani(jj))*loc_rani(jj)+...
                (loc_rani(jj)*loc_rani(jj)-1)/2.0;
        end
        %
        if 0~=wdimM1
            %Maximum of k in a traject */
            rani_norm2=0;
            for jj=1:wdim % starts from 1
                rani_norm2=rani_norm2+loc_rani(jj)*loc_rani(jj);
            end
            %
            tmp=sqrt(tranijDim*(wdim+4*rani_norm2)/CBaseStep24PI2);
            locMaxK=0; % note while, not for
            while (locMaxK<tmp)
                locMaxK=locMaxK+1;
            end
            if MaxK<locMaxK
                fprintf("Error: locMaxK is too large ");
                fprintf("in ran_gene_using_genrand_real1!\n");
                return;
            end
            %
            locMaxKM1=locMaxK-1;
            %
            % Calculating sqrt(an) */
            locSqAn=0;
            for kk=locMaxK:-1:1 % note ":1"
                locSqAn=locSqAn+1.0/(kk*kk);
            end
            %
            locSqAn=sqrt(pi2_6-locSqAn);
            %
            for jj=1:wdim % starts from 1
                % Generating X/k and tilX for jj */
                for kk=1:2:locMaxKM1 % starts from 1
                    while (1)
                        ranu1 = 2*rand(1,1)-1;
                        ranu2 = 2*rand(1,1)-1;
                        caps = ranu1*ranu1 + ranu2*ranu2;
                        if ((1>caps) && (0<caps))
                            break;
                        end
                    end
                    caps = sqrt(-2.0*log(caps)/caps);
                    if(1==loc_rann_empty)
                        X_dv_k(jj,kk) = ranu1*caps;
                        X_dv_k(jj,kk+1) = ranu2*caps;
                    else
                        X_dv_k(jj,kk) = loc_rann_rest;
                        X_dv_k(jj,kk+1) = ranu1*caps;
                        loc_rann_rest = ranu2*caps;
                        loc_rann_empty=0;
                    end
                end
                %---%
                if 0==length(kk) % note this is for MATLAB
                    kk=1;
                end
                %---%
                if ((kk+2)-1)~=locMaxK % note (kk+2)-1
                    if(0==loc_rann_empty)
                        X_dv_k(jj,locMaxK) = loc_rann_rest; % note locMaxK
                        loc_rann_empty = 1;
                    else
                        while (1)
                            ranu1 = 2*rand(1,1)-1;
                            ranu2 = 2*rand(1,1)-1;
                            caps = ranu1*ranu1 + ranu2*ranu2;
                            if ((1>caps) && (0<caps))
                                break;
                            end
                        end
                        caps = sqrt(-2.0*log(caps)/caps);
                        X_dv_k(jj,locMaxK) = ranu1*caps; % note locMaxK
                        loc_rann_rest = ranu2*caps;
                        loc_rann_empty = 0;
                    end
                end
                %
                for kk=1:locMaxK % starts from 1
                    X_dv_k(jj,kk)=X_dv_k(jj,kk)/(kk); % note kk
                end
                tilX(jj)=0;
                for kk=locMaxK:-1:1 % note locMaxK and ":1"
                    tilX(jj)=tilX(jj)+X_dv_k(jj,kk);
                end
                % Completed generating X/k and tilX for jj */
                %%%
                %
                % Generating Y for jj */
                for kk=1:2:locMaxKM1 % starts from 1
                    while (1)
                        ranu1 = 2*rand(1,1)-1;
                        ranu2 = 2*rand(1,1)-1;
                        caps = ranu1*ranu1 + ranu2*ranu2;
                        if ((1>caps) && (0<caps))
                            break;
                        end
                    end
                    caps = sqrt(-2.0*log(caps)/caps);
                    if 1==loc_rann_empty
                        Y(jj,kk) = ranu1*caps;
                        Y(jj,kk+1) = ranu2*caps;
                    else
                        Y(jj,kk) = loc_rann_rest;
                        Y(jj,kk+1) = ranu1*caps;
                        loc_rann_rest = ranu2*caps;
                        loc_rann_empty=0;
                    end
                end
                %---%
                if 0==length(kk) % note this is for MATLAB
                    kk=1;
                end
                %---%
                if ((kk+2)-1)~=locMaxK % note (kk+2)-1
                    if (0==loc_rann_empty)
                        Y(jj,locMaxK) = loc_rann_rest; % note locMaxK
                        loc_rann_empty = 1;
                    else
                        while (1)
                            ranu1 = 2*rand(1,1)-1;
                            ranu2 = 2*rand(1,1)-1;
                            caps = ranu1*ranu1 + ranu2*ranu2;
                            if ((1>caps) && (0<caps))
                                break;
                            end
                        end
                        caps = sqrt(-2.0*log(caps)/caps);
                        Y(jj,locMaxK) = ranu1*caps; % note locMaxK
                        loc_rann_rest = ranu2*caps;
                        loc_rann_empty = 0;
                    end
                end
                % Completed generating Y for jj */
                %%%
            end % end of the loop for jj
            %
            % Generating G */
            for jj=1:2:gdimM1 % starts from 1
                while (1)
                    ranu1 = 2*rand(1,1)-1;
                    ranu2 = 2*rand(1,1)-1;
                    caps = ranu1*ranu1 + ranu2*ranu2;
                    if ((1>caps) && (0<caps))
                        break;
                    end
                end
                caps = sqrt(-2.0*log(caps)/caps);
                if 1==loc_rann_empty
                    G(jj) = ranu1*caps;
                    G(jj+1) = ranu2*caps;
                else
                    G(jj) = loc_rann_rest;
                    G(jj+1) = ranu1*caps;
                    loc_rann_rest = ranu2*caps;
                    loc_rann_empty=0;
                end
            end % End of the loop for jj
            % % %
            if 0==length(jj) % note this is for MATLAB
                jj=1;
            end
            % % %
            if ((jj+2)-1)~=gdim % note (jj+2)-1
                if 0==loc_rann_empty
                    G(gdim) = loc_rann_rest; % note gdim
                    loc_rann_empty = 1;
                else
                    while (1)
                        ranu1 = 2*rand(1,1)-1;
                        ranu2 = 2*rand(1,1)-1;
                        caps = ranu1*ranu1 + ranu2*ranu2;
                        if ((1>caps) && (0<caps))
                            break;
                        end
                    end
                    caps = sqrt(-2.0*log(caps)/caps);
                    G(gdim) = ranu1*caps; % note gdim
                    loc_rann_rest = ranu2*caps;
                    loc_rann_empty = 0;
                end
            end
            % Completed generatig G */
            %
            % Calculating tilZ */
            for jj=1:wdim % starts from 1
                for ll=1:jj-1 % starts from 1, but note ll<jj
                    tilZ(jj,ll)=0;
                    for kk=1:locMaxK % starts from 1
                        tilZ(jj,ll)= tilZ(jj,ll)+X_dv_k(jj,kk)*Y(ll,kk);
                    end
                end
                for ll=jj+1:wdim
                    tilZ(jj,ll-1)=0;
                    for kk=1:locMaxK % starts from 1
                        tilZ(jj,ll-1)= tilZ(jj,ll-1)+...
                            X_dv_k(jj,kk)*Y(ll,kk);
                    end
                end
            end
            % Completed calculating tilZ */
            %
            % Calculating xi_dv_sqstep */
            for jj=1:wdim % starts from 1
                xi_dv_sqstep(jj)=0;
                for ll=jj+1:wdim
                    %kk=jj*wdim-jj*(jj+1)/2+ll-jj-1;
                    kk=(jj-1)*wdim-(jj-1)*jj/2+ll-1-(jj); % note jj-1, ll-1 
                    xi_dv_sqstep(jj)=xi_dv_sqstep(jj)+...
                        loc_rani(ll)*G(kk+1); % note k+1
                end
                for ll=1:jj-1 % starts from 1, but note ll<jj
                    %kk=ll*wdim-ll*(ll+1)/2+jj-ll-1;
                    kk=(ll-1)*wdim-(ll-1)*ll/2+jj-1-(ll); % note jj-1, ll-1
                    xi_dv_sqstep(jj)=xi_dv_sqstep(jj)-...
                        loc_rani(ll)*G(kk+1); % note k+1 
                end
            end
            %
            % Calculating tilA/h */
            tmp=1+sqrt(1+rani_norm2);
            for jj=1:wdimM1 % starts from 1
                for ll=jj+1:wdim
                    %kk=jj*wdim-jj*(jj+1)/2+ll-jj-1;
                    kk=(jj-1)*wdim-(jj-1)*jj/2+ll-1-(jj); % note jj-1, ll-1
                    correction_term=... % note kk+1
                        locSqAn*(G(kk+1)+(loc_rani(ll)*xi_dv_sqstep(jj)-...
                        loc_rani(jj)*xi_dv_sqstep(ll))/tmp);
                    tilZ(jj,ll-1)=...
                        ((tilZ(jj,ll-1)-tilZ(ll,jj))/sq2+...
                        tilX(jj)*loc_rani(ll)-tilX(ll)*loc_rani(jj)+...
                        correction_term...
                        )/sq2pi;
                    tilZ(ll,jj)=-tilZ(jj,ll-1);
                end
            end
            % Completed calculating tilA/h */
            %
            % Calculating Transpose(Iij)/h */
            for ll=1:wdimM1 % starts from 1
                tmp_ibase_index=ibase_index+(ll-1)*wdimM1; % note ll-1
                for jj=ll+1:wdim
                    tranij(tmp_ibase_index+jj-1)=...
                        tranij(tmp_ibase_index+jj-1)+...
                        (randI(ibase+jj)*sq_rate_step+...
                        rani(ibase+jj)-loc_rani(jj))*loc_rani(ll)+...
                        loc_rani(jj)*loc_rani(ll)/2+tilZ(jj,ll);
                end
            end
            for ll=2:wdim % starts from 2
                tmp_ibase_index=ibase_index+(ll-1)*wdimM1; % note ll-1
                for jj=1:ll-1 % starts from 1, but note jj<ll
                    tranij(tmp_ibase_index+jj)=...
                        tranij(tmp_ibase_index+jj)+...
                        (randI(ibase+jj)*sq_rate_step+...
                        rani(ibase+jj)-loc_rani(jj))*loc_rani(ll)+...
                        loc_rani(jj)*loc_rani(ll)/2+tilZ(jj,ll-1);
                end
            end
            %
            if eps>=abs(step-(tmpStep+base_step))
                for ll=1:wdimM1 % starts from 1
                    tmp_ibase_index=ibase_index+(ll-1)*wdimM1; % note ll-1
                    for jj=ll+1:wdim
                        tranij(tmp_ibase_index+jj-1)=...
                            (tranij(tmp_ibase_index+jj-1)-...
                            randI(ibase+jj)*...
                            sq_rate_step*rani(ibase+ll))/rate_step;
                    end
                end
                for ll=2:wdim % starts from 2
                    tmp_ibase_index=ibase_index+(ll-1)*wdimM1; % note ll-1
                    for jj=1:ll-1 % starts from 1, but note jj<ll
                        tranij(tmp_ibase_index+jj)=...
                            (tranij(tmp_ibase_index+jj)-...
                            randI(ibase+jj)*...
                            sq_rate_step*rani(ibase+ll))/rate_step;
                    end
                end
                % Completed calculating Transpose(Iij)/h */
            end
            %
            %%%
            ibase_index=ibase_index+tranijDim;
        end % end of if for 0~=wdimM1
        %
        if eps>=abs(step-(tmpStep+base_step))
            for jj=1:wdim % starts from 1
                ran_diag(ibase+jj)=(ran_diag(ibase+jj)-...
                    randI(ibase+jj)*...
                    sq_rate_step*rani(ibase+jj))/rate_step;
                rani(ibase+jj)=rani(ibase+jj)/sq_rate_step;
            end
            % Completed generating delW/sqrt(h) and Iii/h. */
        end
        %%%
        ibase=ibase+wdim;
    end % end of the foop for ii
    %%%
    tmpStep =tmpStep+base_step; % Update on tmpStep */
    
end % end of while for tmpStep


