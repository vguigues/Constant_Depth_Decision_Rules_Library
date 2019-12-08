
%Nb(m): number of thermal units in subsystem m
%Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
%Tendancy_Inflows(t,m) is ct(m) for t=1,...,TabS(m),m=1,...,NS.
%sigmaInflows(t,m): s.d. of \xi_t(m)
%Probabilities(t,j): probability of realization Inflow_Noises(t,:,j) for period t+1,t=1,...,T-1,scenario j, sum(Probabilities(t,:))=1.

function [Time,Zsups,Zinfs,Iter,DynamicSubi,DynamicSubj,DynamicValij,Alphas,Thetas]=sddp_hydro_interstage_dependent(Nb,talpha,T,S,NS,N,gamma,M,Inflow_Noises,Probabilities,Demand,Thermal_Costs,vminT,vmaxT,vminH,vmaxH,xmax,x0,lambda,evap,phi,TabS,TabP,tol,Tendancy_Inflows,Recent_History_Inflows,sigmaInflows,Itermax)


zsup=inf;
zinf=-inf;
Iter=1;

Zsups=[];
Zinfs=[];

Alphas=cell(1,T-1);
Thetas=cell(1,T-1);

FixedSubi=[];
FixedSubj=[];
FixedValij=[];

Nb_Thermal=sum(Nb);

for i=1:NS
    FixedSubi=[FixedSubi,i,i,i];
    FixedSubj=[FixedSubj,Nb_Thermal+i,Nb_Thermal+NS+i,Nb_Thermal+2*NS+i];
    FixedValij=[FixedValij,1,1,1];
end

for i=1:NS
    FixedSubi=[FixedSubi,(NS+i)*ones(1,Nb(i)+1)];
    FixedSubj=[FixedSubj,Nb_Thermal+NS+i,[sum(Nb(1:i-1))+1:1:sum(Nb(1:i))]];
    FixedValij=[FixedValij,ones(1,Nb(i)+1)];
end

DynamicSubi=cell(1,T-1);
DynamicSubj=cell(1,T-1);
DynamicValij=cell(1,T-1);

taille=zeros(T,1);
StateOrder=zeros(T,NS);
for t=1:T
    taille(T-t+1)=NS;
    for m=1:NS
        periode(T-t+1,m)=mod(T-t+1,TabS(m));
        if (periode(T-t+1,m)==0)
            periode(T-t+1,m)=TabS(m);
        end
        if (t==1)
            StateOrder(T,m)=TabP{1,m}(periode(T,m));
        else
            StateOrder(T-t+1,m)=max(TabP{1,m}(periode(T-t+1,m)),StateOrder(T-t+2,m)-1);
        end
        taille(T-t+1)=taille(T-t+1)+StateOrder(T-t+1,m);
    end
end

Cum_Probas=zeros(T-1,M+1);
for t=1:T-1
    Cum_Probas(t,:)=[0,cumsum(Probabilities(t,:))];
end

Iter=1;
End_Algo=1;
Costs=zeros(N,1);

tic

while End_Algo
    %Iter
    %Forward pass
    Trial_States=zeros(T-1,NS,S);
    Sampled_Inflows=zeros(T,NS,S);
    for s=1:S
        ct=0;
        for t=1:T
            if (t==1)
                Sampled_Inflows(1,:,s)=Recent_History_Inflows(1,:);
            else
                Alea_Uniform=rand;
                [~,Index] = histc(Alea_Uniform,Cum_Probas(t-1,:));
                if (Alea_Uniform==1)
                    Sampled_Noise=Inflow_Noises(t-1,:,M);
                else
                    Sampled_Noise=Inflow_Noises(t-1,:,Index);
                end
                for m=1:NS
                    Sampled_Inflows(t,m,s)=Tendancy_Inflows(periode(t,m),m);
                    for j=1:TabP{1,m}(periode(t,m))
                        if ((t-j)<=1)
                            Sampled_Inflows(t,m,s)=Sampled_Inflows(t,m,s)+phi{1,m}{1,periode(t,m)}(j)*Recent_History_Inflows(-t+j+2,m);
                        else
                            Sampled_Inflows(t,m,s)=Sampled_Inflows(t,m,s)+phi{1,m}{1,periode(t,m)}(j)*Sampled_Inflows(t-j,m,s);
                        end
                    end
                    Sampled_Inflows(t,m,s)=Sampled_Inflows(t,m,s)+sigmaInflows(periode(t,m),m)*Sampled_Noise(m);
                    if (Sampled_Inflows(t,m,s)<0)
                        disp('Negative inflows');
                        pause
                    end
                end
            end
            
            blx=[];
            bux=[];
            for m=1:NS
                blx=[blx;vminT{t,m}];
                bux=[bux;vmaxT{t,m}];
            end
            blx=[blx;(lambda(t,:).*xmax(t,:))';vminH(t,:)';zeros(NS,1)];
            bux=[bux;xmax(t,:)';vmaxH(t,:)';inf*ones(NS,1)];
            if ((t<T)&&(Iter>1))
                blx=[blx;-inf];
                bux=[bux;inf];
            end
            
            
            if ((t==T)||(Iter==1))
                nbvar=Nb_Thermal+3*NS;
                c=zeros(nbvar,1);
                blc=zeros(2*NS,1);
                buc=inf*ones(2*NS,1);
            else
                nbvar=1+Nb_Thermal+3*NS;
                c=zeros(nbvar,1);
                c(nbvar)=1;
                blc=zeros(2*NS+S*(Iter-1),1);
                buc=inf*ones(2*NS+S*(Iter-1),1);
            end
            
            if ((t<T)&&(Iter>1))
                subi=[FixedSubi,DynamicSubi{1,t}];
                subj=[FixedSubj,DynamicSubj{1,t}];
                valij=[FixedValij,DynamicValij{1,t}];
            else
                subi=[FixedSubi];
                subj=[FixedSubj];
                valij=[FixedValij];
            end
            for m=1:NS
                
                c(sum(Nb(1:(m-1)))+1:sum(Nb(1:m)))=Thermal_Costs{t,m};
                
                if (t>1)
                    blc(m)=gamma(t,m)*Sampled_Inflows(t,m,s)+Trial_States(t-1,m,s)-evap(t,m);
                    buc(m)=gamma(t,m)*Sampled_Inflows(t,m,s)+Trial_States(t-1,m,s)-evap(t,m);
                else
                    blc(m)=gamma(t,m)*Sampled_Inflows(t,m,s)+x0(m)-evap(t,m);
                    buc(m)=gamma(t,m)*Sampled_Inflows(t,m,s)+x0(m)-evap(t,m);
                end
                
                blc(m+NS)=-(1-gamma(t,m))*Sampled_Inflows(t,m,s)+Demand(t,m);
                
            end
            
            if (t<T)
                for kp=1:(Iter-1)
                    for sp=1:S
                        buc(2*NS+(kp-1)*S+sp)=inf;
                        blc(2*NS+(kp-1)*S+sp)=Thetas{1,t}((kp-1)*S+sp);
                        compteur=1;
                        for m=1:NS
                            for cp=1:StateOrder(t+1,m)
                                if (t+1-cp>=1)
                                    blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Sampled_Inflows(t+1-cp,m,s);
                                else
                                    blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Recent_History_Inflows(-t+1+cp,m);
                                end
                            end
                            compteur=compteur+StateOrder(t+1,m);
                        end
                    end
                end
            end
            
            clear prob;
            prob.c=c;
            prob.blc=blc;
            prob.buc=buc;
            prob.blx=blx;
            prob.bux=bux;
            if ((t==T)||(Iter==1))
                prob.a=sparse(subi,subj,valij,2*NS,nbvar);
            else
                prob.a=sparse(subi,subj,valij,2*NS+(Iter-1)*S,nbvar);
            end
            [r,res]=mosekopt('minimize echo(0) statuskeys(1)',prob);
            if (res.sol.itr.solsta==5)
                disp('Unfeasible problem in forward pass');
                pause
            end
            sol=res.sol.bas.xx;
            if (t<T)
                Trial_States(t,:,s)=sol(1+Nb_Thermal:Nb_Thermal+NS)';
                ct=ct+sol'*c-sol(nbvar);
            else
                ct=ct+sol'*c;
            end
        end %for t
        if ((Iter-1)*S+s<=N)
            Costs((Iter-1)*S+s)=ct;
        else
            Costs=[Costs(2:N);ct];
        end
    end %for s
    
    if (Iter*S>=N)
        Mean_Cost=mean(Costs);
        Sigma_Cost=sqrt(var(Costs));
        zsup=Mean_Cost+talpha*Sigma_Cost/sqrt(N);
        Zsups=[Zsups,zsup];
    end
    
    %Backward phase of SDDP
    for ta=1:T-1
        t=T-ta+1;
        for s=1:S
            Thetas{1,t-1}((Iter-1)*S+s)=0;
            beta=zeros(1,NS);
            Alphas{1,t-1}=[Alphas{1,t-1};zeros(1,taille(t))];
            for j=1:M
                if (t==1)
                    Inflowst=Recent_History_Inflows(1,:);
                else
                    Sampled_Noise=Inflow_Noises(t-1,:,j);
                    for m=1:NS
                        Inflowst(m)=Tendancy_Inflows(periode(t,m),m);
                        for p=1:TabP{1,m}(periode(t,m))
                            if ((t-p)<=1)
                                Inflowst(m)=Inflowst(m)+phi{1,m}{1,periode(t,m)}(p)*Recent_History_Inflows(-t+p+2,m);
                            else
                                Inflowst(m)=Inflowst(m)+phi{1,m}{1,periode(t,m)}(p)*Sampled_Inflows(t-p,m,s);
                            end
                        end
                        Inflowst(m)=Inflowst(m)+sigmaInflows(periode(t,m),m)*Sampled_Noise(m);
                        if (Inflowst(m)<0)
                            disp('Negative inflows');
                            pause
                        end
                    end
                end
                
                blx=[];
                bux=[];
                for m=1:NS
                    blx=[blx;vminT{t,m}];
                    bux=[bux;vmaxT{t,m}];
                end
                blx=[blx;(lambda(t,:).*xmax(t,:))';vminH(t,:)';zeros(NS,1)];
                bux=[bux;xmax(t,:)';vmaxH(t,:)';inf*ones(NS,1)];
                if (t<T)
                    blx=[blx;-inf];
                    bux=[bux;inf];
                end
                
                
                if (t==T)
                    nbvar=Nb_Thermal+3*NS;
                    c=zeros(nbvar,1);
                    blc=zeros(2*NS,1);
                    buc=inf*ones(2*NS,1);
                else
                    nbvar=1+Nb_Thermal+3*NS;
                    c=zeros(nbvar,1);
                    c(nbvar)=1;
                    blc=zeros(2*NS+S*Iter,1);
                    buc=inf*ones(2*NS+S*Iter,1);
                end
                
                if (t==T)
                    subi=[FixedSubi];
                    subj=[FixedSubj];
                    valij=[FixedValij];
                else
                    subi=[FixedSubi,DynamicSubi{1,t}];
                    subj=[FixedSubj,DynamicSubj{1,t}];
                    valij=[FixedValij,DynamicValij{1,t}];
                end
                
                for m=1:NS
                    
                    c(sum(Nb(1:(m-1)))+1:sum(Nb(1:m)))=Thermal_Costs{t,m};
                                    
                    if (t>1)
                        blc(m)=gamma(t,m)*Inflowst(m)+Trial_States(t-1,m,s)-evap(t,m);
                        buc(m)=gamma(t,m)*Inflowst(m)+Trial_States(t-1,m,s)-evap(t,m);
                    else
                        blc(m)=gamma(t,m)*Inflowst(m)+x0(m)-evap(t,m);
                        buc(m)=gamma(t,m)*Inflowst(m)+x0(m)-evap(t,m);
                    end
                    
                    blc(m+NS)=-(1-gamma(t,m))*Inflowst(m)+Demand(t,m);
                    
                end
                
                if (t<T)
                    for kp=1:Iter
                        for sp=1:S
                            blc(2*NS+(kp-1)*S+sp)=Thetas{1,t}((kp-1)*S+sp);
                            compteur=1;
                            for m=1:NS
                                for cp=1:StateOrder(t+1,m)
                                    if (t+1-cp<=1)
                                        blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Recent_History_Inflows(-t+1+cp,m);
                                    elseif (cp==1)
                                        blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Inflowst(m);
                                    else
                                        blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Sampled_Inflows(t+1-cp,m,s);
                                    end
                                end
                                compteur=compteur+StateOrder(t+1,m);
                            end
                        end
                    end
                end
                
                clear prob;
                prob.c=c;
                prob.blc=blc;
                prob.buc=buc;
                prob.blx=blx;
                prob.bux=bux;
                if (t==T)
                    prob.a=sparse(subi,subj,valij,2*NS,nbvar);
                else
                    prob.a=sparse(subi,subj,valij,2*NS+Iter*S,nbvar);
                end
                [r,res]=mosekopt('minimize echo(0) statuskeys(1)',prob);
                sol=res.sol.bas.xx;
                solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
                if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
                    disp('Unfeasible primal problem in backward pass');
                    pause
                elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
                    disp('Primal infinite optimal value in backward pass');
                    pause
                end
                
                ThetaAux=sol'*c;
                dual=res.sol.bas.slc;
                dualb=res.sol.bas.suc;
                lambda1=dual(1:NS)-dualb(1:NS);
                lambda2=dual(NS+1:2*NS);
                ThetaAux=ThetaAux-lambda1'*Trial_States(t-1,:,s)';
                beta=beta+Probabilities(t-1,j)*lambda1';
                compteur=1;
                compteurbis=1;
                AlphasAux=zeros(1,taille(t));
                for m=1:NS
                    Aux=(lambda1(m)*gamma(t,m)-lambda2(m)*(1-gamma(t,m)))*(phi{1,m}{1,periode(t,m)});
                    AlphasAux(compteur:compteur+TabP{1,m}(periode(t,m))-1)=Aux;
                    if (t<T)
                        for kp=1:Iter
                            for sp=1:S
                                AlphasAux(compteur:(compteur+TabP{1,m}(periode(t,m))-1))=AlphasAux(compteur:(compteur+TabP{1,m}(periode(t,m))-1))+Alphas{1,t}((kp-1)*S+sp,compteurbis)*(phi{1,m}{1,periode(t,m)})*(dual(2*NS+(kp-1)*S+sp));
                                for w=2:StateOrder(t+1,m)
                                    AlphasAux(compteur+w-2)=AlphasAux(compteur+w-2)+Alphas{1,t}((kp-1)*S+sp,compteurbis+w-1)*(dual(2*NS+(kp-1)*S+sp));
                                end
                            end
                        end
                        compteurbis=compteurbis+StateOrder(t+1,m);
                    end
                    for Ind=1:StateOrder(t,m)
                        if (t-Ind>=1)
                           ThetaAux=ThetaAux-AlphasAux(compteur+Ind-1)*Sampled_Inflows(t-Ind,m,s);
                        else
                           ThetaAux=ThetaAux-AlphasAux(compteur+Ind-1)*Recent_History_Inflows(-t+Ind+2,m); 
                        end
                    end
                    compteur=compteur+StateOrder(t,m);
                end
                Thetas{1,t-1}((Iter-1)*S+s)=Thetas{1,t-1}((Iter-1)*S+s)+ThetaAux*Probabilities(t-1,j);
                Alphas{1,t-1}((Iter-1)*S+s,:)=Alphas{1,t-1}((Iter-1)*S+s,:)+Probabilities(t-1,j)*AlphasAux;
            end %for j
            DynamicSubi{1,t-1}=[DynamicSubi{1,t-1},(2*NS+(Iter-1)*S+s)*ones(1,NS+1)];
            DynamicSubj{1,t-1}=[DynamicSubj{1,t-1},Nb_Thermal+3*NS+1,[Nb_Thermal+1:1:Nb_Thermal+NS]];
            DynamicValij{1,t-1}=[DynamicValij{1,t-1},1,-beta];
        end %for s
    end %for ta
    %zinf
    
    t=1;
    Inflowst=Recent_History_Inflows(1,:);
    
    blx=[];
    bux=[];
    for m=1:NS
        blx=[blx;vminT{1,m}];
        bux=[bux;vmaxT{1,m}];
    end
    blx=[blx;(lambda(1,:).*xmax(1,:))';vminH(1,:)';zeros(NS,1)];
    bux=[bux;xmax(1,:)';vmaxH(1,:)';inf*ones(NS,1)];
    blx=[blx;-inf];
    bux=[bux;inf];
    
    nbvar=1+Nb_Thermal+3*NS;
    c=zeros(nbvar,1);
    c(nbvar)=1;
    blc=zeros(2*NS+S*Iter,1);
    buc=inf*ones(2*NS+S*Iter,1);
    
    subi=[FixedSubi,DynamicSubi{1,1}];
    subj=[FixedSubj,DynamicSubj{1,1}];
    valij=[FixedValij,DynamicValij{1,1}];
    
    for m=1:NS
        
        c(sum(Nb(1:(m-1)))+1:sum(Nb(1:m)))=Thermal_Costs{1,m};
        
        c(Nb_Thermal+2*NS+m)=0;
        c(Nb_Thermal+NS+m)=0;
        
        blc(m)=gamma(1,m)*Inflowst(m)+x0(m)-evap(1,m);
        buc(m)=gamma(1,m)*Inflowst(m)+x0(m)-evap(1,m);
        
        blc(m+NS)=-(1-gamma(1,m))*Inflowst(m)+Demand(1,m);
        
    end
    
    for kp=1:Iter
        for sp=1:S
            blc(2*NS+(kp-1)*S+sp)=Thetas{1,t}((kp-1)*S+sp);
            compteur=1;
            for m=1:NS
                for cp=1:StateOrder(t+1,m)
                    blc(2*NS+(kp-1)*S+sp)=blc(2*NS+(kp-1)*S+sp)+Alphas{1,t}((kp-1)*S+sp,compteur+cp-1)*Recent_History_Inflows(-t+1+cp,m);
                end
                compteur=compteur+StateOrder(t+1,m);
            end
        end
    end
    
    clear prob;
    prob.c=c;
    prob.blc=blc;
    prob.buc=buc;
    prob.blx=blx;
    prob.bux=bux;
    prob.a=sparse(subi,subj,valij,2*NS+Iter*S,nbvar);
    [r,res]=mosekopt('minimize echo(0) statuskeys(1)',prob);
    sol=res.sol.bas.xx;
    solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
    if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
        disp('Unfeasible first stage problem');
        pause
    elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
        disp('Primal infinite optimal value');
        pause
    end
    
    %Update zinf
    
    zinf=c'*sol;
    Zinfs=[Zinfs,zinf];
    %zinf
    %zsup
    if (Iter==Itermax)
        End_Algo=0;
    elseif ((Iter*S>=N)&&(abs((zsup-zinf)/zsup)<=tol))
       End_Algo=0;
    end
    Iter=Iter+1;
end
Iter=Iter-1;
Time=toc;