


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix A^{t tau} is given in sparse format by 
%subi_a{1,t}{1,tau}, subj_a{1,t}{1,tau}, valij_a{1,t}{1,tau}. 
%beta{1,t}{1,s}(ell,i) is beta^{t i}_{s scenario(ell)} of the paper where 
%scenario(ell) is the scenario numbered ell.
%cost{1,t} is the cost vector for stage t
%probabilities(t-1,k) is the probability of kth realization for stage t=2,..,T.
%ds(t) is the number of realizations for stage t
%T is the number of stages
%A^{t tau} has size ps(t)*qs(t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sol: optimal solution 
%opt_value: optimal value
%out=0: optimal solution found; out=1: infeasible; out=2: primal infinite
%opt value
%time1: time to store the data for the optimization problem
%time2: time to solve the optimization problem
%nvars: number of variables of the problem
%counter: number of constraints of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sol,opt_value,out,time1,time2,nbvars,counter]=solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)

subi=[];
subj=[];
valij=[];
blc=[];
buc=[];

tic

ds=[ones(depth-1,1);ds];

productsu=zeros(T,depth);
for s=1:T
    productsu(s,depth)=ds(s+depth-1);
    for i=1:(depth-1)
        productsu(s,depth-i)=productsu(s,depth-i+1)*ds(s+depth-i-1);
    end
end

productsz=zeros(T,depth-1);
for s=1:T
    productsz(s,depth-1)=ds(s+depth-2);
    for i=1:(depth-2)
        productsz(s,depth-1-i)=productsz(s,depth-i)*ds(s+depth-2-i);
    end
end

Iu=cell(1,T);
Iv=cell(1,T);
Iz=cell(1,T);

nbvars=0;

for t=1:T
    Iv{1,t}=cell(1,t);
    for s=1:t
        Iv{1,t}{1,s}=zeros(productsu(s,1),ps(t));
        for ell=1:productsu(s,1)
            Iv{1,t}{1,s}(ell,:)=[nbvars+1:nbvars+ps(t)];
            nbvars=nbvars+ps(t);
        end
    end
end


for t=1:T
    Iu{1,t}=cell(1,t);
    for s=1:t
        Iu{1,t}{1,s}=zeros(productsu(s,1),qs(t));
        for ell=1:productsu(s,1)
            Iu{1,t}{1,s}(ell,:)=[nbvars+1:nbvars+qs(t)];
            nbvars=nbvars+qs(t);
        end
    end
end


for t=2:T
    Iz{1,t}=cell(1,t);
    for s=2:t
        Iz{1,t}{1,s}=zeros(productsz(s,1),ps(t));
        for ell=1:productsz(s,1)
            Iz{1,t}{1,s}(ell,:)=[nbvars+1:nbvars+ps(t)];
            nbvars=nbvars+ps(t);
        end
    end
end

blx=-inf*ones(nbvars,1);
bux=inf*ones(nbvars,1);
bux(Iv{1,1}{1,1}(1,:))=zeros(ps(1),1);

counter=1;
for t=1:T
    for s=1:t
        for i=1:ps(t)
            for ell=1:productsu(s,1)
                %We need the index of v_{s, ell}^{t i}    
                blc=[blc;-betas{1,t}{1,s}(ell,i)];
                buc=[buc;-betas{1,t}{1,s}(ell,i)];
     
                subi=[subi,counter];
                subj=[subj,Iv{1,t}{1,s}(ell,i)];
                valij=[valij,1];
                
                for tau=s:t
                    [index]=find(subi_a{1,t}{1,tau}==i);
                    for k=1:length(index)
                        subi=[subi,counter];
                        subj=[subj,Iu{1,tau}{1,s}(ell,subj_a{1,t}{1,tau}(index(k)))];
                        valij=[valij,-valij_a{1,t}{1,tau}(index(k))];
                    end
                end
                counter=counter+1;
            end
        end
    end
end


for t=2:T
    for i=1:ps(t)
        for ell=1:productsu(t)
            [etas]=index_to_scenario(ell,productsu(t,:),depth);
            ellz=scenario_to_index(etas(1:depth-1),productsz(t,:)',depth-1);
            
            subi=[subi,counter];
            subj=[subj,Iz{1,t}{1,t}(ellz,i)];
            valij=[valij,1];
            
            subi=[subi,counter];
            subj=[subj,Iv{1,t}{1,t}(ell,i)];
            valij=[valij,-1];
            
            blc=[blc;0];
            buc=[buc;inf];
            
            counter=counter+1;
        end
    end
end

for t=3:T
    for s=3:t
        for i=1:ps(t)
            for ell=1:productsu(s-1)
                [etas]=index_to_scenario(ell,productsu(s-1,:),depth);
                etasleft=etas(1:depth-1);
                etasright=etas(2:depth);
                  
                ellzl=scenario_to_index(etasleft,productsz(s-1,:)',depth-1);
                ellzr=scenario_to_index(etasright,productsz(s,:)',depth-1);
                 
                subi=[subi,counter];
                subj=[subj,Iz{1,t}{1,s-1}(ellzl,i)];
                valij=[valij,1];
                    
                subi=[subi,counter];
                subj=[subj,Iz{1,t}{1,s}(ellzr,i)];
                valij=[valij,-1];
                   
                subi=[subi,counter];
                subj=[subj,Iv{1,t}{1,s-1}(ell,i)];
                valij=[valij,-1];
                
                blc=[blc;0];
                buc=[buc;inf];
                
                counter=counter+1; 
            end
        end
    end
end

for t=2:T
    for i=1:ps(t)
        subi=[subi,counter];
        subj=[subj,Iz{1,t}{1,2}(1,i)];
        valij=[valij,1];
        subi=[subi,counter];
        subj=[subj,Iv{1,t}{1,1}(1,i)];
        valij=[valij,1];
        blc=[blc;-inf];
        buc=[buc;0];
        counter=counter+1;
    end
end

prob.blx=blx;
prob.bux=bux;
prob.blc=blc;
prob.buc=buc;
prob.a = sparse(subi,subj,valij,counter-1,nbvars);
obj=zeros(nbvars,1);
for t=1:T
    for s=1:t
        for ell=1:productsu(s,1)
            [etas]=index_to_scenario(ell,productsu(s,:),depth);
            aux=1;
            for k=max(1,depth-s+2):depth
                  aux=aux*probabilities(s-depth+k-1,etas(k));
            end
            for i=1:qs(t)
                obj(Iu{1,t}{1,s}(ell,i))=aux*cost{1,t}(i);
            end
        end
    end
end
prob.c=obj;
time1=toc;

tic
[~,res]=mosekopt('minimize echo(0)',prob);
sol=res.sol.bas.xx;      
solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
    disp('Unfeasible problem');
    out=1;
    opt_value=inf;
elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
    disp('Primal infinite optimal value');
    out=2;
    opt_value=-inf;
else
    out=0;
    opt_value=prob.c'*sol;
end
time2=toc;
counter=counter-1;