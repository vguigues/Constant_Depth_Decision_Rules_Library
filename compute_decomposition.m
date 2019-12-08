

%For the model Z_t(m)=\sum_{j=1}^{Lag{1,m}(t)} phinp{1,m}{1,t}(j) Z_{t-j}(m)+eta_t(m)
%we first compute alpha, beta, and theta such that, for a given t,
%Z_{t+j}(m)=\sum_{ell=0}^{Lagmax(m)} alpha{el1,m}(j+Lagmax(m)+1,ell+1)Z_{t-ell}(m)
%+\sum_{ell=1}^{j} beta{m,j}(ell) eta_{t+l}(m)
%We then compute (output of the function) the coefficients alpha, beta, and
%theta such that
%\xi_{t+j}(m)=\theta(j,m)+\sum_{l=0}^{Lagmax(m)} alpha{1,m}(j,l+1)\xi_{t-l}(m)
%+\sum_{l=1}^{j} beta{m,j}(l) eta_{t+l}(m)
%To be called with t=1 if observations are available until time 1.

function [alpha,beta,theta]=compute_decomposition(mus,sigmas,phinp,M,Lags,t,T,Lagmax)

alpha=cell(1,M);
beta=cell(M,T-t);
theta=zeros(T-t,M);

for i=1:M
    alpha{1,i}=zeros(T-t+Lagmax(i)+1,Lagmax(i)+1);
    for p=1:T-t
        beta{i,p}=zeros(1,p);
    end
end

for i=1:M
    beta{i,1}(1)=1;
    for p=0:Lagmax(i) 
        for k=0:Lagmax(i)
	        if (k==p)
	            alpha{1,i}(Lagmax(i)-p+1,k+1)=1;
            else
 	            alpha{1,i}(Lagmax(i)-p+1,k+1)=0;
            end
        end
    end
end

for i=1:M
    for p=1:T-t
        for k=0:Lagmax(i)
	        for j=1:Lags{1,i}(t+p)
                
                alpha{1,i}(p+Lagmax(i)+1,k+1)=alpha{1,i}(p+Lagmax(i)+1,k+1)+(phinp{1,i}{1,t+p}(j))*(alpha{1,i}(p+Lagmax(i)+1-j,k+1));
            end
        end
        for k=1:p-1
            for j=1:min(Lags{1,i}(t+p),p-k)
			    beta{i,p}(k)=beta{i,p}(k)+(phinp{1,i}{1,t+p}(j))*(beta{i,p-j}(k));
            end
        end
        beta{i,p}(p)=1;
    end
end

for i=1:M
    for p=1:T-t
        theta(p,i)=mus(t+p,i);
        for k=0:Lagmax(i)
            if (t-k<=0)
            alpha{1,i}(p,k+1)=(alpha{1,i}(p+Lagmax(i)+1,k+1))*(sigmas(t+p,i))/(sigmas(t-k+12,i));        
            theta(p,i)=theta(p,i)-(alpha{1,i}(p,k+1))*(mus(t-k+12,i));
            else
            alpha{1,i}(p,k+1)=(alpha{1,i}(p+Lagmax(i)+1,k+1))*(sigmas(t+p,i))/(sigmas(t-k,i));        
            theta(p,i)=theta(p,i)-(alpha{1,i}(p,k+1))*(mus(t-k,i));    
            end
        end
        for k=1:p
            beta{i,p}(k)=(beta{i,p}(k))*(sigmas(t+p,i));
        end
    end
end
