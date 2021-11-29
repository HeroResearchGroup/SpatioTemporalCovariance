function [L,S,u,sT,v] = robust_PCA(X,lambdaL, lambdaS,method,toep_flg)
%ROBUST_PCA Does robust PCA as in Nadakuditi et al 2014 (SSP). Proximal
%gradient method.
%
%Author: Kristjan Greenewald
%Date: 3/31/14
%

thresh = 1e-4;
%Step size
tau = .3;%.1%Find out how to set this!!!
r = lambdaL; %Need to choose for OptShrink.

c = size(X,2)/size(X,1);
q = min(size(X));
S = sparse(zeros(size(X)));
M = X; %Initialize.
err = inf;
Mprev = zeros(size(M));
while err > thresh %not converged
    %SVT
    
    switch method
        case 'SVT'
            sT = svd(M - S,'econ');
    
            %sT = diag(s);
            sT = sT - lambdaL*tau;
            
            sT(sT<0) = 0;
            r = find(sT>0,1,'last');
            if isempty(r)
                r = 1;
            end
            sT = sT(1:r);
            if sT(1) == 0
                sT(1) = 1e-10;
            end
            [u,~,v] = svds(M-S,max(1,r));
        case 'OptShrink'
            sT = svd(M - S,'econ');
            [u,~,v] = svds(M-S,r);
            %sT = diag(s);
            for i = 1:r
                sumV = 1/(q-r)*sum(sT(i)./(sT(i)^2-sT(r+1:q).^2));
                D = sumV*(c*sumV+(1-c)/sT(i));
                sumVprime = 1/(q-r)*sum((-sT(i)^2-sT(r+1:q).^2)./(sT(i)^2-sT(r+1:q).^2).^2);
                Dprime = (c*sumV+(1-c)/sT(i))*(sumVprime) + (c*sumVprime-(1-c)/sT(i)^2)*sumV;
                sT(i) = -2*D/Dprime;
            end
            %sT(r+1:end) = 0;
            sT = sT(1:r);
            
    end
    %sT(sT<0) = 0;
    L = u*diag(sparse(sT))*v';
    %Soft
    Sabs = abs(M-L)-lambdaS*tau;
    if toep_flg
        int = (size(X,1)+1)/2;
        for i = -int+1:int-1
            Sabs(i+int,:) = abs(M(i+int,:)-L(i+int,:))-lambdaS/sqrt(int-abs(i))*tau;
        end
    end
    Sabs(Sabs<0) = 0;
    Sabs = sparse(Sabs);
    Ssign = sign(M-L);
    S = sparse(Ssign.*Sabs);
    
    %Update M.
    M = L + S - tau*(L+S-X);
    
    
    err = sqrt(sum(sum((M-Mprev).^2))/sum(sum(M.^2)));
    Mprev = M;
end

