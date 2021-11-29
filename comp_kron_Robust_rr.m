function [T,S,kron_mat,err_full] = comp_kron_Robust_rr(corrmat,m,int,lambdaL,lambdaS,toep_flg,optSVD)
%Function to do LS Robust PCA-based Kronecker decomposition of empirical
%Sample covariance matrix: corrmat = T (x) S.
%m = S factor dimension, int = time length
%T: Time covariance, S: spatial.
%num_comp: nuclear norm penalty. 
%LambdaS sparsity penalization.
%kron(T,S) is the Kronecker component of the robust estimate
%OptSVD: Controls which RobustPCA method to use (see robust_PCA.m)
%
%Options:
%
%
%toep_flg: if nonzero, constrain estimate to be block Toeplitz (for stationary temporal
%processes).
%
%
%
%NOTE: After this, due to numerics may need to do additional regularization such as
%projection onto the set of positive semidefinite matrices.
%
%Author: Kristjan Greenewald
%Date: 12/31/12
%3/30/15
%
nonZeroVec = NaN;

n = 1;
num_comp = lambdaL;
basis_flg{1} = 0;

%Kronecker
%Rearrange
if isnan(nonZeroVec)
    Sflag = 1;
fun = @(x)x.data(:)';
Sig1 = blockproc(corrmat,[m*n m*n],fun);
Sig = reshape(Sig1',m^2*n^2,int^2)';
nonZeroVec = 1:m^2*n^2;
else
    Sflag = 0;
    Sig = corrmat;
end

%Remove the zero columns.
%Sig = Sig(:,nonZeroVec);



% if length(diag_flg) == 2 && diag_flg(2) == 1 && diag_flg(1)~= 1
%     indmt = reshape(1:int^2,int,int);
%     indmtt = triu(indmt);
%     indvc = indmtt(:);
%     indvc = indvc(indvc ~= 0);
%     Sig = sqrt(2)*Sig(indvc,:);
%     Sig(diag(indmt),:) = 1/sqrt(2)*Sig(diag(indmt),:);
% else
    
    
% end

    
%else %Standard Kronecker.
    %%
    if basis_flg{1}
            Tmat = basis_flg{2}(:);
    end
    if toep_flg
        
        %PARAMETER: Length of toeplitz extension (determines weights). Increase if
        %doing block banded estimation.
        tlength = int;
        
        
        %Combine rows(temporal) based on tlength length Toeplitz constraint. Hence scaling.
        SigToep = zeros(2*int-1,size(Sig,2));
        
        %Get temporal diag row numbers.
        
        
        for i = -int+1:int-1%Loop through temporal diagonals.
            %Find diag locations            
            indi = (i <= 0)*(abs(i)+1) + (i > 0)*(1+i*int);
            for ii = 1:int-abs(i)-1
                indi = [indi indi(end)+1+int];
            end
%             if indi(end) > int^2;
%                 indi = indi(1:end-1);
%             end
            
            %Combine rows
            SigToep(i+int,:) = (tlength-abs(i))/(int-abs(i))*1/sqrt(tlength-abs(i))*sum(Sig(indi,:),1);%*(int-abs(i));
            
            if basis_flg{1}
                U(i+int,1) = (tlength-abs(i))/(int-abs(i))*1/sqrt(tlength-abs(i))*sum(Tmat(indi));%*(int-abs(i));
            end
        end
        Sig = SigToep;
    else
        if basis_flg{1}
            U = Tmat;
        end
    end
    %SVD% Low rank approx of Sig.

    if basis_flg{1} == 0 %Usual
        %[U,s,V] = svds(Sig,num_comp);%8.3e5
        lambdaL = num_comp;
        if Sflag
            [~,Sp,U,s,V] = robust_PCA(Sig,lambdaL, lambdaS,optSVD,toep_flg);
        else
        [~,Sp,U,s,V] = robust_PCA(Sig,lambdaL, lambdaS,optSVD,toep_flg); %optSVD: 'SVT','OptShrink'
        end
        s = diag((s));
        
    else %Predetermined T vector.
        s = diag([1 zeros(1,size(Sig,1))]);
        U = U/norm(U,2);
        %numcomp = 1 for guaranteed psdness.
        %if toep_flg
            V = (U'*Sig)';
        %end
    end
     
%      if nuke_flg(1) %PRLS nuclear norm regularization.
%          diagS = diag(s) - nuke_flg(2)/2;
%          diagS(diagS < 0) = 0;
%          
%          if size(s,1) > size(s,2)
%             s = [diag(diagS); zeros(size(s,1)-size(s,2),size(s,2))];
%          elseif size(s,2) > size(s,1)
%              s = [diag(diagS), zeros(size(s,1),size(s,2)-size(s,1))];
%          end
%          
%      end
%     
%     
%     
%     %Get kronecker components
%     uvec = zeros(1,int^2);
%     vvec = zeros(1,(m*n)^2);
%     uvec(indvc) = U(:,1)';
%     vvec(:) = V(:,1)';
%     T = reshape(uvec,int,int)*sqrt(s(1,1));
%     S = reshape(vvec,m*n,m*n)*sqrt(s(1,1));
%     
%     err = (sum(diag(s).^2) - s(1,1)^2)/s(1,1)^2;
%     
%     kron_mat = kron(T,S);
    
    err_full = 0;
    %if num_comp > 1
%         K1 = blockproc(kron_mat,[m*n m*n],fun);
%         Kmat= reshape(K1',(m*n)^2,int^2)';
%         %SVD again
%         [U,s,V] = svds(Sig-Kmat,num_comp-1);
        

        %Form kron_mat
        if Sflag
            kron_mat = zeros(size(corrmat));
        else
            
        kron_mat = sparse(size(corrmat));
        end
        %Sflag = 0;
        for i = 1:nnz(s)
            
            %uvec = zeros(1,int^2);
            vvec = sparse(1,(m*n)^2);
            indvc = 1:int^2;
            uvec = U(:,i)';
            if Sflag
                vvec(nonZeroVec) = V(:,i)';
            end
            
            if toep_flg
                tv = uvec;%uvec(r,:);
                T{i} = zeros(int);
                Spp = zeros(int^2,size(Sp,2));
                for ii = -int+1:int-1
                    weight = sqrt(tlength-abs(ii));
                    T{i} = T{i} + tv(ii+int)/weight*diag(ones(int-abs(ii),1),ii);
                    inxxx = (diag(ones(int-abs(ii),1),ii));
                    inxxx = find(inxxx(:)>0);
                    for iii = 1:length(inxxx)
                    Spp(inxxx(iii),:) = Sp(ii+int,:)/weight;
                    end
                end
                T{i} = T{i}*sqrt(s(i,i));
                
            else %Standard.
                T{i} = reshape(uvec,int,int)*sqrt(s(i,i));
            end
            %T{i} = .5*(T{i} + T{i}');
            
            if Sflag
            S{i} = reshape(vvec,m*n,m*n)*sqrt(s(i,i));
            else
                S{i} = 0;
            end
            %S{i} = .5*(S{i}+S{i}');
            if Sflag
                kron_mat = kron_mat + kron(T{i}',S{i});
            else
                kron_mat = kron_mat;% + kron(T{i}',S{i});
            end
        end
    %end
%end

%Make sure result is symmetric.
%kron_mat = .5*(kron_mat + kron_mat');
if toep_flg
    
Sp = Spp;
end
if Sflag 
    %if ~toep_flg %ONLY WORKS IF NOT TOEPLITZ.
        SS = reshape(full(Sp)',m^2*int,int);
        funC = @(x) reshape(x.data,m,m)';
        Sp_mat = blockproc(SS,[m^2 1],funC);
        %Sp_mat = Sp;
    
    %else
        
    %end
    kron_mat = kron_mat+Sp_mat;
end



