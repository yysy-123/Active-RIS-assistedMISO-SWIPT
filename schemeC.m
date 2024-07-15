% No RIS
function [obj,flag]=schemeC(L,K,M,PBS,Rth,Eth, ...
    sigmaI,sigmaE,lamda,PSta,PBSCir,epsilon,Lmax,gd,hd)
        obj=zeros(1,2);
        f=sqrt(PBS/(L+K)/M)*ones(M,L);
        w=sqrt(PBS/(L+K)/M)*ones(M,K);
        flag=0;
        %% step1:optimize w and f
        iteD1=0;
        objD1=zeros(1,100);
        mu1=0;
        H=zeros(M,M,K);
        G=zeros(M,M,L);
        W=zeros(M,M,K);
        F=zeros(M,M,L);
        for k=1:K
            H(:,:,k)=hd(:,k)*hd(:,k)';
            W(:,:,k)=w(:,k)*w(:,k)';
        end
        for l=1:L
            G(:,:,l)=gd(:,l)*gd(:,l)';
            F(:,:,l)=f(:,l)*f(:,l)';
        end
        while 1
            ite1=0;
            obj1=zeros(1,100);
            while 1
                %optimize t
                t=zeros(1,K);
                for k=1:K
                    sum1=0;sum2=0;
                    for j=1:K
                        if j~=k
                            sum1=sum1+real(trace(H(:,:,k)*W(:,:,j)));
                        end
                    end
                    for i=1:L
                        sum2=sum2+real(trace(H(:,:,k)*F(:,:,i)));
                    end
                    t(k)=real(lamda*real(sum1+sum2)+1)^(-1);
                end
                %optimize W and F
                cvx_begin quiet
                variable W(M,M,K) complex semidefinite
                variable F(M,M,L) complex semidefinite
                expression phi(1,K)
                expression E1(1,L)
                expression Rth_right(1,K)
                expressions Ptotal_1 PowBS 
                for k=1:K
                    sumw=0;sum_w=0;sumf=0;
                    for j=1:K
                        sumw=sumw+real(trace(H(:,:,k)*W(:,:,j)));
                        if j~=k
                            sum_w=sum_w+real(trace(H(:,:,k)*W(:,:,j)));
                        end
                    end
                    for i=1:L
                        sumf=sumf+real(trace(H(:,:,k)*F(:,:,i)));
                    end
                    phi(k)=(log(real(lamda*(sumw+sumf)+1))-real(t(k))*real(lamda*(sum_w+sumf)+1)+log(real(t(k)))+1)/log(2);          
                    Rth_right(k)=lamda*(sum_w+sumf)+1;
                    PowBS=PowBS+real(trace(W(:,:,k)));
                end
                for l=1:L
                    PowBS=PowBS+real(trace(F(:,:,l)));
                    sumgw=0;sumgf=0;
                    for j=1:K
                       sumgw=sumgw+real(trace(G(:,:,l)*W(:,:,j))); 
                    end
                    for i=1:L
                       sumgf=sumgf+real(trace(G(:,:,l)*F(:,:,i)));
                    end
                    E1(l)=sumgw+sumgf+sigmaE;
                end
                Ptotal_1=PowBS+PSta+PBSCir;
                maximize sum(phi)-mu1*Ptotal_1
                subject to
                PowBS<=PBS;
                for k=1:K
                    lamda*real(trace(H(:,:,k)*W(:,:,k)))>=(2^Rth-1)*Rth_right(k);
                end
                for l=1:L
                    1e10*E1(l)>=1e10*Eth;
                end
                cvx_end
                if isnan(W(1,1,1)) || isnan(F(1,1,1))
                    flag=1;
                    break
                end
                %judge convergence
                ite1=ite1+1;
                obj1(ite1)=sum(phi)/Ptotal_1;
                if ite1>1 && (abs(obj1(ite1)-obj1(ite1-1))/abs(obj1(ite1-1))<=epsilon || ite1>=Lmax)
                    break
                end
            end%ending of lemma
            if flag==1
                break
            end
            mu1=sum(phi)/Ptotal_1;
            iteD1=iteD1+1;
            objD1(iteD1)=sum(phi)/Ptotal_1;
            if iteD1>1 && (abs(objD1(iteD1)-objD1(iteD1-1))/abs(objD1(iteD1-1))<=epsilon || iteD1>=Lmax)
                break
            end
        end%ending of sub1
        if flag==0
            for k=1:K
                [w1,eigw,w2]=svd(W(:,:,k));
                w(:,k)=sqrt(eigw(1,1))*w1(:,1);
            end
            for l=1:L
                [f1,eigf,f2]=svd(F(:,:,l));
                f(:,l)=sqrt(eigf(1,1))*f1(:,1);
            end
        end
        %% step3:judge convergence
        Rate=zeros(1,K);pow_BS=0;
        for k=1:K
            sum_int=0;             
            for j=1:K
                if j~=k
                    sum_int=sum_int+norm(hd(:,k)'*w(:,j))^2;
                end
            end
            for i=1:L
                sum_int=sum_int+norm(hd(:,k)'*f(:,i))^2;
            end
            sum_int=sum_int+sigmaI;
            Rate(k)=log2(1+norm(hd(:,k)'*w(:,k))^2/sum_int);
            pow_BS=pow_BS+norm(w(:,k))^2;
        end
        for l=1:L
            pow_BS=pow_BS+norm(f(:,l))^2;
        end
        Ptotal=pow_BS+PSta+PBSCir;
        obj(1)=sum(Rate)/Ptotal;
end