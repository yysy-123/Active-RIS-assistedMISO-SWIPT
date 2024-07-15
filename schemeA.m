% active RIS
function [obj,ite,flag]=schemeA(L,K,M,N,PBS,PI,Rth,Eth,sigmaI,sigmaE,deltaI,Amax,lamda,PDyn,PSta,PBSCir,epsilon,Lmax,R,gr,gd,hr,hd)
        ite=0;
        obj=zeros(1,100);
        f=sqrt(PBS/(L+K)/M)*ones(M,L);
        w=sqrt(PBS/(L+K)/M)*ones(M,K);
        v=Amax*ones(N,1);
        flag=0;
        while 1
            %% step1:optimize w and f
            iteD1=0;
            objD1=zeros(1,100);
            mu1=0;
            H=zeros(M,M,K);
            G=zeros(M,M,L);
            W=zeros(M,M,K);
            F=zeros(M,M,L);
            Y=R'*diag(v')*diag(v)*R;
            for k=1:K
                H(:,:,k)=(hd(:,k)'+hr(:,k)'*diag(v)*R)'*(hd(:,k)'+hr(:,k)'*diag(v)*R);
                W(:,:,k)=w(:,k)*w(:,k)';
            end
            for l=1:L
                G(:,:,l)=(gd(:,l)'+gr(:,l)'*diag(v)*R)'*(gd(:,l)'+gr(:,l)'*diag(v)*R);
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
                                sum1=sum1+trace(H(:,:,k)*W(:,:,j));
                            end
                        end
                        for i=1:L
                            sum2=sum2+trace(H(:,:,k)*F(:,:,i));
                        end
                        t(k)=real(lamda*(sum1+sum2+deltaI*norm(hr(:,k)'*diag(v))^2)+1)^(-1);
                    end
                    %optimize W and F
                    cvx_begin quiet
                    variable W(M,M,K) complex semidefinite
                    variable F(M,M,L) complex semidefinite
                    expression phi(1,K)
                    expression E1(1,L)
                    expression Rth_right(1,K)
                    expressions Ptotal_1 PowBS PowRIS
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
                        phi(k)=(log(real(lamda*(sumw+sumf+deltaI*norm(hr(:,k)'*diag(v))^2)+1))-real(t(k))*real(lamda*(sum_w+sumf+deltaI*norm(hr(:,k)'*diag(v))^2)+1)+log(real(t(k)))+1)/log(2);          
                        Rth_right(k)=lamda*real(sum_w+sumf+deltaI*norm(hr(:,k)'*diag(v))^2)+1;
                        PowBS=PowBS+real(trace(W(:,:,k)));
                        PowRIS=PowRIS+real(trace(Y*W(:,:,k)));
                    end
                    for l=1:L
                        PowBS=PowBS+real(trace(F(:,:,l)));
                        PowRIS=PowRIS+real(trace(Y*F(:,:,l)));
                        sumgw=0;sumgf=0;
                        for j=1:K
                           sumgw=sumgw+real(trace(G(:,:,l)*W(:,:,j))); 
                        end
                        for i=1:L
                           sumgf=sumgf+real(trace(G(:,:,l)*F(:,:,i)));
                        end
                        E1(l)=sumgw+sumgf+deltaI*norm(gr(:,l)'*diag(v))^2+sigmaE;
                    end
                    PowRIS=PowRIS+deltaI*norm(diag(v),'fro')^2;
                    Ptotal_1=PowBS+PowRIS+N*PDyn+PSta+PBSCir;
                    maximize sum(phi)-mu1*Ptotal_1
                    subject to
                    PowBS<=PBS;
                    PowRIS<=PI;
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
            if flag==1
                break
            end
            for k=1:K
                [w1,eigw,w2]=svd(W(:,:,k));
                w(:,k)=sqrt(eigw(1,1))*w1(:,1);
            end
            for l=1:L
                [f1,eigf,f2]=svd(F(:,:,l));
                f(:,l)=sqrt(eigf(1,1))*f1(:,1);
            end
            %% step2:optimize v
            T=zeros(N,N,K);
            T_=zeros(N,N,K);
            t=zeros(N,K);
            t_=zeros(N,K);
            d=zeros(1,K);
            d_=zeros(1,K);
            Hr=zeros(N,N,K);
            Gr=zeros(N,N,L);
            S=zeros(N,N,L);
            s=zeros(N,L);
            s_=zeros(1,L);
            U=zeros(N,N,L);
            u=zeros(N,L);
            u_=zeros(1,L);
            Rw=zeros(N,N);
            Rf=zeros(N,N);
            b_2=zeros(1,K);
            v_=v;
            powBS=0;
            for k=1:K
                powBS=powBS+norm(w(:,k))^2;
                suma=0;sumb=0;
                for j=1:K
                    suma=suma+w(:,j)*w(:,j)';
                    if j~=k
                        sumb=sumb+w(:,j)*w(:,j)';
                    end
                end
                for i=1:L
                    suma=suma+f(:,i)*f(:,i)';
                    sumb=sumb+f(:,i)*f(:,i)';
                end
                T(:,:,k)=diag(hr(:,k)')*R*suma*R'*diag(hr(:,k));
                T_(:,:,k)=diag(hr(:,k)')*R*sumb*R'*diag(hr(:,k));
                t(:,k)=diag(hr(:,k)')*R*suma*hd(:,k);
                t_(:,k)=diag(hr(:,k)')*R*sumb*hd(:,k);
                d(k)=hd(:,k)'*suma*hd(:,k);
                d_(k)=hd(:,k)'*sumb*hd(:,k);
                Hr(:,:,k)=diag(hr(:,k)')*diag(hr(:,k));
                Rw=Rw+diag(R*w(:,k))*diag(w(:,k)'*R');
                b_2(k)=lamda*real(v'*(T_(:,:,k)+deltaI*Hr(:,:,k))*v+v'*t_(:,k)+t_(:,k)'*v+d_(k))+1;
            end
            for l=1:L
                powBS=powBS+norm(f(:,l))^2;
                sumc=0;sumd=0;
                for j=1:K
                    sumc=sumc+w(:,j)*w(:,j)';
                end
                for i=1:L
                    sumd=sumd+f(:,i)*f(:,i)';
                end
                S(:,:,l)=diag(gr(:,l)')*R*sumc*R'*diag(gr(:,l));
                s(:,l)=diag(gr(:,l)')*R*sumc*gd(:,l);
                s_(l)=gd(:,l)'*sumc*gd(:,l);
                U(:,:,l)=diag(gr(:,l)')*R*sumd*R'*diag(gr(:,l));
                u(:,l)=diag(gr(:,l)')*R*sumd*gd(:,l);
                u_(l)=gd(:,l)'*sumd*gd(:,l);
                Gr(:,:,l)=diag(gr(:,l)')*diag(gr(:,l));
                Rf=Rf+diag(R*f(:,l))*diag(f(:,l)'*R');
            end
            iteD2=0;
            objD2=zeros(1,100);
            mu2=0;
            while 1
                ite2=0;
                obj2=zeros(1,100);
                while 1
                    %optimize v
                    cvx_begin quiet
                    variable v(N,1) complex
                    variable a1(1,K)
                    variable a2(1,K)
                    variable b1(1,K)
                    variable b2(1,K)
                    expression Ptotal_2 
                    expression E2(1,L)
                    Ptotal_2=real(quad_form(v,Rw+Rf+deltaI*eye(N))+powBS+N*PDyn+PSta+PBSCir);
                    for l=1:L
                        E2(l)=real(2*real(v_'*(S(:,:,l)+U(:,:,l)+deltaI*Gr(:,:,l))*v)-real(v_'*(S(:,:,l)+U(:,:,l)+deltaI*Gr(:,:,l))*v_)+v'*(s(:,l)+u(:,l))+(s(:,l)'+u(:,l)')*v+s_(l)+u_(l)+sigmaE);
                    end
                    maximize sum(a1)-sum(a2)-mu2*Ptotal_2
                    subject to
                    real(quad_form(v,Rw+Rf+deltaI*eye(N)))<=PI;
                    for n=1:N
                        norm(v(n))<=Amax;
                    end
                    for l=1:L
                        1e10*E2(l)>=1e10*Eth;
                    end
                    for k=1:K
                        a1(k)-a2(k)>=Rth;
                        log(b1(k))>=log(2)*a1(k);
                        lamda*(2*real(v_'*(T(:,:,k)+deltaI*Hr(:,:,k))*v)-real(v_'*(T(:,:,k)+deltaI*Hr(:,:,k))*v_)+real(v'*t(:,k)+t(:,k)'*v+d(k)))+1>=b1(k);
                        lamda*(quad_form(v,T_(:,:,k)+deltaI*Hr(:,:,k))+real(v'*t_(:,k)+t_(:,k)'*v+d_(k)))+1<=b2(k);
                        log(b_2(k))+(b2(k)-b_2(k))/b_2(k)<=log(2)*a2(k);
                    end
                    cvx_end
                    if isnan(v(1))
                        flag=1;
                        break
                    end
                    %update
                    v_=v;
                    b_2=b2;
                    %judge convergence
                    ite2=ite2+1;
                    obj2(ite2)=(sum(a1)-sum(a2))/Ptotal_2;
                    if ite2>1 && (abs(obj2(ite2)-obj2(ite2-1))/abs(obj2(ite2-1))<=epsilon || ite2>=Lmax)
                        break
                    end
                end%ending of sca
                if flag==1
                    break
                end
                mu2=(sum(a1)-sum(a2))/Ptotal_2;
                iteD2=iteD2+1;
                objD2(iteD2)=(sum(a1)-sum(a2))/Ptotal_2;
                if iteD2>1 && (abs(objD2(iteD2)-objD2(iteD2-1))/abs(objD2(iteD2-1))<=epsilon || iteD2>=Lmax)
                    break
                end
            end%ending of sub2
            if flag==1
                break
            end
            v=conj(v);
            %% step3:judge convergence
            Rate=zeros(1,K);pow_BS=0;pow_RIS=0;   
            for k=1:K
                sum_int=0;             
                for j=1:K
                    if j~=k
                        sum_int=sum_int+norm((hd(:,k)'+hr(:,k)'*diag(v)*R)*w(:,j))^2;
                    end
                end
                for i=1:L
                    sum_int=sum_int+norm((hd(:,k)'+hr(:,k)'*diag(v)*R)*f(:,i))^2;
                end
                sum_int=sum_int+deltaI*norm(hr(:,k)'*diag(v))^2+sigmaI;
                Rate(k)=log2(1+norm((hd(:,k)'+hr(:,k)'*diag(v)*R)*w(:,k))^2/sum_int);
                pow_BS=pow_BS+norm(w(:,k))^2;
                pow_RIS=pow_RIS+norm(diag(v)*R*w(:,k))^2;
            end
            for l=1:L
                pow_BS=pow_BS+norm(f(:,l))^2;
                pow_RIS=pow_RIS+norm(diag(v)*R*f(:,l))^2;
            end
            Ptotal=pow_BS+pow_RIS+deltaI*norm(diag(v),'fro')^2+N*PDyn+PSta+PBSCir;
            ite=ite+1;
            obj(ite)=sum(Rate)/Ptotal;
            if ite>1 && (abs(obj(ite)-obj(ite-1))/abs(obj(ite-1))<=epsilon || ite>=Lmax)
                break
            end
        end%ending of AO
end