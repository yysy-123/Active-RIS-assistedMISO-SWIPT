function [R,gr,gd,hr,hd]=Channel_gene(M,N,K,L,C0,d0,IU,EU,BS,RIS,kr)        
        R=sqrt(C0*(d0/norm(RIS-BS)^2.2))*(sqrt(kr/(kr+1))*ones(N,M)+sqrt(1/(kr+1))*(sqrt(1/2)*randn(N,M)+sqrt(-1/2)*randn(N,M)));%BS-RIS
        gr=zeros(N,L);%RIS-EU
        gd=zeros(M,L);%BS-EU
        hr=zeros(N,K);%RIS-IU
        hd=zeros(M,K);%BS-IU
        for l=1:L
            gr(:,l)=sqrt(C0*(d0/norm(RIS-EU(l,:))^2.6))*(sqrt(kr/(kr+1))*ones(N,1)+sqrt(1/(kr+1))*(sqrt(1/2)*randn(N,1)+sqrt(-1/2)*randn(N,1)));
            gd(:,l)=sqrt(C0*(d0/norm(BS-EU(l,:))^3.5))*(sqrt(1/2)*randn(M,1)+sqrt(-1/2)*randn(M,1));
        end
        for k=1:K
            hr(:,k)=sqrt(C0*(d0/norm(RIS-IU(k,:))^2.6))*(sqrt(kr/(kr+1))*ones(N,1)+sqrt(1/(kr+1))*(sqrt(1/2)*randn(N,1)+sqrt(-1/2)*randn(N,1)));
            hd(:,k)=sqrt(C0*(d0/norm(BS-IU(k,:))^3.5))*(sqrt(1/2)*randn(M,1)+sqrt(-1/2)*randn(M,1));
        end
end