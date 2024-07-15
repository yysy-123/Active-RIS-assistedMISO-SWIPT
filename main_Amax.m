LL=600;
EE=5;
eeA=zeros(LL,EE);%active RIS
eeB=zeros(LL,EE);%passive RIS
eeC=zeros(LL,EE);%w/o RIS
eeD=zeros(LL,EE);%passive RIS with random phase
eeE=zeros(LL,EE);%AF relay
for e=1:EE
    BS=[0,0];
    RIS=[5,5];
    L=3;%number of EUs
    K=3;%number of IUs
    kr=10;%rician factor
    PBS=db2pow(30)*1e-3;%maximum transmit power at BS
    PI=db2pow(30)*1e-3;%maximum transmit power at RIS    
    N=24;%number of elements
    M=5;%number of antennas
    Rth=5;%QoS for IUs
    Eth=db2pow(-48)*1e-3;%QoS for EUs
    C0=db2pow(-30);%pass loss
    d0=1;%reference distance
    sigmaI=db2pow(-80)*1e-3;%AWGN at each IU
    sigmaE=db2pow(-80)*1e-3;%AWGN at each EU
    deltaI=db2pow(-70)*1e-3;%noise caused by RIS
    Amax=2*e;
    AFmax=10;
    lamda=1/sigmaI;
    PDyn=0.3e-3;PSta=100e-3;PBSCir=1.2;
    Lmax=20;
    epsilon=1e-2;
    parfor i=1:LL
        IU=IU_Location(50,10,5,K);
        EU=EU_Location(20,10,5,L);
        [R,gr,gd,hr,hd]=Channel_gene(M,N,K,L,C0,d0,IU,EU,BS,RIS,kr);    
        %% -----------------OPTIMIZATION----------------------------
        [objA,iteA,flagA]=schemeA(L,K,M,N,PBS,PI,Rth,Eth,sigmaI,sigmaE,deltaI,Amax,lamda,PDyn,PSta,PBSCir,epsilon,Lmax,R,gr,gd,hr,hd);
        eeA(i,e)=objA(iteA+flagA);
        disp(['A:',num2str(eeA(i,e)),',i=',num2str(i),',e=',num2str(e)]);
        [objB,iteB,flagB]=schemeB(L,K,M,N,PBS,Rth,Eth,sigmaI,sigmaE,lamda,PDyn,PSta,PBSCir,epsilon,Lmax,R,gr,gd,hr,hd);
        eeB(i,e)=objB(iteB+flagB);
        disp(['B:',num2str(eeB(i,e)),',i=',num2str(i),',e=',num2str(e)]);
        [objC,flagC]=schemeC(L,K,M,PBS,Rth,Eth,sigmaI,sigmaE,lamda,PSta,PBSCir,epsilon,Lmax,gd,hd);
        eeC(i,e)=objC(1+flagC);
        disp(['C:',num2str(eeC(i,e)),',i=',num2str(i),',e=',num2str(e)]);
        [objD,iteD,flagD]=schemeD(L,K,M,N,PBS,Rth,Eth,sigmaI,sigmaE,lamda,PDyn,PSta,PBSCir,epsilon,Lmax,R,gr,gd,hr,hd);
        eeD(i,e)=objD(iteD+flagD);
        disp(['D:',num2str(eeD(i,e)),',i=',num2str(i),',e=',num2str(e),',iter=',num2str(iteD)]);
        [objE,iteE,flagE]=schemeE(L,K,M,N,PBS,PI,Rth,Eth,sigmaI,sigmaE,deltaI,AFmax,lamda,PAF_relay,epsilon,Lmax,R,gr,gd,hr,hd);
        eeE(i,e)=objE(iteE+flagE);
        disp(['E:',num2str(eeE(i,e)),',i=',num2str(i),',e=',num2str(e),',iter=',num2str(iteE)]);
    end%ending of montecarlo
end%ending of sampling
dataA=zeros(1,EE);
dataB=zeros(1,EE);
dataC=zeros(1,EE);
dataD=zeros(1,EE);
dataE=zeros(1,EE);
for e=1:EE
    dataA(e)=sum(eeA(:,e))/sum(eeA(:,e)~=0);
    dataB(e)=sum(eeB(:,e))/sum(eeB(:,e)~=0);
    dataC(e)=sum(eeC(:,e))/sum(eeC(:,e)~=0);
    dataD(e)=sum(eeD(:,e))/sum(eeD(:,e)~=0);
    dataE(e)=sum(eeE(:,e))/sum(eeE(:,e)~=0);
end
x=2:2:(EE*2);
hold on;axis on;grid on;
plot(x,dataA,'-d','markersize',8,'linewidth',1.2);
plot(x,dataB,'-s','markersize',8,'linewidth',1.2);
plot(x,dataC,'-o','markersize',8,'linewidth',1.2);
plot(x,dataD,'-*','markersize',8,'linewidth',1.2);
plot(x,dataE,'-p','markersize',8,'linewidth',1.2);
legend('active RIS','passive RIS','w/o RIS','passive RIS with random phase','AF relay');
xlabel('Maximum amplitude at the active RIS');
ylabel('Energy efficiency (Mbits/J)');
set(gca,'xtick',x,'fontsize',12);
xlim([x(1) x(end)]);













