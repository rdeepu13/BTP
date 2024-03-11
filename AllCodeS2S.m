%13/8/2021

clc;
clear;
close all;
%rng('default')



 format long;


    T=1;
M=1;
% sigThdeg=1.5;
timespan=(30*M);

if rem(timespan,T)>0
    
    remainder=rem(timespan,T);
    timespan=(timespan-remainder);
end

% knot2=0.000514444;%1knot=0.000514444km/sec-------------------
knot2=0.0308667;%1knot=0.0308667km/min

%
%---------------------S2S In Min--------------------------------
xt1=4.9286;
yt1=.8420;
osm=0.15433;
tsm=0.123466;
tarc=-140*pi/180;
xori = S2SInMin(timespan);
q=1.944*10^-6;%-------------------------in km^2/min^3
sigThdeg=1.5;
sigR=2;%in km
%-------------------------------------------------------------------
%}
%{
%---------------------S2S In Min--------------------------------
xt1=7;
yt1=7;
osm=0.15433;
tsm=0.463;%15 knot
tarc=-135.4*pi/180;
xori = S2SInMinHighlyNonlinear(timespan);
q=1.944*10^-6;%-------------------------in km^2/min^3
sigThdeg=2;
sigR=4;%in km
%-------------------------------------------------------------------
%}

HRH=zeros(4,4,timespan);
SumP0=zeros(4,4);
MCrun=10;

err=zeros(MCrun,timespan/T);
errv=zeros(MCrun,timespan/T);
errekf=zeros(MCrun,timespan/T);
errekfv=zeros(MCrun,timespan/T);

ctcckf=0;
ctccqkf=0;
ctcukf=0;
ctcghf=0;

ctcekf=0;

cerrckf=zeros(MCrun,timespan/T);
cerrckfv=zeros(MCrun,timespan/T);
cerrcqkf=zeros(MCrun,timespan/T);
cerrcqkfv=zeros(MCrun,timespan/T);
cerrukf=zeros(MCrun,timespan/T);
cerrukfv=zeros(MCrun,timespan/T);
cerrghf=zeros(MCrun,timespan/T);
cerrghfv=zeros(MCrun,timespan/T);

cerrekf=zeros(MCrun,timespan/T);
cerrekfv=zeros(MCrun,timespan/T);
err2=zeros(MCrun,timespan/T);
err2v=zeros(MCrun,timespan/T);
ctc2=0;
cerr2=zeros(MCrun,timespan/T);
cerr2v=zeros(MCrun,timespan/T);

for run=1:MCrun
    
    F=[1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
%     q=1.944*10^-6/(60^3);%-------------------------in km^2/sec^3
    Q=[T^3/3 0 T^2/2 0; 0 T^3/3 0 T^2/2; T^2/2 0 T 0; 0 T^2/2 0 T]*q;

     xo=zeros(4,timespan/T);
     k=1;
     for i=1:T:timespan
         xo(:,k)=xori(:,i);
         k=k+1;
     end
    
    xk=zeros(4,timespan/T);
    xk(:,1)=[xt1-xo(1,1)
        yt1-xo(2,1)
        tsm*sin(tarc)-xo(3,1)
        tsm*cos(tarc)-xo(4,1)];

    U=zeros(4,timespan/T);
    %v=zeros(4,timespan+1);
    z=zeros(1,timespan/T);
    zdegree=zeros(1,timespan/T);
    bz=zeros(2,timespan/T+1);
    
    sigTh=(pi*sigThdeg)/180;%Standard Deviation of theta
    R=(sigTh)^2;
    
    z(1,1)=atan2(xk(1,1),xk(2,1))+0+sqrt(R)*randn(1);
    zdegree(1,1)=atan2(xk(1,1),xk(2,1))*180/pi;
    bz(:,1)=[xk(1,1) xk(2,1)]';
    
    %wk=zeros(1,timespan-1);
    k=1;
    for t=1:timespan/T-1
        
        U(:,t)=[xo(1,t+1)-xo(1,t)-(xo(3,t)*T)
            xo(2,t+1)-xo(2,t)-(xo(4,t)*T)
            xo(3,t+1)-xo(3,t)
            xo(4,t+1)-xo(4,t)];
        v=(mvnrnd(zeros(4,1),Q))';
        xk(:,t+1)=F*xk(:,t)-U(:,t)+v;
        wk=sqrt(R)*randn(1);
        zdegree(1,t+1)=atan2(xk(1,t+1),xk(2,t+1))*180/pi;
        z(1,t+1)=atan2(xk(1,t+1),xk(2,t+1))+wk;
        bz(:,t+1)=[xk(1,t+1) xk(2,t+1)]';
        k=k+1;
        
    end
    xtrack=zeros(4,timespan/T);
    
    for t=1:timespan/T
        xtrack(:,t)=xo(:,t)+xk(:,t);
    end

    %-------------------------Initial Covariance and Relative Position-------------------------------

    for t=1:timespan/T
        x=xk(1,t);
        y=xk(2,t);
        H=[y/(x^2+y^2) -x/(x^2+y^2) 0 0];
        HRH(:,:,t)=HRH(:,:,t)+(H'*R^-1*H);
    end

    
    sigTh = (pi*sigThdeg)/180;%in radian
%     sigR=2;%in km
    sigS=2*knot2;  %2knot
    sigC=pi/sqrt(12);
    

     
        r=(sqrt((xk(1,1))^2+(xk(2,1))^2))+((sigR)*randn(1));
        th0=z(1,1);
        %c=pi+th0;%traj 4
        c=pi+th0;%traj 3
        s=(tsm)+(sigS*randn(1));

     Pxx=(r*sigTh*cos(th0))^2+(sigR*sin(th0))^2;
     Pyy=(r*sigTh*sin(th0))^2+(sigR*cos(th0))^2;
     Pxy=(sigR^2-(r*sigTh)^2)*sin(th0)*cos(th0);
     Pyx=Pxy;
     Pxxdot=(s*sigC*cos(c))^2+(sigS*sin(c))^2;
     Pyydot=(s*sigC*sin(c))^2+(sigS*cos(c))^2;
     Pxydot=(sigS^2-(s*sigC)^2)*sin(c)*cos(c);
     Pyxdot=Pxydot;

    P0=[Pxx Pxy 0 0
        Pyx Pyy 0 0 
        0 0 Pxxdot Pxydot
        0 0 Pyxdot Pyydot];
 
    SumP0=SumP0+P0;
    x0=[r*sin(th0)
        r*cos(th0)
        (s*sin(c))-xo(3,1)
        (s*cos(c))-xo(4,1)];
    nx=4; 
%     P0=[1.24394312931410,1.79262898307086,0,0;1.79262898307086,2.83401583431552,0,0;0,0,0.000122851897844758,-7.92183133556663e-05;0,0,-7.92183133556663e-05,5.25847767657967e-05];
%     x0=[5.81506617751481;8.94030680268175;-0.0111063134815636;-0.00350447435491635];
    %
    %------------------------------------------------------
    for method=1:4
        if method==1%ckf
                n1 = 1;
                [ep, W] = cqkf_p(nx,n1);
                [~, col] = size(ep);
        end
        if method ==2%cqkf
                n1 = 2;
                [ep, W] = cqkf_p(nx,n1);
                [~, col] = size(ep);
        end
        if method==3%ukf
            kappa = -1; %nx+kappa = 3
            ep = [zeros(nx, 1) sqrt(nx+kappa)*eye(nx) -sqrt(nx+kappa)*eye(nx)];
            W = [kappa/(nx+kappa) 1/(2*(nx+kappa))*ones(1, 2*nx)];
            [~, col] = size(ep);
        end
        if method==4%ghf
            alpha = 3; %number of ponts alpha^nx
            [ep, W] = ghf_p(alpha);
            [~, col] = size(ep);
        end
    xh=zeros(4,timespan/T);
    zh=zeros(1,timespan/T);
    xn=zeros(1,timespan/T);
    yn=zeros(1,timespan/T);
    Pc=P0;
%      display(method)
    xh(:,1)=x0;
%     display(xk(:,1))
    
    xn(1)=xh(1,1)+xo(1,1);
    yn(1)=xh(2,1)+xo(2,1);
    
    
    
    err(run,1)=((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2);
    errv(run,1)=((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2);
    abserr(run,1)=sqrt(((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2));
    trueval(run,1)=sqrt(xk(1,1)^2+xk(2,1)^2);
    abserrv(run,1)=sqrt(((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2));
    truevalv(run,1)=sqrt(xk(3,1)^2+xk(4,1)^2);
    abserrc(run,1)=abs(atan2(xk(1,1),xk(2,1))-atan2(xh(1,1),xh(2,1)));
    truevalc(run,1)=abs(atan2(xk(1,1),xk(2,1)));
%     display(err(1,1))
%      display(errv(1,1))
%       display(abserr(1,1))
%       display(trueval(1,1))
%       display(abserrv(1,1))
%       display(truevalv(1,1))
%      display(abserrc(1,1))
%      display(truevalc(1,1))
    
    for t=2:timespan/T

        %---Time update ----
        % System is linear, so we can use here Kalman filter
        Xpri = F*xh(:,t-1) - U(:,t-1);
        Ppri = F*Pc*F' + Q;
                
        S = chol(Ppri,'lower');
        kiipred = zeros(nx, col); zpred = zeros(1, col);
        for i = 1:col
            kiipred(:, i) = S*ep(:, i) + Xpri;
            zpred(i) = atan2(kiipred(1, i), kiipred(2, i));
        end
       
        zhat = 0;
        for i = 1:col
            zhat = zhat + W(i)*zpred(i);
        end
        
        Pzz = 0;
        Pxz = zeros(nx, 1);
        for i = 1:col
            Pzz = Pzz + W(i)*(zpred(i) - zhat)^2;
            Pxz = Pxz + W(i)*(kiipred(:, i) - Xpri)*(zpred(i) - zhat)';
        end
        Pzz = Pzz + R;
        Wk = Pxz*((Pzz)^-1);
        xh(:,t) = Xpri + (Wk*(z(:,t) - zhat));
        Pc = Ppri - (Wk*Pzz*Wk');
        
        xn(t)=xh(1,t)+xo(1,t);
        yn(t)=xh(2,t)+xo(2,t);
        
        err(run,t)=((xh(1,t)-xk(1,t))^2+(xh(2,t)-xk(2,t))^2);
        errv(run,t)=((xh(3,t)-xk(3,t))^2+(xh(4,t)-xk(4,t))^2);
        abserr(run,t)=sqrt(((xh(1,t)-xk(1,t))^2+(xh(2,t)-xk(2,t))^2));
        trueval(run,t)=sqrt(xk(1,t)^2+xk(2,t)^2);
        abserrv(run,t)=sqrt(((xh(3,t)-xk(3,t))^2+(xh(4,t)-xk(4,t))^2));
        truevalv(run,t)=sqrt(xk(3,t)^2+xk(4,t)^2);
        abserrc(run,t)=abs(atan2(xk(1,t),xk(2,t))-atan2(xh(1,t),xh(2,t)));
        truevalc(run,t)=abs(atan2(xk(1,t),xk(2,t)));
        
    end
  
    if (err(run,(timespan/T)))<=1 && method==1
        ctcckf=ctcckf+1;
        cerrckf(ctcckf,:)=err(run,:);
        cerrckfv(ctcckf,:)=errv(run,:); 
        cabserrckf(ctcckf,:)=abserr(run,:);
        ctruevalckf(ctcckf,:)=trueval(run,:);
        cabserrckfv(ctcckf,:)=abserrv(run,:);
        ctruevalckfv(ctcckf,:)=truevalv(run,:);
        cabserrckfc(ctcckf,:)=abserrc(run,:);
        ctruevalckfc(ctcckf,:)=truevalc(run,:);
       
    end
    if (err(run,(timespan/T)))<=1 && method==2
        ctccqkf=ctccqkf+1;
        cerrcqkf(ctccqkf,:)=err(run,:);
        cerrcqkfv(ctccqkf,:)=errv(run,:);  
        cabserrcqkf(ctccqkf,:)=abserr(run,:);
        ctruevalcqkf(ctccqkf,:)=trueval(run,:);
        cabserrcqkfv(ctccqkf,:)=abserrv(run,:);
        ctruevalcqkfv(ctccqkf,:)=truevalv(run,:);
        cabserrcqkfc(ctccqkf,:)=abserrc(run,:);
        ctruevalcqkfc(ctccqkf,:)=truevalc(run,:);
        
    end
    if (err(run,(timespan/T)))<=1 && method==3
        ctcukf=ctcukf+1;
        cerrukf(ctcukf,:)=err(run,:);
        cerrukfv(ctcukf,:)=errv(run,:);  
        cabserrukf(ctcukf,:)=abserr(run,:);
        ctruevalukf(ctcukf,:)=trueval(run,:);
        cabserrukfv(ctcukf,:)=abserrv(run,:);
        ctruevalukfv(ctcukf,:)=truevalv(run,:);
        cabserrukfc(ctcukf,:)=abserrc(run,:);
        ctruevalukfc(ctcukf,:)=truevalc(run,:);
    end
    if (err(run,(timespan/T)))<=1 && method==4
        ctcghf=ctcghf+1;
        cerrghf(ctcghf,:)=err(run,:);
        cerrghfv(ctcghf,:)=errv(run,:);  
        cabserrghf(ctcghf,:)=abserr(run,:);
        ctruevalghf(ctcghf,:)=trueval(run,:);
        cabserrghfv(ctcghf,:)=abserrv(run,:);
        ctruevalghfv(ctcghf,:)=truevalv(run,:);
        cabserrghfc(ctcghf,:)=abserrc(run,:);
        ctruevalghfc(ctcghf,:)=truevalc(run,:);
    end
    end
    %------------------------------EKF--------------------------------
    xh=zeros(4,timespan/T);
    xekf=zeros(1,timespan/T);
    yekf=zeros(1,timespan/T);
%      display('EKF')
    xh(:,1)=x0;
%     display(xk(:,1))
    P=P0;
    xekf(1)=xh(1,1)+xo(1,1);
    yekf(1)=xh(2,1)+xo(2,1);
    errekf(run,1)=((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2);
    errekfv(run,1)=((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2);
    abserr(run,1)=sqrt(((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2));
    trueval(run,1)=sqrt(xk(1,1)^2+xk(2,1)^2);
    abserrv(run,1)=sqrt(((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2));
    truevalv(run,1)=sqrt(xk(3,1)^2+xk(4,1)^2);
    abserrc(run,1)=abs(atan2(xk(1,1),xk(2,1))-atan2(xh(1,1),xh(2,1)));
    truevalc(run,1)=abs(atan2(xk(1,1),xk(2,1)));
%     
%       display(errekf(1,1))
%      display(errekfv(1,1))
%       display(abserr(1,1))
%       display(trueval(1,1))
%       display(abserrv(1,1))
%       display(truevalv(1,1))
%      display(abserrc(1,1))
%      display(truevalc(1,1))
    for t=2:timespan/T
        %time update
        xpri = F*xh(:,t-1) - U(:,t-1);
        ppri = F*P*F' + Q;
        
        zhat = atan2(xpri(1), xpri(2));
        
        Hjac = [xpri(2)/(xpri(1)^2+xpri(2)^2), -xpri(1)/(xpri(1)^2+xpri(2)^2), 0, 0]; 
        
        Pzz = Hjac*ppri*Hjac' + R; 
        
        Pxz = ppri*Hjac'; 
        
        K = Pxz/Pzz;
        
        xh(:, t) = xpri + K*(z(t) - zhat);
        P = ppri - K*Pzz*K';    
        
        xekf(t)=xh(1,t)+xo(1,t);
        yekf(t)=xh(2,t)+xo(2,t);
        
        errekf(run,t)=((xh(1,t)-xk(1,t))^2+(xh(2,t)-xk(2,t))^2);
        errekfv(run,t)=((xh(3,t)-xk(3,t))^2+(xh(4,t)-xk(4,t))^2);
        abserr(run,t)=sqrt(((xh(1,t)-xk(1,t))^2+(xh(2,t)-xk(2,t))^2));
        trueval(run,t)=sqrt(xk(1,t)^2+xk(2,t)^2);
        abserrv(run,t)=sqrt(((xh(3,t)-xk(3,t))^2+(xh(4,t)-xk(4,t))^2));
        truevalv(run,t)=sqrt(xk(3,t)^2+xk(4,t)^2);
        abserrc(run,t)=abs(atan2(xk(1,t),xk(2,t))-atan2(xh(1,t),xh(2,t)));
        truevalc(run,t)=abs(atan2(xk(1,t),xk(2,t)));
        
    end
    if (errekf(run,(timespan/T)))<=1
        ctcekf=ctcekf+1;
        cerrekf(ctcekf,:)=errekf(run,:);
        cerrekfv(ctcekf,:)=errekfv(run,:);
        cabserrekf(ctcekf,:)=abserr(run,:);
        ctruevalekf(ctcekf,:)=trueval(run,:);
        cabserrekfv(ctcekf,:)=abserrv(run,:);
        ctruevalekfv(ctcekf,:)=truevalv(run,:);
        cabserrekfc(ctcekf,:)=abserrc(run,:);
        ctruevalekfc(ctcekf,:)=truevalc(run,:);
    end
        %-----------------------------------SRF--------------------------------------
    
%     h=[1 0 0 0
%         0 1 0 0];
%     xh=zeros(4,timespan/T);
%     xn=zeros(1,timespan/T);
%     yn=zeros(1,timespan/T);
%     xh(:,1)=x0;
%     xn(1)=xh(1,1)+xo(1,1);
%     yn(1)=xh(2,1)+xo(2,1);
%     xh(:,1)=x0;
% %     P0=[1.3232    1.8286         0         0
% %         1.8286    2.7508         0         0
% %              0         0    0.0003   -0.0002
% %              0         0   -0.0002    0.0001];
%     P=P0;
%     
%     err2(run,1)=((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2);
%     err2v(run,1)=((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2);
%     abserr(run,t)=sqrt(((xh(1,1)-xk(1,1))^2+(xh(2,1)-xk(2,1))^2));
%     trueval(run,t)=sqrt(xk(1,1)^2+xk(2,1)^2);
%     abserrv(run,t)=sqrt(((xh(3,1)-xk(3,1))^2+(xh(4,1)-xk(4,1))^2));
%     truevalv(run,t)=sqrt(xk(3,1)^2+xk(4,1)^2);
%     abserrc(run,t)=abs(atan2(xk(1,1),xk(2,1))-atan2(xh(1,1),xh(2,1)));
%     truevalc(run,t)=abs(atan2(xk(1,1),xk(2,1)));
%     ERR(:,1)=xk(:,1)-xh(:,1);
% %     display('SRF')
%     %      display(err2(1,1))
% %      display(err2v(1,1))
% %       display(abserr(1,1))
% %       display(trueval(1,1))
% %       display(abserrv(1,1))
% %       display(truevalv(1,1))
% %      display(abserrc(1,1))
% %      display(truevalc(1,1))
%     
%     arrdelta=zeros(1,timespan/T-1);
%     zt=zeros(1,timespan/T-1);
%     rho=zeros(1,timespan/T-1);
%     rho1=zeros(1,timespan/T-1);
%     checking=zeros(1,timespan/T-1);
%     gama=zeros(1,timespan/T-1);
%     delta=zeros(1,timespan/T-1);
%     Qmmag=zeros(1,timespan/T-1);
%     for t=1:timespan/T-1
%         xhat=(F*xh(:,t))-U(:,t);
%         Pnew=(F*P*F')+Q;
%     
%         %Qm=(sigTh^2)*(xhat(1)^2+xhat(2)^2+Pnew(1,1)+Pnew(2,2));
%          Qm=(sigTh^2)*(xhat(1)^2+xhat(2)^2+Pnew(1,1)+Pnew(2,2))*eye(2,2);
%         %Qm=sigTh^2*(sqrt(xhat(1,1)^2+xhat(2,1)^2))^2*eye(2,2);
% 
%         V=(h*Pnew*(h'))+Qm;
%         Kt=(Pnew*(h')*(V)^-1);
%         
%         %wkk=(mvnrnd(zeros(2,1),Qm))';
%         %bk=bz(:,t+1)+wkk;
% %         bk=[xhat(1,1) xhat(2,1)]';
% %           bk=[xh(1,t) xh(2,t)]';
%           
%         
%           bk=[sin(z(t+1))
%               cos(z(t+1))];
%         zt(1,t)=(((bk')*(V)^-1*bk)^(-0.5))*bk'*(V)^-1*(h*xhat);
%         
%         checking=(bk'*(V)^-1*(h*xhat));
%         %display(bk'*(V)^-1*(h*xhat))
%         
%         Fn=0.5*erfc(-zt(1,t)/sqrt(2));
%         %display(Fn)
%         
%         rho1(1,t) = zt(1,t)+1/(zt(1,t)+exp(-zt(1,t)^2/2)/(sqrt(2*pi)*Fn));
%         
%         if zt(1,t)>=37
%             rho(1,t)=zt(1,t);
%         else
%             
%             rho(1,t) = (zt(1,t)+ (sqrt(2*pi)*((zt(1,t)^2) + 1)*exp(0.5*(zt(1,t)^2))*Fn))/(1+ (sqrt(2*pi)*zt(1,t)*exp(0.5*(zt(1,t)^2))* Fn));
%         end
%         gama(1,t)=((bk'*(V)^-1*bk)^(-0.5))*rho(1,t);
%         delta(1,t)=((bk'*(V)^-1*bk)^(-1))*(2+(zt(1,t)*rho(1,t))-(rho(1,t)^2));%----------------
%         %delta(1,t)=((bk'*(V)^-1*bk)^(-1))*(2+(zt(1,t)*rho(1,t))-(rho(1,t)^2))+10;
%         
%         xh(:,t+1)=((eye(4)-(Kt*h))*xhat)+(gama(1,t)*Kt*bk);
%         P=((eye(4)-(Kt*h))*Pnew)+(delta(1,t)*Kt*bk*(bk')*(Kt'));
%         
%         err2(run,t+1)=((xh(1,t+1)-xk(1,t+1))^2+(xh(2,t+1)-xk(2,t+1))^2);
%         err2v(run,t+1)=((xh(3,t+1)-xk(3,t+1))^2+(xh(4,t+1)-xk(4,t+1))^2);
%         abserr(run,t+1)=sqrt(((xh(1,t+1)-xk(1,t+1))^2+(xh(2,t+1)-xk(2,t+1))^2));
%         trueval(run,t+1)=sqrt(xk(1,t+1)^2+xk(2,t+1)^2);
%         abserrv(run,t+1)=sqrt(((xh(3,t+1)-xk(3,t+1))^2+(xh(4,t+1)-xk(4,t+1))^2));
%         truevalv(run,t+1)=sqrt(xk(3,t+1)^2+xk(4,t+1)^2);
%         abserrc(run,t+1)=abs(atan2(xk(1,t+1),xk(2,t+1))-atan2(xh(1,t+1),xh(2,t+1)));
%         truevalc(run,t+1)=abs(atan2(xk(1,t+1),xk(2,t+1)));
%         
%         xn(t+1)=xh(1,t+1)+xo(1,t+1);
%         yn(t+1)=xh(2,t+1)+xo(2,t+1);
% 
%         arrdelta(t)=delta(1,t);
%     end
%     
%     if (err2(run,(timespan/T)))<=1
%         
%         ctc2=ctc2+1;
%         cerr2(ctc2,:)=err2(run,:);
%         cerr2v(ctc2,:)=err2v(run,:);
%         cabserrsrf(ctc2,:)=abserr(run,:);
%         ctruevalsrf(ctc2,:)=trueval(run,:);
%         cabserrsrfv(ctc2,:)=abserrv(run,:);
%         ctruevalsrfv(ctc2,:)=truevalv(run,:);
%         cabserrsrfc(ctc2,:)=abserrc(run,:);
%         ctruevalsrfc(ctc2,:)=truevalc(run,:);
%         %{
%         figure(run+10)
%         plot(xtrack(1,:),xtrack(2,:),'c','Linewidth',1.5);
%         hold on;
%         plot(xn,yn,'g','Linewidth',1.5);
%         hold on;
%         plot(xo(1,:),xo(2,:),'r','Linewidth',1.5);
%         legend('Target True Trajectory','Target Estimated Trajectory','Observer Trajectory');
%         xlabel('x axis in km');
%         ylabel('y axis in km');
%         title('Trajectory Plot SRF');
%     else
%         figure(run+10)
%         plot(xo(1,:),xo(2,:),'r','Linewidth',1.5);
%         hold on;
%         plot(xtrack(1,:),xtrack(2,:),'c');
%         hold on;
%         plot(xn,yn,'r');
%         hold on;
%         title('Trajectory Plot SRF');
%         %}
%     end
fprintf('MCrun:%d completed\n',run)
%     disp(run)
    
end

%------------------------------------------------------
sumerrckf=zeros(1,timespan/T);
sumerrckfv=zeros(1,timespan/T);
sumabserrckf=zeros(1,timespan/T);
sumtruevalckf=zeros(1,timespan/T);
sumabserrckfv=zeros(1,timespan/T);
sumtruevalckfv=zeros(1,timespan/T);
sumabserrckfc=zeros(1,timespan/T);
sumtruevalckfc=zeros(1,timespan/T);
for k=1:ctcckf
    sumerrckf(1,:)=sumerrckf(1,:)+cerrckf(k,:);
    sumerrckfv(1,:)=sumerrckfv(1,:)+cerrckfv(k,:);
    sumabserrckf(1,:)=sumabserrckf(1,:)+cabserrckf(k,:);
    sumtruevalckf(1,:)=sumtruevalckf(1,:)+ctruevalckf(k,:);
    sumabserrckfv(1,:)=sumabserrckfv(1,:)+cabserrckfv(k,:);
    sumtruevalckfv(1,:)=sumtruevalckfv(1,:)+ctruevalckfv(k,:);
    sumabserrckfc(1,:)=sumabserrckfc(1,:)+cabserrckfc(k,:);
    sumtruevalckfc(1,:)=sumtruevalckfc(1,:)+ctruevalckfc(k,:);
end
rmserrckf=zeros(1,timespan/T);
rmserrckfv=zeros(1,timespan/T);
relerrckf=zeros(1,timespan/T);
relerrckfv=zeros(1,timespan/T);
relerrckfc=zeros(1,timespan/T);
for t=1:timespan/T
    rmserrckf(1,t)=sqrt(sumerrckf(1,t)/ctcckf);
    rmserrckfv(1,t)=sqrt(sumerrckfv(1,t)/ctcckf);
    relerrckf(1,t)=sumabserrckf(1,t)/sumtruevalckf(1,t);
    relerrckfv(1,t)=sumabserrckfv(1,t)/sumtruevalckfv(1,t);
    relerrckfc(1,t)=sumabserrckfc(1,t)/sumtruevalckfc(1,t);
end

CKFTrackLoss=((MCrun-ctcckf)/MCrun)*100;
display(CKFTrackLoss)
%--------------------------------------------------------
sumerrcqkf=zeros(1,timespan/T);
sumerrcqkfv=zeros(1,timespan/T);
sumabserrcqkf=zeros(1,timespan/T);
sumtruevalcqkf=zeros(1,timespan/T);
sumabserrcqkfv=zeros(1,timespan/T);
sumtruevalcqkfv=zeros(1,timespan/T);
sumabserrcqkfc=zeros(1,timespan/T);
sumtruevalcqkfc=zeros(1,timespan/T);
for k=1:ctccqkf
    sumerrcqkf(1,:)=sumerrcqkf(1,:)+cerrcqkf(k,:);
    sumerrcqkfv(1,:)=sumerrcqkfv(1,:)+cerrcqkfv(k,:);
    sumabserrcqkf(1,:)=sumabserrcqkf(1,:)+cabserrcqkf(k,:);
    sumtruevalcqkf(1,:)=sumtruevalcqkf(1,:)+ctruevalcqkf(k,:);
    sumabserrcqkfv(1,:)=sumabserrcqkfv(1,:)+cabserrcqkfv(k,:);
    sumtruevalcqkfv(1,:)=sumtruevalcqkfv(1,:)+ctruevalcqkfv(k,:);
    sumabserrcqkfc(1,:)=sumabserrcqkfc(1,:)+cabserrcqkfc(k,:);
    sumtruevalcqkfc(1,:)=sumtruevalcqkfc(1,:)+ctruevalcqkfc(k,:);
end
rmserrcqkf=zeros(1,timespan/T);
rmserrcqkfv=zeros(1,timespan/T);
relerrcqkf=zeros(1,timespan/T);
relerrcqkfv=zeros(1,timespan/T);
relerrcqkfc=zeros(1,timespan/T);
for t=1:timespan/T
    rmserrcqkf(1,t)=sqrt(sumerrcqkf(1,t)/ctccqkf);
    rmserrcqkfv(1,t)=sqrt(sumerrcqkfv(1,t)/ctccqkf);
    relerrcqkf(1,t)=sumabserrcqkf(1,t)/sumtruevalcqkf(1,t);
    relerrcqkfv(1,t)=sumabserrcqkfv(1,t)/sumtruevalcqkfv(1,t);
    relerrcqkfc(1,t)=sumabserrcqkfc(1,t)/sumtruevalcqkfc(1,t);
end
CQKFTrackLoss=((MCrun-ctccqkf)/MCrun)*100;
display(CQKFTrackLoss)
%------------------------------------------------------
sumerrukf=zeros(1,timespan/T);
sumerrukfv=zeros(1,timespan/T);
sumabserrukf=zeros(1,timespan/T);
sumtruevalukf=zeros(1,timespan/T);
sumabserrukfv=zeros(1,timespan/T);
sumtruevalukfv=zeros(1,timespan/T);
sumabserrukfc=zeros(1,timespan/T);
sumtruevalukfc=zeros(1,timespan/T);
for k=1:ctcukf
    sumerrukf(1,:)=sumerrukf(1,:)+cerrukf(k,:);
    sumerrukfv(1,:)=sumerrukfv(1,:)+cerrukfv(k,:);
    sumabserrukf(1,:)=sumabserrukf(1,:)+cabserrukf(k,:);
    sumtruevalukf(1,:)=sumtruevalukf(1,:)+ctruevalukf(k,:);
    sumabserrukfv(1,:)=sumabserrukfv(1,:)+cabserrukfv(k,:);
    sumtruevalukfv(1,:)=sumtruevalukfv(1,:)+ctruevalukfv(k,:);
    sumabserrukfc(1,:)=sumabserrukfc(1,:)+cabserrukfc(k,:);
    sumtruevalukfc(1,:)=sumtruevalukfc(1,:)+ctruevalukfc(k,:);
end
rmserrukf=zeros(1,timespan/T);
rmserrukfv=zeros(1,timespan/T);
relerrukf=zeros(1,timespan/T);
relerrukfv=zeros(1,timespan/T);
relerrukfc=zeros(1,timespan/T);
for t=1:timespan/T
    rmserrukf(1,t)=sqrt(sumerrukf(1,t)/ctcukf);
    rmserrukfv(1,t)=sqrt(sumerrukfv(1,t)/ctcukf);
    relerrukf(1,t)=sumabserrukf(1,t)/sumtruevalukf(1,t);
    relerrukfv(1,t)=sumabserrukfv(1,t)/sumtruevalukfv(1,t);
    relerrukfc(1,t)=sumabserrukfc(1,t)/sumtruevalukfc(1,t);
end
UKFTrackLoss=((MCrun-ctcukf)/MCrun)*100;
display(UKFTrackLoss)
%------------------------------------------------------
sumerrghf=zeros(1,timespan/T);
sumerrghfv=zeros(1,timespan/T);
sumabserrghf=zeros(1,timespan/T);
sumtruevalghf=zeros(1,timespan/T);
sumabserrghfv=zeros(1,timespan/T);
sumtruevalghfv=zeros(1,timespan/T);
sumabserrghfc=zeros(1,timespan/T);
sumtruevalghfc=zeros(1,timespan/T);
for k=1:ctcghf
    sumerrghf(1,:)=sumerrghf(1,:)+cerrghf(k,:);
    sumerrghfv(1,:)=sumerrghfv(1,:)+cerrghfv(k,:);
    sumabserrghf(1,:)=sumabserrghf(1,:)+cabserrghf(k,:);
    sumtruevalghf(1,:)=sumtruevalghf(1,:)+ctruevalghf(k,:);
    sumabserrghfv(1,:)=sumabserrghfv(1,:)+cabserrghfv(k,:);
    sumtruevalghfv(1,:)=sumtruevalghfv(1,:)+ctruevalghfv(k,:);
    sumabserrghfc(1,:)=sumabserrghfc(1,:)+cabserrghfc(k,:);
    sumtruevalghfc(1,:)=sumtruevalghfc(1,:)+ctruevalghfc(k,:);
end
rmserrghf=zeros(1,timespan/T);
rmserrghfv=zeros(1,timespan/T);
relerrghf=zeros(1,timespan/T);
relerrghfv=zeros(1,timespan/T);
relerrghfc=zeros(1,timespan/T);
for t=1:timespan/T
    rmserrghf(1,t)=sqrt(sumerrghf(1,t)/ctcghf);
    rmserrghfv(1,t)=sqrt(sumerrghfv(1,t)/ctcghf);
    relerrghf(1,t)=sumabserrghf(1,t)/sumtruevalghf(1,t);
    relerrghfv(1,t)=sumabserrghfv(1,t)/sumtruevalghfv(1,t);
    relerrghfc(1,t)=sumabserrghfc(1,t)/sumtruevalghfc(1,t);
end
GHFTrackLoss=((MCrun-ctcghf)/MCrun)*100;
display(GHFTrackLoss)
%--------------------------------------------------------
sumerrekf=zeros(1,timespan/T);
sumerrekfv=zeros(1,timespan/T);
sumabserrekf=zeros(1,timespan/T);
sumtruevalekf=zeros(1,timespan/T);
sumabserrekfv=zeros(1,timespan/T);
sumtruevalekfv=zeros(1,timespan/T);
sumabserrekfc=zeros(1,timespan/T);
sumtruevalekfc=zeros(1,timespan/T);
for k=1:ctcekf
    sumerrekf(1,:)=sumerrekf(1,:)+cerrekf(k,:);
    sumerrekfv(1,:)=sumerrekfv(1,:)+cerrekfv(k,:);
    sumabserrekf(1,:)=sumabserrekf(1,:)+cabserrekf(k,:);
    sumtruevalekf(1,:)=sumtruevalekf(1,:)+ctruevalekf(k,:);
    sumabserrekfv(1,:)=sumabserrekfv(1,:)+cabserrekfv(k,:);
    sumtruevalekfv(1,:)=sumtruevalekfv(1,:)+ctruevalekfv(k,:);
    sumabserrekfc(1,:)=sumabserrekfc(1,:)+cabserrekfc(k,:);
    sumtruevalekfc(1,:)=sumtruevalekfc(1,:)+ctruevalekfc(k,:);
end
rmserrekf=zeros(1,timespan/T);
rmserrekfv=zeros(1,timespan/T);
relerrekf=zeros(1,timespan/T);
relerrekfv=zeros(1,timespan/T);
relerrekfc=zeros(1,timespan/T);
for t=1:timespan/T
    rmserrekf(1,t)=sqrt(sumerrekf(1,t)/ctcekf);
    rmserrekfv(1,t)=sqrt(sumerrekfv(1,t)/ctcekf);
    relerrekf(1,t)=sumabserrekf(1,t)/sumtruevalekf(1,t);
    relerrekfv(1,t)=sumabserrekfv(1,t)/sumtruevalekfv(1,t);
    relerrekfc(1,t)=sumabserrekfc(1,t)/sumtruevalekfc(1,t);
end
EKFTrackLoss=((MCrun-ctcekf)/MCrun)*100;
display(EKFTrackLoss)
%---------------------------------------------------------
% sumerr2=zeros(1,timespan/T);
% sumerr2v=zeros(1,timespan/T);
% sumabserrsrf=zeros(1,timespan/T);
% sumtruevalsrf=zeros(1,timespan/T);
% sumabserrsrfv=zeros(1,timespan/T);
% sumtruevalsrfv=zeros(1,timespan/T);
% sumabserrsrfc=zeros(1,timespan/T);
% sumtruevalsrfc=zeros(1,timespan/T);
% for k=1:ctc2
%     sumerr2(1,:)=sumerr2(1,:)+cerr2(k,:);
%     sumerr2v(1,:)=sumerr2v(1,:)+cerr2v(k,:);
%     sumabserrsrf(1,:)=sumabserrsrf(1,:)+cabserrsrf(k,:);
%     sumtruevalsrf(1,:)=sumtruevalsrf(1,:)+ctruevalsrf(k,:);
%     sumabserrsrfv(1,:)=sumabserrsrfv(1,:)+cabserrsrfv(k,:);
%     sumtruevalsrfv(1,:)=sumtruevalsrfv(1,:)+ctruevalsrfv(k,:);
%     sumabserrsrfc(1,:)=sumabserrsrfc(1,:)+cabserrsrfc(k,:);
%     sumtruevalsrfc(1,:)=sumtruevalsrfc(1,:)+ctruevalsrfc(k,:);
% end
% rmserr2=zeros(1,timespan/T);
% rmserr2v=zeros(1,timespan/T);
% relerrsrf=zeros(1,timespan/T);
% relerrsrfv=zeros(1,timespan/T);
% relerrsrfc=zeros(1,timespan/T);
% for t=1:timespan/T
%     rmserr2(1,t)=sqrt(sumerr2(1,t)/(ctc2));
%     rmserr2v(1,t)=sqrt(sumerr2v(1,t)/(ctc2));
%     relerrsrf(1,t)=sumabserrsrf(1,t)/sumtruevalsrf(1,t);
%     relerrsrfv(1,t)=sumabserrsrfv(1,t)/sumtruevalsrfv(1,t);
%     relerrsrfc(1,t)=sumabserrsrfc(1,t)/sumtruevalsrfc(1,t);
% end

% figure(4)
% plot(1:timespan/T,rmserr2,'r');
% SRFTrackLoss=((MCrun-ctc2)/MCrun)*100;
% display(SRFTrackLoss)
% disp(T)








%--------------------CRLB--------------------------
%
ExpP0=SumP0/MCrun;

J=(ExpP0)^-1;
%J=P0^-1;
CRLB=zeros(1,timespan/T);
CRLBv=zeros(1,timespan/T);
J1=inv(J);
CRLB(1)=sqrt((J1(1,1))+(J1(2,2)));
CRLBv(1)=sqrt((J1(3,3))+(J1(4,4)));
D11=F'*(inv(Q))*F;
D12=-F'*(inv(Q));
D21=D12';

for t=2:timespan/T
    D22=inv(Q)+(HRH(:,:,t)/MCrun);
    J=D22-(D21*(inv(J+D11))*D12);
    J1=inv(J);
    CRLB(t)=sqrt((J1(1,1))+(J1(2,2)));
    CRLBv(t)=sqrt((J1(3,3))+(J1(4,4)));
end
%M=900/timespan;
figure(2);
plot(1:(timespan/T), CRLB,'k','Linewidth',1.5);
hold on;
hold on;
% plot(1:(timespan/T),rmserr2,'c','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrckf,'r','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrcqkf,'m','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrukf,'b','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrghf,'color',[0.8 0.8 1],'Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrekf,'g','Linewidth',1.5);
ylabel('RMS error of position in km');
xlabel('Time in min');
title('RMS Position Error Plot');
legend('CRLB','CKF','CQKF','UKF','GHF','EKF');

figure(12);
plot(1:(timespan/T), CRLBv,'k','Linewidth',1.5);
hold on;
hold on;
% plot(1:(timespan/T),rmserr2v,'c','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrckfv,'r','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrcqkfv,'m','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrukfv,'b','Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrghfv,'color',[0.8 0.8 1],'Linewidth',1.5);
hold on;
plot(1:timespan/T,rmserrekfv,'g','Linewidth',1.5);
ylabel('RMS error of velocity in km/min');
xlabel('Time in min');
title('RMS Velocity Error Plot');
legend('CRLB','CKF','CQKF','UKF','GHF','EKF');

xtrack=zeros(1,timespan/T);
ytrack=zeros(1,timespan/T);
for t=1:timespan/T
    xtrack(1,t)=(xo(1,t)+xk(1,t));
    ytrack(1,t)=(xo(2,t)+xk(2,t));
end

figure(100)
plot(1:timespan/T,relerrckf*100,'r','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrcqkf*100,'b','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrukf*100,'g','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrekf*100,'color',[0.8 0.8 1],'Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrghf*100,'c','Linewidth',1.5);
hold on;
% plot(1:timespan/T,relerrsrf*100,'m','Linewidth',1.5);
hold on;
ylabel('Relative error of range in %');
xlabel('Time in min');
title('Relative Error of Range Plot');
legend('CKF','CQKF','UKF','EKF','GHF');

figure(101)
plot(1:timespan/T,relerrckfv*100,'r','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrcqkfv*100,'b','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrukfv*100,'g','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrekfv*100,'color',[0.8 0.8 1],'Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrghfv*100,'c','Linewidth',1.5);
hold on;
% plot(1:timespan/T,relerrsrfv*100,'m','Linewidth',1.5);
hold on;
ylabel('Relative error of velocity in %');
xlabel('Time in min');
title('Relative Error in Velocity Plot');
legend('CKF','CQKF','UKF','EKF','GHF');

figure(102)
plot(1:timespan/T,relerrckfc*100,'r','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrcqkfc*100,'b','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrukfc*100,'g','Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrekfc*100,'color',[0.8 0.8 1],'Linewidth',1.5);
hold on;
plot(1:timespan/T,relerrghfc*100,'c','Linewidth',1.5);
hold on;
% plot(1:timespan/T,relerrsrfc*100,'m','Linewidth',1.5);
hold on;
ylabel('Relative error of course in %');
xlabel('Time in min');
title('Relative Error in Course Plot');
legend('CKF','CQKF','UKF','EKF','GHF');

%%
figure(5)
plot(xo(1,:),xo(2,:),'r-.','Linewidth',1.5);
hold on;
plot(xtrack,ytrack,'b--','Linewidth',1.5);
 text(xo(1,1),xo(2,1),'Start');
text(xtrack(1),ytrack(1),'Start');
% plot(xo(1,13),xo(2,13),'k*','MarkerSize', 10,'Linewidth',1)
% plot(xo(1,17),xo(2,17),'k+','MarkerSize', 10,'Linewidth',1.5)
plot(xo(1,15),xo(2,15),'k*','MarkerSize', 10,'Linewidth',1)
% plot(xtrack(1,15),xtrack(2,15),'k+','MarkerSize', 10,'Linewidth',1)
 for i=100:100:timespan
     text(xo(1,i),xo(2,i),int2str(i));
     text(xtrack(i),ytrack(i),int2str(i));
 end
 legend({'Own-ship','Target','Observer Maneuvers'},'fontsize',12);
% legend({'Own-ship','Target','Observer Maneuver Starts','Observer Maneuver Ends'},'fontsize',12);
% title('True Trajectory');
xlabel('x-axis in km','fontsize',16);
ylabel('y-axis in km','fontsize',16);
% ylim([-2,7])
% ratez=zeros(1,timespan/T-1);
% for k=2:timespan/T
%     ratez(1,k-1)=z(1,k)-z(1,k-1);
% end
% zcrop=zeros(1,timespan/T);
% for k=1:timespan/T
%     zcrop(1,k)=z(1,k);
% end

