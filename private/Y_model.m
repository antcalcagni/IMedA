function [results] = Y_model(xc,xl,MC,ML,yc,yl,plotFN,maxIter,eps,txt,startingValues)

%% Preparing variables
n=size(xc,1);
m=size(MC,2);
one=ones(n,1);
iter=0;fail=0;
if sum(xl)~=0
    X=[xc xl];
else
    X=[xc];
end

%% Generating Starting Values
if isstruct(startingValues)==0
    ac = random('normal',0,1,1,1);
    al = random('normal',0,1,1,1);
    dl = random('normal',0,1,1,1);
    By = random('normal',0,2,size(X,2),1);
    GammaC = random('normal',0,1,m,1);
    GammaL = random('normal',0,1,m,1);
else
    ac = startingValues.ac;
    al = startingValues.al;
    dl = startingValues.dl;
    By = startingValues.By;
    GammaC = startingValues.Gamma(1:m);
    GammaL = startingValues.Gamma(m+1:end);
end

XC=xc;
XL=xl;

%% ALS loop
for i=1:maxIter
    ac_old=ac;al_old=al;dl_old=dl;By_old=By;GammaC_old=GammaC;GammaL_old=GammaL;
    
    %%% Computing alternating estimations
    
    al = (1/n) * one'*(yl - one*ac*dl - X*By*dl - MC*GammaC*dl - ML*GammaL*dl);
    
    dl = inv(...
        ac*one'*(one*ac + 2*X*By + 2*MC*GammaC + 2*ML*GammaL) +...
        By'*X'*(X*By + 2*MC*GammaC + 2*ML*GammaL) + ...
        GammaC'*MC'*(MC*GammaC + 2*ML*GammaL) + ...
        GammaL'*ML'*ML*GammaL) * (...
        ac*one'*(yl - one*al) + ...
        By'*X'*(yl - one*al) + ...
        GammaC'*MC'*(yl - one*al) + ...
        GammaL'*ML'*(yl - one*al) );
        
    ac = (1/(n*(1+dl^2))) * (...
        one'*(yc-X*By-MC*GammaC-ML*GammaL) + ...
        one'*(yl-one*al-X*By*dl-MC*GammaC*dl-ML*GammaL*dl)*dl );
    
%     By = inv(X'*X + dl*X'*X*dl) * (...
%         (X'*yc - X'*one*ac - X'*MC*GammaC - X'*ML*GammaL) + ...
%         (X'*yl*dl - X'*one*al*dl - X'*one*ac*dl*dl - X'*MC*GammaC*dl*dl - X'*ML*GammaL*dl*dl) );
    
    Byl = By(2);
    Byc = inv(XC'*XC + dl*XC'*XC*dl) * (...
        XC'*(yc - one*ac - xl*Byl - MC*GammaC - ML*GammaL) + ...
        XC'*(yl - one*al - one*ac*dl - xl*Byl*dl - MC*GammaC*dl - ML*GammaL*dl)*dl );
    
    Byl = inv(XL'*XL + dl*XL'*XL*dl) * (...
        XL'*(yc - one*ac - xc*Byc - MC*GammaC - ML*GammaL) + ...
        XL'*(yl - one*al - one*ac*dl - xc*Byc*dl - MC*GammaC*dl - ML*GammaL*dl)*dl );    
    By = [Byc;Byl];
    
    GammaC = inv(MC'*MC + dl*MC'*MC*dl) * (...
        MC'*(yc - one*ac - X*By - ML*GammaL) + ...
        MC'*(yl - one*al - one*ac*dl - X*By*dl - ML*GammaL*dl)*dl );    

    GammaL = inv(ML'*ML + dl*ML'*ML*dl) * (...
        ML'*(yc - one*ac - X*By - MC*GammaC) + ...
        ML'*(yl - one*al - one*ac*dl - X*By*dl - MC*GammaC*dl)*dl );    
    
    Gamma = [GammaC;GammaL];Gamma_old = [GammaC_old;GammaL_old];
   
    
    %%% Computing convergence criteria
    theta=[al;dl;ac;By;Gamma];theta_old=[al_old;dl_old;ac_old;By_old;Gamma_old];
    %[theta theta_old]
    Table_theta(i,:) = [norm(ac_old-ac)^2 norm(al_old-al)^2 norm(dl_old-dl)^2 norm(By_old-By)^2 norm(Gamma_old-Gamma)^2];
    D2_theta(i) = norm(theta_old-theta)^2;
    mean_theta(i) = mean(theta_old-theta);
    var_theta(i) = var(theta_old-theta);
    
    %%% Stopping rule
    if i>15 || i>maxIter
        delta(i)=D2_theta(i-5)-D2_theta(i);
        if abs(delta(i))<=eps
            iter=i;
            break
        end
    end
    
end

%% Plotting convergence
switch plotFN
    case 1
    %close all
    figure();
    subplot(2,2,1:2);plot(D2_theta,'r');ylim([-1 max(D2_theta)]);title('D^2(theta\_old-theta)');
    subplot(2,2,3);plot(mean_theta);title('mean(theta\_old-theta)')
    subplot(2,2,4);plot(var_theta);ylim([-0.5 max(var_theta)]);title('var(theta\_old-theta)')
end

%% Computing estimated components
yc_star = one*ac + X*By + MC*GammaC + ML*GammaL;
yl_star = one*al + (one*ac + X*By + MC*GammaC + ML*GammaL)*dl;

%% Goodness of fit
try
    TSSy = norm(yc-one*mean(yc))^2 + norm(yl-one*mean(yl))^2;
    RSSy = norm(yc-yc_star)^2 + norm(yl-yl_star)^2;
    ESSy = TSSy-RSSy;
    Ry = 1-(RSSy/TSSy);
    RMSE = [sqrt(mean((yc - yc_star).^2)) sqrt(mean((yl - yl_star).^2))];
    %[TSSy RSSy ESSy Ry RMSE]
end

%% Checking ALS convergence
if iter==0 || Ry < 0
    fail=1;
    TSSy=0;RSSy=0;ESSy=0;Ry=0;
    if txt==1 
        disp('* Error: MODEL Y: ALS did not converge at a global/local minimum *') 
    end
else
    if txt==1 
        disp('* MODEL Y: ALS converged at a global/local minimum *') 
    end
end

%% Saving final data
results.fail = fail;
results.TotIter = iter;
results.R2 = Ry;
results.RMSE = RMSE;
results.convergence.table = Table_theta;
results.convergence.D2_theta = D2_theta;
results.convergence.mean_theta = mean_theta;
results.convergence.var_theta = var_theta;
results.convergence.delta = delta;
results.pars.al = al;
results.pars.ac = ac;
results.pars.dl = dl;
results.pars.By = By;
results.pars.Gamma = Gamma;
results.pars.yc_star = yc_star;
results.pars.yl_star = yl_star;





end