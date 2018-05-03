function [results,boot] = IMedA_ALS(xc,xl,MC,ML,yc,yl,plot,maxIter,eps,txt,bootstrap,startingValues)
close all
warning off
boot=0;
err=0;
disp('==============================');disp('** ALS IMedA started **');disp('==============================');

%% Model Y and M estimations
disp('  ');disp('* Computing model M');
[modelM] = M_model(xc,xl,MC,ML,plot(1),maxIter(1),eps(1),txt,startingValues);
disp('  ');disp('* Computing model Y');
[modelY] = Y_model(xc,xl,MC,ML,yc,yl,plot(2),maxIter(2),eps(2),txt,startingValues);

%% Decomposition of effects
if modelM.fail==0 && modelY.fail==0
    EFFECTS = compute_effects(modelM.pars,modelY.pars);
end

%% Variance decomposition 
if modelM.fail==0 && modelY.fail==0
    Gamma = reshape(modelY.pars.Gamma,size(MC,2),2);
    one=ones(size(MC,2),1);
    
    % Decomposition of var(yc)
    c1 = covar(yc,xc*modelY.pars.By(1)); %DEc_c
    c2 = covar(yc,xl*modelY.pars.By(2)); %DEc_s
    c3 = covar(yc,xc*(modelM.pars.Bm(1,:) .* Gamma(:,1)')*one); %IEc_cc
    c4 = covar(yc,xc*(modelM.pars.Bm(1,:) .* Gamma(:,2)' .* diag(modelM.pars.Dl)')*one); %IEc_cs
    c5 = covar(yc,xl*(modelM.pars.Bm(2,:) .* Gamma(:,1)')*one); %IEc_sc
    c6 = covar(yc,xl*(modelM.pars.Bm(2,:) .* Gamma(:,2)' .* diag(modelM.pars.Dl)')*one); %IEc_ss
    c7 = covar(yc,(yc-modelY.pars.yc_star)); %resid
    c8 = covar(yc,(MC-modelM.pars.MC_star)*Gamma(:,1)); %resid
    c9 = covar(yc,(ML-modelM.pars.ML_star)*Gamma(:,2)); %resid

    % Decomposition of var(yl)
    s1 = covar(yl,modelY.pars.dl*xc*modelY.pars.By(1)); %DEs_c
    s2 = covar(yl,modelY.pars.dl*xl*modelY.pars.By(2)); %DEs_s
    s3 = covar(yl,modelY.pars.dl*xc*(modelM.pars.Bm(1,:) .* Gamma(:,1)')*one); %IEs_cc
    s4 = covar(yl,modelY.pars.dl*xc*(modelM.pars.Bm(1,:) .* Gamma(:,2)' .* diag(modelM.pars.Dl)')*one); %IEs_cs
    s5 = covar(yl,modelY.pars.dl*xl*(modelM.pars.Bm(2,:) .* Gamma(:,1)')*one); %IEs_sc
    s6 = covar(yl,modelY.pars.dl*xl*(modelM.pars.Bm(2,:) .* Gamma(:,2)' .* diag(modelM.pars.Dl)')*one); %IEs_ss
    s7 = covar(yl,(yl-modelY.pars.yl_star)); %resid
    s8 = covar(yl,modelY.pars.dl*(MC-modelM.pars.MC_star)*Gamma(:,1)); %resid
    s9 = covar(yl,modelY.pars.dl*(ML-modelM.pars.ML_star)*Gamma(:,2)); %resid

    % Computing decomposition weights 
    %[c1 c2 c3 c4 c5 c6 c7 c8 c9 sum([c1 c2 c3 c4 c5 c6 c8 c9]) sum([c1 c2 c3 c4 c5 c6 c7 c8 c9]) var(yc)]
    %[s1 s2 s3 s4 s5 s6 s7 s8 s9 sum([s1 s2 s3 s4 s5 s6 s8 s9]) sum([s1 s2 s3 s4 s5 s6 s7 s8 s9]) var(yl)]
    u1 = (c1+c2+s1+s2)/(var(yc)+var(yl)); %DE
    u2 = (c3+c4+c5+c6+s3+s4+s5+s6)/(var(yc)+var(yl)); %IE
    u3  = (c7+c8+c9+s7+s8+s9)/(var(yc)+var(yl)); %resid

    w1=(c1+c2+s1+s2)/(var(yc)+var(yl));
    w2=(c3+c4+c5+c6+s3+s4+s5+s6)/(var(yc)+var(yl));
    w3=(c8+c9+s8+s9)/(var(yc)+var(yl));

    t = u1+u2;
    t1 = u1/(u1+u2);
    t2 = u2/(u1+u2);

    FACTORS.DE = w1/(w1+w2+w3);
    FACTORS.IE = w2/(w1+w2+w3);
    FACTORS.other = w3/(w1+w2+w3);
    FACTORS.t = t;
    FACTORS.t1 = t1;
    FACTORS.t2 = t2;
end

%% Descriptive messages
if modelM.fail==1 || modelY.fail==1
    err=1;
    if txt==1
        disp('  ');disp('ATTENTION: Mediation coefficients cannot be computed - both M and Y convergences are required')
    end
    EFFECTS=0;
    FACTORS=0;
end

%% Saving results
results.M = modelM;
results.Y = modelY;
results.EFFECTS = EFFECTS;
results.FACTORS = FACTORS;

disp('  ');disp('** ALS IMedA routine finished **');

if bootstrap(1)==1 && err==0
    [boot] = IMedA_bootstrap(xc,xl,MC,ML,yc,yl,bootstrap(2),maxIter,eps,0)
else
    if err==1
       disp('  ');disp('ATTENTION: Bootstrap cannot be run - both M and Y convergences are required')
    end
end

warning on
end