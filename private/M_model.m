function [results] = M_model(xc,xl,MC,ML,plotFN,maxIter,eps,txt,startingValues)

%% Preparing variables
n=size(xc,1);
m=size(MC,2);
One=ones(n,m);
Im=eye(m,m);
iter=0;fail=0;
if sum(xl)~=0 || sum(ML(:))~=0
    X=[xc xl];
    M=[MC ML];
else
    X=[xc];
    M=[MC];
end

%% Generating Starting Values
if isstruct(startingValues)==0
    Ac = diag(diag(random('normal',0,1,m,m)));
    Al = diag(diag(random('normal',0,1,m,m)));
    Dl = diag(diag(random('normal',0,1.5,m,m)));
    Bm = random('normal',0,2,size(X,2),m);
else
    Ac = startingValues.Ac;
    Al = startingValues.Al;
    Dl = startingValues.Dl;
    Bm = startingValues.Bm;
end

%% ALS loop
for i=1:maxIter
    Ac_old=Ac;Al_old=Al;Dl_old=Dl;Bm_old=Bm;
    
    %%% Computing alternating estimations
    Al = inv(diag(diag(kron(Im,One'*One)))) * diag(diag(kron(Im,One)'*vec(ML-X*Bm*Dl-One*Ac*Dl)));
    Al = diag(diag(reshape(Al,m,m)));
    
    Dl = inv(diag(diag(...
        kron(Im,Bm'*X'*X*Bm) + kron(Im,Bm'*X'*One*Ac) + kron(Im,Ac'*One'*One*Ac) + kron(Im,Ac'*One'*X*Bm) ))) * diag(diag(...
        kron(Im,X*Bm)'*vec(ML-One*Al) + kron(Im,One*Ac)'*vec(ML-One*Al) ));
    Dl = diag(diag(reshape(Dl,m,m)));
    
    Bm = inv(kron(Im,X'*X) + kron(Dl^2,X'*X)) * (kron(Im,X)'*vec(MC-One*Ac) + kron(Dl,X)'*vec(ML-One*Al-One*Ac*Dl) );
    Bm = reshape(Bm,size(X,2),m);
    
    Ac = inv(diag(diag(kron(Im,One'*One) + kron(Dl^2,One'*One)))) * diag(diag(...
        kron(Im,One)'*vec(MC-X*Bm) + kron(Dl,One)'*vec(ML-One*Al-X*Bm*Dl) ));
    Ac = diag(diag(reshape(Ac,m,m)));
    
    Al = inv(diag(diag(kron(Im,One'*One)))) * diag(diag(kron(Im,One)'*vec(ML-X*Bm*Dl-One*Ac*Dl)));
    Al = diag(diag(reshape(Al,m,m)));

    %%% Computing convergence criteria
    theta=[diag(Al);diag(Dl);diag(Ac);Bm(:)];
    theta_old=[diag(Al_old);diag(Dl_old);diag(Ac_old);Bm_old(:)];
    Table_theta(i,:) = [norm(Ac_old-Ac)^2 norm(Al_old-Al)^2 norm(Dl_old-Dl)^2 norm(Bm_old-Bm)^2];
    D2_theta(i) = norm(theta_old-theta)^2;
    mean_theta(i) = mean(theta_old-theta);
    var_theta(i) = var(theta_old-theta);
    
    %%% Stopping rule
        if i>15 || i>maxIter
            delta(i)=D2_theta(i-10)-D2_theta(i);
            if abs(delta(i))<=eps
                iter=i;
                break
            end
        end
end

%% Plotting convergence
switch plotFN
    case 1
    figure();
    subplot(2,2,1:2);plot(D2_theta,'r');ylim([-1 max(D2_theta)]);title('D^2(theta\_old-theta)');
    subplot(2,2,3);plot(mean_theta);title('mean(theta\_old-theta)')
    subplot(2,2,4);plot(var_theta);ylim([-0.5 max(var_theta)]);title('var(theta\_old-theta)')
end

%% Computing estimated components
MC_star = One*Ac + X*Bm;
ML_star = One*Al + (One*Ac + X*Bm)*Dl;

%% Goodness of fit
try
    TSSm = norm(MC-One*diag(mean(MC)))^2 + norm(ML-One*diag(mean(ML)))^2;
    RSSm = norm(MC-MC_star)^2 + norm(ML-ML_star)^2;
    ESSm = TSSm-RSSm;
    Rm = 1-(RSSm/TSSm);
    RMSE = [sqrt(mean((MC(:)-MC_star(:)).^2)) sqrt(mean((ML(:)-ML_star(:)).^2))];
end

%% Checking ALS convergence
if iter==0 || Rm < 0
    fail=1;
    TSSy=0;RSSy=0;ESSy=0;Rm=0;
    if txt==1 
        disp('* Error: MODEL M: ALS did not converge at a global/local minimum *') 
    end
else
    if txt==1 
        disp('* MODEL M: ALS converged at a global/local minimum *') 
    end
end

%% Saving final data
results.fail = fail;
results.TotIter = iter;
results.R2 = Rm;
results.RMSE = RMSE;
results.convergence.table = Table_theta;
results.convergence.D2_theta = D2_theta;
results.convergence.mean_theta = mean_theta;
results.convergence.var_theta = var_theta;
results.convergence.delta = delta;
results.pars.Al = Al;
results.pars.Ac = Ac;
results.pars.Dl = Dl;
results.pars.Bm = Bm;
results.pars.MC_star = MC_star;
results.pars.ML_star = ML_star;






end