function [EFFECTS] = compute_effects(parsM,parsY)

%% Preparing variables
Dl = parsM.Dl;m=size(Dl,2);
Bm = parsM.Bm;
By = parsY.By;
Gamma = reshape(parsY.Gamma,m,2);
dl = parsY.dl;

% Direct Effects
EFFECTS.DE_c = By(1) + By(1)*dl;
EFFECTS.DE_s = By(2) + By(2)*dl;

% Exploded indirect effects (centers)
EFFECTS.IE_cc = (Bm(1,:) .* Gamma(:,1)')+(Bm(1,:) .* Gamma(:,1)')*dl;
EFFECTS.IE_cs = (Bm(1,:) .* Gamma(:,2)' .* diag(Dl)')+(Bm(1,:) .* Gamma(:,2)' .* diag(Dl)')*dl;

% Exploded indirect effects (ranges)
EFFECTS.IE_sc = (Bm(2,:) .* Gamma(:,1)')+(Bm(2,:) .* Gamma(:,1)')*dl;
EFFECTS.IE_ss = (Bm(2,:) .* Gamma(:,2)' .* diag(Dl)')+(Bm(2,:) .* Gamma(:,2)' .* diag(Dl)')*dl;

% Aggregate indirect effect
EFFECTS.IE_c = EFFECTS.IE_cc + EFFECTS.IE_cs;
EFFECTS.IE_s = EFFECTS.IE_sc + EFFECTS.IE_ss;

% Create effects table
table(1,1:m) = [EFFECTS.DE_c zeros(1,m-1)];
table(2,1:m) = [EFFECTS.DE_s zeros(1,m-1)];
table(3,1:m) = [EFFECTS.IE_c zeros(1,m-size(EFFECTS.IE_c,2))];
table(4,1:m) = [EFFECTS.IE_s zeros(1,m-size(EFFECTS.IE_s,2))];
EFFECTS.table=table;


end