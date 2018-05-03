function [CI] = BCA_CI(s_boot,trim,alpha)

parStar = nanmean(s_boot);
s_boot(isnan(s_boot))=[];

if length(s_boot) > 2
    if nanmean(s_boot) == NaN
        CI = [0 0];
    end
    if nanmean(s_boot) ~= 0
        s_boot(isnan(s_boot)) = [];
        Sboot = sort(s_boot);
        q=length(Sboot);
        b = sum(Sboot < parStar)/q;
        a = sum((trimmean(Sboot,trim) - Sboot).^3) / (6 * (sum((trimmean(Sboot,trim) - Sboot).^2)^1.5));
        zb = norminv(b,0,1);
        zalpha2 = norminv((alpha/2),0,1);
        zalpha21 = norminv(1-(alpha/2),0,1);
        r = round(q * normcdf(zb + ((zb + zalpha2) / (1 - a * (zb + zalpha2))),0,1));
        if r<=0
            r=1;
        end
        s = round(q * normcdf(zb + ((zb + zalpha21) / (1 - a * (zb + zalpha21))),0,1));
        CI = [Sboot(r) Sboot(s)];
    else
        CI = [0 0];
    end
else
    CI = [0 0];
end