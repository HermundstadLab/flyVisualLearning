function AOV = performANOVA(rsL,rsNL,rsK90C,rsK90,rsK96C,rsK96,sL,sNL,sK90C,sK90,sK96C,sK96,p)
% PERFORMANOVA performs an ANOVA on the collection of datasets

FactorNames = {"genotype","laser","kir","pattern"};

%construct factors
m = 1;
for i=1:numel(rsL.fly)
    for j=1:numel(rsL.fly(i).trial)
        PI(m,j) = sL.fly(i).trial(j).PI;
        ninds(m,j) = numel(sL.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "WT";
    laser(m,1)    = "on";
    kir(m,1)      = "not expressed";
    pattern(m,1)  = string(rsL.fly(i).punishedBar);
    m = m+1;
end
for i=1:numel(rsNL.fly)
    for j=1:numel(rsNL.fly(i).trial)
        PI(m,j) = sNL.fly(i).trial(j).PI;
        ninds(m,j) = numel(sNL.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "WT";
    laser(m,1)    = "off";
    kir(m,1)      = "not expressed";
    pattern(m,1)  = string(rsNL.fly(i).punishedBar);
    m = m+1;
end
for i=1:numel(rsK90C.fly)
    for j=1:numel(rsK90C.fly(i).trial)
        PI(m,j) = sK90C.fly(i).trial(j).PI;
        ninds(m,j) = numel(sK90C.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "SS00090";
    laser(m,1)    = "on";
    kir(m,1)      = "not expressed";
    pattern(m,1)  = string(rsK90C.fly(i).punishedBar);
    m = m+1;
end
for i=1:numel(rsK90.fly)
    for j=1:numel(rsK90.fly(i).trial)
        PI(m,j) = sK90.fly(i).trial(j).PI;
        ninds(m,j) = numel(sK90.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "SS00090";
    laser(m,1)    = "on";
    kir(m,1)      = "expressed";
    pattern(m,1)  = string(rsK90.fly(i).punishedBar);
    m = m+1;
end
for i=1:numel(rsK96C.fly)
    for j=1:numel(rsK96C.fly(i).trial)
        PI(m,j) = sK96C.fly(i).trial(j).PI;
        ninds(m,j) = numel(sK96C.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "SS00096";
    laser(m,1)    = "on";
    kir(m,1)      = "not expressed";
    pattern(m,1)  = string(rsK96C.fly(i).punishedBar);
    m = m+1;
end
for i=1:numel(rsK96.fly)
    for j=1:numel(rsK96.fly(i).trial)
        PI(m,j) = sK96.fly(i).trial(j).PI;
        ninds(m,j) = numel(sK96.fly(i).trial(j).PIinds);
    end
    genotype(m,1) = "SS00096";
    laser(m,1)    = "on";
    kir(m,1)      = "expressed";
    pattern(m,1)  = string(rsK96.fly(i).punishedBar);
    m = m+1;
end

PInaive = sum(ninds(:,p.naiveTrials ).*PI(:,p.naiveTrials ),2)./sum(ninds(:,p.naiveTrials ),2);
PIprobe = sum(ninds(:,p.probeTrials2).*PI(:,p.probeTrials2),2)./sum(ninds(:,p.probeTrials2),2);
dPI = PIprobe-PInaive;

tbl = table(genotype, laser, pattern, kir, dPI, VariableNames=["genotype","laser","pattern","kir","dPI"]);
tbl = table(genotype, laser, pattern, kir, dPI, VariableNames=["genotype","laser","pattern","kir","dPI"]);
if isMATLABReleaseOlderThan("R2022b") 
    error('This function depends on anova.m, which was introduced to MATLAB starting in their 2022b release.');
else
    AOV = anova(tbl,"dPI ~ genotype + laser + pattern + kir + pattern:kir");
end
    

end
