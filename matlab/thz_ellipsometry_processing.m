%% FePS3 THz ellipsometry - processing script
%
% Processes the three measured angle-of-incidence datasets:
% windowing, FFT, per-bin harmonic fit (Eq. 2.10), bootstrap
% uncertainties, and the isotropic three-layer inversion at
% 70 deg (Eq. 2.15). Plotting code is omitted.

clear; clc; close all;
%% ===== Datasets =====
baseDir = fullfile(pwd, 'data', 'FePS3');   % one sub-folder per angle of incidence
DS(1).name = '60 DEGREE';  DS(1).AOI = 60;
DS(2).name = '70 DEGREE';  DS(2).AOI = 70;
DS(3).name = '75 DEGREE';  DS(3).AOI = 75;
invertAOI = 70;
sampleName = 'FePS3';
outFolder = fullfile(baseDir, 'thesis_figures');
%% ===== Parameters =====
c_light_cm = 2.99792458e10;
c_SI = 299792458;
preScanDist = 3.0;
postScanDist = 3.0;
max_freq_cm = 100;
freq_res_cm = 0.3;
diagMin_THz = 0.15;  diagMax_THz = 1.20;
coreMin_THz = 0.40;  coreMax_THz = 0.70;
dropAngles = [0 360];
d_FePS3_um = 80;
n_Si = 3.418;
sampleSize_mm = 3.0;
spot1THz_mm = 3.44;
nBoot = 300;
nBootNK = 60;
rng(1);
PsiIsoLimit = 45;
if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% ===== Process each AOI =====
for q = 1:numel(DS)
    folderPath = fullfile(baseDir, DS(q).name);
    files = dir(fullfile(folderPath,'*.txt'));
    files = files(~contains(lower({files.name}),'_fft'));
    if isempty(files), error('No .txt files in %s', folderPath); end
    nF = numel(files);
    pkTime = zeros(1,nF); pkAmp = zeros(1,nF);
    rawT0 = cell(1,nF); rawE0 = cell(1,nF); ang0 = zeros(1,nF);
    for k = 1:nF
        dat = readmatrix(fullfile(files(k).folder,files(k).name), ...
            'FileType','text','CommentStyle','#');
        t = dat(:,1); E = dat(:,2);
        ok = isfinite(t) & isfinite(E); t = t(ok); E = E(ok);
        [amp,pk] = max(abs(E));
        rawT0{k} = t; rawE0{k} = E; pkTime(k) = t(pk); pkAmp(k) = amp;
        [~,bn] = fileparts(files(k).name);
        ang0(k) = str2double(regexp(bn,'\d+','match','once'));
    end
    [~,kRef] = max(pkAmp);
    t_ref = pkTime(kRef);
    spec = cell(1,nF); wref = [];
    for k = 1:nF
        [tw,Ew] = WindowTDSDataMenlo(rawT0{k}, rawE0{k}, t_ref, ...
            preScanDist, postScanDist);
        [ta,Ea,~] = AppodizeTDSData(tw, Ew, t_ref);
        [w,Y] = CalcTDSFFT_Menlo_New(ta, Ea, max_freq_cm, freq_res_cm);
        spec{k} = Y(:);
        if isempty(wref), wref = w(:); end
    end
    ECfull = complex(zeros(numel(wref), nF));
    for k = 1:nF, ECfull(:,k) = spec{k}; end
    [angS, si] = sort(ang0);
    ECfull = ECfull(:,si);
    rawT0 = rawT0(si); rawE0 = rawE0(si);
    keep = ~ismember(angS, dropAngles);
    ang = angS(keep);
    EC = ECfull(:,keep);
    F = fitHarmonic(EC, ang);
    rho = (F.alpha + 1) ./ (F.alpha - 1);
    f_THz = (wref * c_light_cm / 1e12).';
    DS(q).f_THz = f_THz;
    DS(q).mDiag = f_THz >= diagMin_THz & f_THz <= diagMax_THz;
    DS(q).mCore = f_THz >= coreMin_THz & f_THz <= coreMax_THz;
    DS(q).ang = ang; DS(q).EC = EC; DS(q).nAng = numel(ang);
    DS(q).a0 = F.a0; DS(q).as = F.as; DS(q).bs = F.bs;
    DS(q).alpha = F.alpha; DS(q).beta = F.beta; DS(q).R2 = F.R2;
    DS(q).rho = rho;
    DS(q).Psi = atand(abs(rho));
    DS(q).Delta = mod(rad2deg(angle(rho)), 360);
    DS(q).t_ref = t_ref;
    DS(q).rawT = rawT0(keep); DS(q).rawE = rawE0(keep);

    % ---- bootstrap over analyzer angles ----
    mC = DS(q).mCore;
    nb = numel(ang);
    ECc = EC(mC,:);
    PsiB = nan(nBoot, sum(mC));
    for b = 1:nBoot
        idx = randi(nb, 1, nb);
        Fb = fitHarmonic(ECc(:,idx), ang(idx));
        PsiB(b,:) = atand(abs((Fb.alpha + 1) ./ (Fb.alpha - 1)));
    end
    DS(q).PsiCI = prctile(PsiB, [16 84], 1);
    DS(q).PsiSE = std(PsiB, 0, 1);
    DS(q).rhoBoot = [];
    if DS(q).AOI == invertAOI
        DS(q).rhoBoot = complex(nan(nBootNK, sum(mC)));
        for b = 1:nBootNK
            idx = randi(nb, 1, nb);
            Fb = fitHarmonic(ECc(:,idx), ang(idx));
            DS(q).rhoBoot(b,:) = (Fb.alpha + 1) ./ (Fb.alpha - 1);
        end
    end
end

%% ===== Inversion at 70 deg, with bootstrap propagation =====
qi = find([DS.AOI] == invertAOI, 1);
theta = deg2rad(DS(qi).AOI);
d_SI = d_FePS3_um*1e-6;
fCore = DS(qi).f_THz(DS(qi).mCore);
rhoC = DS(qi).rho(DS(qi).mCore);
nC = numel(fCore);
fsOpts = optimoptions('fsolve','Display','off', ...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);
starts = [3.5 0.5; 2.5 0.2; 4.5 1.0; 3.0 2.0; 4.0 0.1];
[n_f, k_f] = invertBandContinuous(fCore, rhoC, n_Si, theta, d_SI, c_SI, ...
    starts, fsOpts);

%% ===================== local functions =====================
function F = fitHarmonic(E, angles)
    A = deg2rad(angles(:));
    X = [ones(numel(A),1), cos(2*A), sin(2*A)];
    cr = X \ real(E)';
    ci = X \ imag(E)';
    F.a0 = cr(1,:) + 1i*ci(1,:);
    F.as = cr(2,:) + 1i*ci(2,:);
    F.bs = cr(3,:) + 1i*ci(3,:);
    F.alpha = F.as ./ F.a0;
    F.beta = F.bs ./ F.a0;
    F.model = F.a0.' + F.as.'*cos(2*A).' + F.bs.'*sin(2*A).';
    F.R2 = (1 - sum(abs(E-F.model).^2,2) ./ sum(abs(E-mean(E,2)).^2,2)).';
end

function [nOut, kOut] = invertBandContinuous(fC, rhoC, n3, theta, d, c, ...
    starts, opts)
    N = numel(fC);
    nOut = nan(1,N); kOut = nan(1,N);
    carry = starts(1,:);
    for i = 1:N
        if ~isfinite(rhoC(i)), continue; end
        fun = @(x) rhoResid(x, n3, theta, d, fC(i)*1e12, c, rhoC(i));
        trials = [carry; starts];
        best = []; bestScore = Inf;
        for s = 1:size(trials,1)
            try
                [xs, fv] = fsolve(fun, trials(s,:), opts);
                if xs(1) > 1 && xs(2) >= 0 && norm(fv) < 1e-6
                    if i == 1, sc = norm(xs - starts(1,:));
                    else, sc = norm(xs - carry);
                    end
                    if sc < bestScore, best = xs; bestScore = sc; end
                end
            catch
            end
        end
        if ~isempty(best)
            nOut(i) = best(1); kOut(i) = best(2); carry = best;
        end
    end
end

function F = rhoResid(x, n3, theta1, d, f_Hz, c_SI, rho_meas)
    n2 = x(1) + 1i*x(2);
    sinT1 = sin(theta1); cosT1 = cos(theta1);
    cosT2 = sqrt(1 - (sinT1/n2)^2);
    cosT3 = sqrt(1 - (sinT1/n3)^2);
    r12p = (n2*cosT1 - cosT2)/(n2*cosT1 + cosT2);
    r12s = (cosT1 - n2*cosT2)/(cosT1 + n2*cosT2);
    r23p = (n3*cosT2 - n2*cosT3)/(n3*cosT2 + n2*cosT3);
    r23s = (n2*cosT2 - n3*cosT3)/(n2*cosT2 + n3*cosT3);
    ph = exp(1i*2*(2*pi*f_Hz*d/c_SI)*n2*cosT2);
    r_p = r12p + (1 - r12p^2)*r23p*ph;
    r_s = r12s + (1 - r12s^2)*r23s*ph;
    e = (r_p/r_s) - rho_meas;
    F = [real(e), imag(e)];
end
