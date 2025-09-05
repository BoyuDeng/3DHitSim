dot_muall = [];   dot_seall = [];
up_muall  = [];    up_seall  = [];
cr_muall  = [];    cr_seall  = [];
yw_muall  = [];    yw_seall  = [];
Wvals_all = [];
%%

numFields  = numel(results);          % e.g. 5
numSamples = numel(results{1});       % e.g. 13

% Pre‑allocate (use NaN so bugs are visible) ---------------------------
dot_mu      = NaN(numFields,numSamples);
dot_se      = NaN(numFields,numSamples);
up_mu       = NaN(numFields,numSamples);
up_se       = NaN(numFields,numSamples);
cr_mu       = NaN(numFields,numSamples);
cr_se       = NaN(numFields,numSamples);
yw_mu       = NaN(numFields,numSamples);
yw_se       = NaN(numFields,numSamples);
Wvals       = NaN(numFields,numSamples);

for n = 1:numFields
    % ---------- locate field folder ----------
    fn = mod(n,10); if fn==0, fn = 10; end
    fldPath = sprintf('/Users/boyi/Simulations/hit/output_fields%dU=8',fn);

    % ---------- read velocity fields ----------
    nFiles   = end_num - start_num + 1;
    uField = cell(nFiles,1); vField = cell(nFiles,1); wField = cell(nFiles,1);
    for k = 1:nFiles
        fID = start_num + k - 1;
        fName = fullfile(fldPath,sprintf('field_%05d',fID));
        [~,~,~,~,vars] = read_field_file(fName);
        uField{k} = vars(1); vField{k} = vars(2); wField{k} = vars(3);
    end

    %---------- reference RMS speed ----------
    Uref = calculateRMS(uField(1:min(100,nFiles)),vField(1:min(100,nFiles)),wField(1:min(100,nFiles)));

    % ---------- loop over trajectory samples ----------
    for s = 1:numSamples
        tmp = results{n}{s};        % coeffs + W
        Vs  = V14unlim(t,tmp.coeffs,tmp.W) ./ Uref;
        Vunit = Vs ./ vecnorm(Vs);

        X  = X14unlim(t,tmp.coeffs,tmp.W);
        Vel = get_vel(uField,vField,wField,64,X) ./ Uref;  % 3×T

        % ------------- scalar products -------------
        dotAll = sum(Vel .* Vunit,1);
        [dot_mu(n,s),dot_se(n,s)] = mean_and_se(dotAll);

        % ------------- components -------------
        upAll = Vel(3,:);
        crAll = abs(Vel(1,:));
        ywAll = Vel(2,:);
        [up_mu(n,s),up_se(n,s)] = mean_and_se(upAll);
        [cr_mu(n,s),cr_se(n,s)] = mean_and_se(crAll);
        [yw_mu(n,s),yw_se(n,s)] = mean_and_se(ywAll);

        % ------------- W value -------------
        Wvals(n,s) = tmp.W / Uref;
    end
end
% ================= accumulate over all data sets =====================
dot_muall = [dot_muall; dot_mu];   dot_seall = [dot_seall; dot_se];
up_muall  = [up_muall;  up_mu];    up_seall  = [up_seall;  up_se];
cr_muall  = [cr_muall;  cr_mu];    cr_seall  = [cr_seall;  cr_se];
yw_muall  = [yw_muall;  yw_mu];    yw_seall  = [yw_seall;  yw_se];
Wvals_all = [Wvals_all; Wvals];

%% ================= Grand means & propagated SE =======================
N_all = size(dot_muall,1);

[meanDot, seDot] = propagate_se(dot_muall , dot_seall , N_all);
[meanUp , seUp ] = propagate_se(up_muall  , up_seall  , N_all);
[meanCr , seCr ] = propagate_se(cr_muall  , cr_seall  , N_all);
[meanYw , seYw ] = propagate_se(yw_muall  , yw_seall  , N_all);
meanW  = mean(Wvals_all,1);
seW    = std (Wvals_all,0,1) ./ sqrt(size(Wvals_all,1));

%% ======================= Plotting ====================================
plot_with_err(Gall,meanDot,seDot,'o','<w·u_{unit}>','Dot Product');
plot_with_err(Gall,meanUp ,seUp ,'s','<W_z>','Updraft');
plot_with_err(Gall,meanCr ,seCr ,'d','<|W_x|>','Cross‑wind');
plot_with_err(Gall,meanYw ,seYw ,'^','<W_y>','Y‑wind');
plot_with_err(Gall,meanW  ,seW , 'v','<U>','RMS Speed');

%% ======================= helpers =====================================
function [mu,se] = mean_and_se(x)
    x = x(:); n = numel(x);
    mu = mean(x);
    se = std(x,0) / sqrt(n);
end

function [grandMu,grandSe] = propagate_se(muMat,seMat,N)
    grandMu = mean(muMat,1);
    %   Var(μ̄) = Σ σ_i^2 / N^2  → SE = √(Σσ_i^2) / N
    grandSe = sqrt( sum(seMat.^2,1) ) / N;
end

function plot_with_err(G,mu,se,marker,ylbl,ttl)
    figure; errorbar(G,mu,se,[marker,'-'],'LineWidth',1.5);
    set(gca,'XScale','log'); grid on;
    xlabel('G'); ylabel(ylbl); title(ttl);
end
