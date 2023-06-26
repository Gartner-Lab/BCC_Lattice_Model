%% simulate_sorting
% 22 March 2021
% Jennifer Hu
%
% A lattice-based model similar to Matt Thompson's, but in 3D. It uses the
% interfacial tension sets created by my R code to parameterize the
% gammas, and then allows it to evolve over time.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

datadir = find_datadir();
const = struct('L',1,'M',2,'C',3,'B',4,'X',5, ...
    'LL',1,'LM',2,'LX',3,'MM',4,'MX',5);
gammaTable = readtable(fullfile(datadir, 'interfacial tension sets.csv'));

%% Constants and settings
nsteps = 10000; steps = 1:nsteps;
n_iter = 100;
T0 = 500; Tfs = [0, 10, 50, 100, 200]; nTs = length(Tfs);
% this anneals until ~midway through steps. increase numerator for faster
base = 1-10/nsteps;
% Wnames to anneal
sets = ["GFP","GFP_agarose","H1047R"]; nws = length(sets);

%% Set up for tissue
radii = [5,7,8]; nrs = length(radii);
tissuestats = zeros(nrs*nws*n_iter*nTs*(1+nsteps), 4+2+5);

%% Anneal tissue
tic;
% Save: Wname, // Tf, r, n, iter, // T, E, areas(x5)
nrows = 0;
for radius = radii
    tissue = make_tissue(radius, 0.5, const);
    for w = 1:nws
        set = sets(w); idx = strcmp(gammaTable.Wname, set);
        fprintf(set);
        gammas = zeros(1,5);
        gammas(const.LL) = gammaTable.gamma(idx & strcmp(gammaTable.edgename, 'LL'));
        gammas(const.LM) = gammaTable.gamma(idx & strcmp(gammaTable.edgename, 'LM'));
        gammas(const.LX) = gammaTable.gamma(idx & strcmp(gammaTable.edgename, 'LX'));
        gammas(const.MM) = gammaTable.gamma(idx & strcmp(gammaTable.edgename, 'MM'));
        gammas(const.MX) = gammaTable.gamma(idx & strcmp(gammaTable.edgename, 'MX'));
        for Tf = Tfs
            fprintf('.');
            annealing = (T0-Tf) * (base .^ steps) + Tf;
            for iter = 1:n_iter
                tissue = reset_tissue(tissue);
                tissue.gammas = gammas;
                tissue.E = tissue_energy(tissue);
                tissuestats((nrows+1):(nrows+1+nsteps),1:4) = ...
                    repmat([Tf, radius, tissue.n(const.C), iter], nsteps+1, 1);
                tissuestats(nrows+1,5:6) = [T0, tissue.E];
                tissuestats(nrows+1,7:end) = quantify_tissue(tissue);
                nrows = nrows + 1;
                for step = steps
                    tissue = evolve_tissue(tissue, annealing(step));
                    tissuestats(step+nrows,5:6) = [annealing(step), tissue.E];
                    tissuestats(step+nrows,7:end) = quantify_tissue(tissue);
                end
                nrows = nrows+nsteps;
            end
        end
        fprintf('%s\n', tocstring());
    end
end
%%
tissueTable = array2table(tissuestats(1:nrows,:), 'VariableNames', ...
    {'Tf','r','n','i','T','E','LL','LM','LX','MM','MX'}); 
Wnames = repmat(repelem(sets, nTs*n_iter*(nsteps+1))', nrs, 1);
tissueTable = addvars(tissueTable, Wnames(1:nrows), ...
    'Before', 1, 'NewVariableNames', 'Wname');
writetable(tissueTable, fullfile(datadir, 'trajectories578.csv'));
