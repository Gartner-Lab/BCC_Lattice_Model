datadir = '' % Edit file path;
%datadir = find_datadir();
const = struct('L',1,'M',2,'C',3,'B',4,'X',5,'EL',6, ...
    'LL',1,'LM',2,'LX',3,'MM',4,'MX',5);

%% Quantify boundary sites ------------------------------------------------
% This section makes a bunch of tissues and writes down how many different
% boundary site types there are (classified by total exposed area). I used
% it to calculate multiplicity.
tic
radii = 5:2:19; nrows = 0; data = zeros(length(radii)*2,5);LF=0.5;
for radius = radii
    tissue = make_tissue_cube(radius, LF, const);
    % where are the CX edges? exclude XC.
    CXareas = sum(tissue.edges.areas .* ((tissue.edges.types(:,:,const.C) | tissue.edges.types(:,:,const.X)) ...
        & (tissue.is(:,const.B) & tissue.is(:,const.X)')),2);
    CXareas = round(CXareas(tissue.is(:,const.C)),6,'significant');
    % categorize boundary cells by amount of CX edges
    % Quantifying the number of cells with different ECM edge areas
    [CXareacounts, CXareavals] = groupcounts(CXareas); nadd = length(CXareacounts);
    % tissue stats
    nC = sum(tissue.is(:,const.C)); nB = sum(tissue.is(:,const.B));
    if nrows + nadd > size(data,1)
        data = [data; zeros(nrows+nadd, 5)];
    end
    data((nrows+1):(nrows+nadd),:) = [repmat([radius, nC, nB],nadd,1), CXareavals, CXareacounts];
    nrows = nrows + nadd;
end
dataTable = array2table(data(1:nrows,:),'VariableNames',...
    {'r','n','nB','edgeArea','nC'});
writetable(dataTable, fullfile(datadir, 'bcc_areacounts.csv'));

radii = 5:2:19; nrows = 0; data = zeros(length(radii)*2,5);LF=0.5;
for radius = radii
    tissue = make_tissue_2Dhex(radius, LF, const);
    % where are the CX edges? exclude XC.
    CXareas = sum(tissue.edges.areas .* ((tissue.edges.types(:,:,const.C) | tissue.edges.types(:,:,const.X)) ...
        & (tissue.is(:,const.B) & tissue.is(:,const.X)')),2);
    CXareas = round(CXareas(tissue.is(:,const.C)),6,'significant');
    % categorize boundary cells by amount of CX edges
    % Quantifying the number of cells with different ECM edge areas
    [CXareacounts, CXareavals] = groupcounts(CXareas); nadd = length(CXareacounts);
    % tissue stats
    nC = sum(tissue.is(:,const.C)); nB = sum(tissue.is(:,const.B));
    if nrows + nadd > size(data,1)
        data = [data; zeros(nrows+nadd, 5)];
    end
    data((nrows+1):(nrows+nadd),:) = [repmat([radius, nC, nB],nadd,1), CXareavals, CXareacounts];
    nrows = nrows + nadd;
end
dataTable = array2table(data(1:nrows,:),'VariableNames',...
    {'r','n','nB','edgeArea','nC'});
writetable(dataTable, fullfile(datadir, 'hex_areacounts.csv'));


%% Random sampling multiplicity ------------------------------------------
% Default: R=9 (134 boundary, 147 core) or R=7 (66 boundary, 33 core) or 
% Generate a large number of tissues (1e6) by random placement of LEP and MEP
% Calculate the distribution of edge LEPs
tic
LFs = [0.5, 0.25, 0.75]; % LEP fraction
n_iter = 1e5; % number of tissues to make
radii = [9]; 
andXS = false;
tissuevarnames = {'r','n','LF','nB','iteration', ...
    'nEL','LL_sq','LL_hex','LM_sq','LM_hex','LX_sq','LX_hex',...
    'MM_sq','MM_hex','MX_sq','MX_hex'};

% Save: r, n, LF, iteration // areas(x10)
for radius = radii
   tissue = make_tissue_cube(radius, LFs(1), const);
   for LF = LFs
       tissuestats = zeros(n_iter, 5+11);
       tissuestats(:, 1:5) = ...
            [repmat([radius, tissue.n(const.C), LF,tissue.n(const.B)], n_iter, 1), (1:n_iter)'];
       fprintf('Starting to quantify %d random tissues', n_iter);
       for iter = 1:n_iter 
            tissue = reset_tissue(tissue, LF);
            stats = quantify_tissue(tissue, andXS);
            tissuestats(iter, 6:end) = stats;
            if mod(iter, n_iter/100) == 0
                fprintf('.');
            end
        end
        fprintf('%s\n', tocstring());
        tissueTable = array2table(tissuestats, 'VariableNames', ...
            tissuevarnames); 
        writetable(tissueTable, fullfile(datadir, sprintf('randoms_r%d_LF%g_n%d.csv',radius, LF, n_iter)));

   end    
end

%% Sample tissues with defined LEP boundary occupancy ------------------------------------------
% Make default tissue (R=9, LF=0.5)
tic
n_iter = 1e2;
radii =[7];
LFs = [0.5];
andXS = false;
tissuevarnames = {'r','n','LF','nB','iteration', ...
    'nEL','LL_sq','LL_hex','LM_sq','LM_hex','LX_sq','LX_hex',...
    'MM_sq','MM_hex','MX_sq','MX_hex'};
%tELs = [0.1:0.1:0.9];
tELs = [0];
for radius = radii
    % make tissue
    tissue = make_tissue_cube(radius, LFs(1), const);
    c = tissue.const;
    for LF = LFs
        % Calculate Lb_min and Lb_max
        tissue = reset_tissue(tissue, LF);
        nC = tissue.n(c.C);
        nL = tissue.n(c.L);
        nB = tissue.n(c.B);
        tEL_min = ceil(max(0,nL-nC+nB)*10/nB)/10;
        tEL_max = floor(min(nB, nL)*10/nB)/10;
        % get new range of tELs in accessible range
        tELs = [tEL_min:0.1:tEL_max];
        % intitialize array of #tEL*n_iter rows
        tissuestats = zeros(n_iter*length(tELs), 5+11);
        tissuestats(:, 1:5) = ...
            [repmat([radius, tissue.n(const.C), LF,tissue.n(const.B)], n_iter*length(tELs), 1), (repmat(1:n_iter,1,length(tELs)))'];
        fprintf('Quantify tissues with LF %d', LF);
        fprintf('%s\n', tocstring()); 
        counter = 0;
        for tEL = tELs
            % Make tissues for given tEL  
                fprintf('Quantifying %d random tissues with EL %d', n_iter,tEL);
            for iter = 1:n_iter 
                    tissue = reset_tissue(tissue, LF, tEL);
                    stats = quantify_tissue(tissue, andXS);
                    tissuestats(counter*n_iter+iter, 6:end) = stats;
                    if mod(iter, n_iter/10) == 0
                        fprintf('.');
                    end
            end
            counter = counter + 1;
            fprintf('%s\n', tocstring());    
        end % tEL loop
        tissueTable = array2table(tissuestats, 'VariableNames', ...
                tissuevarnames); 
            writetable(tissueTable, fullfile(datadir, sprintf('targeted_r%d_LF%g_n%d.csv',radius, LF, n_iter)));
    end % LF loop
end % radius loop


%% 2D hexagonal lattice Sample tissues with defined LEP boundary occupancy ------------------------------------------
% Make default tissue (R=5, LF=0.5)
tic
n_iter = 1e3;
radii =[7];
LFs = [0.5];
andXS = false; is2D = true;
tissuevarnames = {'r','n','LF','nB','iteration', ...
    'nEL','LL_area','LM_area','LX_area',...
    'MM_area','MX_area'};
%tELs = [0.1:0.1:0.9];
tELs = [0:0.25:0.75];
for radius = radii
    % make tissue
    tissue = make_tissue_2Dhex(radius, LFs(1), const);
    c = tissue.const;
    for LF = LFs
        % Calculate Lb_min and Lb_max
        tissue = reset_tissue(tissue, LF);
        nC = tissue.n(c.C);
        nL = tissue.n(c.L);
        nB = tissue.n(c.B);
        tEL_min = ceil(max(0,nL-nC+nB)*10/nB)/10;
        tEL_max = floor(min(nB, nL)*10/nB)/10;
        % get new range of tELs in accessible range
        tELs = [tEL_min:0.1:tEL_max];
        % intitialize array of #tEL*n_iter rows
        tissuestats = zeros(n_iter*length(tELs), 5+6);
        tissuestats(:, 1:5) = ...
            [repmat([radius, tissue.n(const.C), LF,tissue.n(const.B)], n_iter*length(tELs), 1), (repmat(1:n_iter,1,length(tELs)))'];
        fprintf('Quantify tissues with LF %d', LF);
        fprintf('%s\n', tocstring()); 
        counter = 0;
        for tEL = tELs
            % Make tissues for given tEL  
                fprintf('Quantifying %d random tissues with EL %d', n_iter,tEL);
            for iter = 1:n_iter 
                    tissue = reset_tissue(tissue, LF, tEL);
                    stats = quantify_tissue(tissue, andXS, is2D);
                    tissuestats(counter*n_iter+iter, 6:end) = stats;
                    if mod(iter, n_iter/10) == 0
                        fprintf('.');
                    end
            end
            counter = counter + 1;
            fprintf('%s\n', tocstring());    
        end % tEL loop
        tissueTable = array2table(tissuestats, 'VariableNames', ...
                tissuevarnames); 
            writetable(tissueTable, fullfile(datadir, sprintf('hex_targeted_r%d_LF%g_n%d.csv',radius, LF, n_iter)));
    end % LF loop
end % radius loop

%% Examples of 2D hexagonal lattice tissues with defined LEP boundary occupancy ------------------------------------------
% Make default tissue (R=5, LF=0.5)
tic
n_iter = 10;
radii =[7];
LFs = [0.55];
andXS = false; is2D = true;
tissuevarnames = {'r','n','LF','nB','iteration', ...
    'nEL','LL_area','LM_area','LX_area',...
    'MM_area','MX_area'};
%tELs = [0.1:0.1:0.9];
tELs = [0:0.1:0.5];
for radius = radii
    % make tissue
    tissue = make_tissue_2Dhex(radius, LFs(1), const);
    c = tissue.const;
    for LF = LFs
        % Calculate Lb_min and Lb_max
        tissue = reset_tissue(tissue, LF);
        nC = tissue.n(c.C);
        nL = tissue.n(c.L);
        nB = tissue.n(c.B);
        tEL_min = ceil(max(0,nL-nC+nB)*10/nB)/10;
        tEL_max = floor(min(nB, nL)*10/nB)/10;
        % get new range of tELs in accessible range
        %tELs = [tEL_min:0.1:tEL_max];
        tELs = [0,0.1,0.5];
        % intitialize array of #tEL*n_iter rows
        tissuestats = zeros(n_iter*length(tELs), 5+6);
        tissuestats(:, 1:5) = ...
            [repmat([radius, tissue.n(const.C), LF,tissue.n(const.B)], n_iter*length(tELs), 1), (repmat(1:n_iter,1,length(tELs)))'];
        fprintf('Quantify tissues with LF %d', LF);
        fprintf('%s\n', tocstring()); 
        counter = 0;
        for tEL = tELs
            % Make tissues for given tEL  
                fprintf('Quantifying %d random tissues with EL %d', n_iter,tEL);
            for iter = 1:n_iter 
                    tissue = reset_tissue(tissue, LF, tEL);
                    stats = quantify_tissue(tissue, andXS, is2D);
                    tissuestats(counter*n_iter+iter, 6:end) = stats;
                    if mod(iter, n_iter/10) == 0
                        fprintf('.');
                    end
                    plot_idx2D(tissue);
                    exportgraphics(gcf, fullfile(datadir,'examples', ...
                        sprintf('hex_r%d_LF%g_EL%g_n%d.pdf',radius, LF, tEL, iter)), 'ContentType', 'vector');
            end
            counter = counter + 1;
            fprintf('%s\n', tocstring());    
        end % tEL loop
        tissueTable = array2table(tissuestats, 'VariableNames', ...
                tissuevarnames); 
            writetable(tissueTable, fullfile(datadir, sprintf('hex_r%d_LF%g_n%d.csv',radius, LF, n_iter)));
    end % LF loop
end % radius loop
 
%% Make scatter plot with centers as circles
c = tissue.const; 
colorL = [255, 185, 15] / 255; colorM = [100, 30, 100] / 255; colorX = [255, 255, 255] / 255;
colors = colorM .* tissue.is(:,c.M) + colorL .* tissue.is(:,c.L) + colorX .* tissue.is(:,c.X);
colors = colors(~tissue.is(:,c.X),:);
figure(1);
clf; 
ax = gca; ax.XColor = 'none'; ax.YColor = 'none'; ax.ZColor = 'none';
    ax.Clipping = 'off';
set(gcf,'position',[20,20,100,100])
cells = tissue.p(~tissue.is(:,c.X),:);
scatter(cells(:,1), cells(:,2), 60,colors,"filled");
pbaspect([1 1 1]); daspect([1 1 1]); 
axis off;
%scatter(tissue.p(:,1), tissue.p(:,2), 60,colors,"filled");
exportgraphics(gcf, fullfile(datadir,'examples', ...
                        sprintf('test_r%d_LF%g_EL%g_n%d.pdf',radius, LF, tEL, iter)), 'ContentType', 'vector');
