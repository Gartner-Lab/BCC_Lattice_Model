function [datadir] = find_datadir()
    if exist('/Volumes', 'dir') > 0
        hddir = fullfile('/Volumes', 'JenniferHD3', 'DAH Data');
        if isfolder(hddir)
            datadir = hddir;
        else
            datadir = fullfile('/Users', 'Vasudha', 'Library','CloudStorage', 'Box-Box', ...
                             'Gartnerlab Data', 'Individual Folders','Jennifer Hu', ...
                             'Analysis','DAH Modeling','DAH Data','BCC','tissues');
        end
    elseif isfolder(fullfile('C:', 'Users', 'Gartner'))
        datadir = fullfile('C:', 'Users', 'Gartner', 'Box', ...
                     'Gartnerlab Data', 'Individual Folders', ...
                     'Jennifer Hu','Data','Lattice');
    end
    assert(exist(datadir, 'dir') > 0, 'Data directory should exist.');
end

Analysis/DAH Modeling/DAH Data/BCC/tissues