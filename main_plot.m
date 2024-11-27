function [out] = main()

switch getenv('ENV')
    case 'IUHPC'
        disp('loading paths (HPC)')
        addpath(genpath('/N/u/brlife/git/jsonlab'))
        addpath(genpath('/N/u/brlife/git/geom3d'))
    case 'VM'
        disp('loading paths (VM)')
        addpath(genpath('/usr/local/jsonlab'))
end



config = loadjson('config.json');
evecs = loadjson(config.evecs);


filelist = dir([config.surfaces '/*.vtk']);
for file = 2:size(filelist)
    sprintf(filelist(file).name)
    [V,F] = read_vtk([config.surfaces '/' filelist(file).name]);
    plot_eigenfunction(V', F, evecs.(filelist(file).name(1:end-4)), 2)
end
end
