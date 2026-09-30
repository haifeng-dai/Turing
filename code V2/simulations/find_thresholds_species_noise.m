function thresholds = find_thresholds_species_noise(topo_type, topo_val, cfg)
% 图 2c 噪声组分消融专用阈值搜索；二分规则与 find_thresholds 一致。
rng(0, 'twister');
topo_base = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'topology');
type_u = upper(topo_type);
if ~isfield(cfg, 'L_intra') || isempty(cfg.L_intra)
    switch type_u
        case 'ER', net_file = fullfile(topo_base, 'ER', sprintf('N%d_p%.3f.mat', cfg.N, topo_val));
        case 'WS', net_file = fullfile(topo_base, 'WS', sprintf('N%d_pr%.2f.mat', cfg.N, topo_val));
        case 'BA', net_file = fullfile(topo_base, 'BA', sprintf('N%d_m%d.mat', cfg.N, topo_val));
        case 'SF', net_file = fullfile(topo_base, 'SF', sprintf('N%d_g%.1f.mat', cfg.N, topo_val));
        otherwise, error('未知拓扑类型：%s', topo_type);
    end
    if ~isfile(net_file), error('找不到网络文件：%s', net_file); end
    data = load(net_file, 'nets');
    cfg.L_intra = data.nets;
end
cfg.dt = 0.001;
cfg.steps = 2;
threshold_A = 0.05;
epsilon = 0.05;
if ~isfield(cfg, 'y_seed') || isempty(cfg.y_seed)
    seed_cfg = cfg;
    seed_cfg.K = size(cfg.L_inter,1);
    cfg.y_seed = load_standard_pattern_seed(seed_cfg);
end
if ~isnumeric(cfg.y_seed) || ~isreal(cfg.y_seed) || ...
        numel(cfg.y_seed) ~= 2*cfg.N*size(cfg.L_inter,1) || any(~isfinite(cfg.y_seed(:)))
    error('find_thresholds_species_noise:InvalidSeed', ...
        '斑图初值长度须为 2*N*K，且所有数值有限。');
end
if ~isfield(cfg, 'noise_species_mask') || numel(cfg.noise_species_mask) ~= 2
    error('find_thresholds_species_noise:MissingMask', ...
        '请设置 noise_species_mask=[u_mask,v_mask]。');
end

low = cfg.sigma_min;
high = cfg.sigma_max;
while (high-low) > epsilon
    mid = (low+high)/2;
    current = cfg;
    current.sigma = mid;
    current.y0 = [];
    current.early_stop = 'forward';
    [~,~,result] = solve_multiplex_species_noise(current);
    if result.A_final > threshold_A
        high = mid;
    else
        low = mid;
    end
end
sigma_f = (low+high)/2;

low = cfg.sigma_min;
high = cfg.sigma_max;
while (high-low) > epsilon
    mid = (low+high)/2;
    current = cfg;
    current.sigma = mid;
    current.y0 = cfg.y_seed;
    current.early_stop = 'backward';
    [~,~,result] = solve_multiplex_species_noise(current);
    if result.A_final > threshold_A
        high = mid;
    else
        low = mid;
    end
end
sigma_b = (low+high)/2;
thresholds = [sigma_f,sigma_b];
end
