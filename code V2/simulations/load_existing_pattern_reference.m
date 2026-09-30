function [state, A_reference, source_file] = load_existing_pattern_reference(results_dir, topology_type, topology_parameter, cfg, A_cut)
% 吸引概率插值的目标斑图：仅读取，不训练、不从高 sigma 做延拓。
% 文件名精度有限，因此目标参数以 MAT 中保存的 cfg 为准。
files = dir(fullfile(results_dir, sprintf('evolution_%s_N%d_K%d_p*_results.mat', ...
    lower(topology_type),cfg.N,cfg.K)));
pattern = ['^evolution_' lower(topology_type) '_N\d+_K\d+_p([0-9.]+)' ...
    '_a[0-9.]+_b[0-9.]+_s[0-9.]+_n[0-9.]+_(fwd|bwd)_results\.mat$'];
% 同参数都存在时优先读取反向斑图，其次读取已经形成斑图的正向终态。
names = {files.name};
is_bwd = contains(names,'_bwd_results.mat');
files = files([find(is_bwd),find(~is_bwd)]);
for i = 1:numel(files)
    token = regexp(files(i).name,pattern,'tokens','once');
    if isempty(token) || abs(str2double(token{1})-topology_parameter)>5.1e-4
        continue;
    end
    path = fullfile(files(i).folder,files(i).name);
    saved = load(path,'Y','cfg');
    if ~isfield(saved,'Y') || ~isnumeric(saved.Y) || isempty(saved.Y) || ...
            ~isreal(saved.Y) || size(saved.Y,2)~=2*cfg.N*cfg.K || ...
            any(~isfinite(saved.Y(end,:))) || ~isfield(saved,'cfg')
        continue;
    end
    fields = {'N','K','alpha','beta','sigma','noise'};
    match = true;
    for j = 1:numel(fields)
        key = fields{j};
        if ~isfield(saved.cfg,key) || ~isnumeric(saved.cfg.(key)) || ...
                ~isscalar(saved.cfg.(key)) || ~isfinite(saved.cfg.(key)) || ...
                abs(saved.cfg.(key)-cfg.(key))>1e-10*max(1,abs(cfg.(key)))
            match = false;
            break;
        end
    end
    if ~match || ~isfield(saved.cfg,'L_intra') || ~iscell(saved.cfg.L_intra) || ...
            numel(saved.cfg.L_intra)<cfg.K || ~isfield(saved.cfg,'L_inter') || ...
            ~isequal(sparse(saved.cfg.L_inter),sparse(cfg.L_inter))
        continue;
    end
    for k = 1:cfg.K
        if ~isequal(sparse(saved.cfg.L_intra{k}),sparse(cfg.L_intra{k}))
            match = false;
            break;
        end
    end
    if ~match, continue; end
    candidate = saved.Y(end,:)';
    NK = cfg.N*cfg.K;
    endpoint_A = sqrt(sum((candidate(1:NK)-5).^2+(candidate(NK+1:end)-10).^2)/NK);
    candidate_A = endpoint_A;
    if isfield(saved.cfg,'A_final')
        candidate_A = saved.cfg.A_final;
    end
    if ~isnumeric(candidate_A) || ~isscalar(candidate_A) || ...
            ~isfinite(candidate_A) || candidate_A<=A_cut || endpoint_A<=A_cut
        continue;
    end
    state = candidate;
    A_reference = candidate_A;
    source_file = path;
    fprintf('[REFERENCE READ] 读取已有目标斑图：%s\n',source_file);
    return;
end
error('load_existing_pattern_reference:MissingReference', ...
    ['严格读取模式：没有匹配的已有目标斑图；不会生成或自动延拓。\n' ...
    '需要 N=%d, K=%d, alpha=%.10g, beta=%.10g, sigma=%.10g, noise=%.10g，' ...
    '且网络矩阵一致、斑图有效。\n请先准备并保存对应 evolution_* 文件（MAT 中需包含 Y、cfg）。'], ...
    cfg.N,cfg.K,cfg.alpha,cfg.beta,cfg.sigma,cfg.noise);
end
