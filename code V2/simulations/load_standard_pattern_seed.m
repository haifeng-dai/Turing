function [y_seed, seed_file] = load_standard_pattern_seed(cfg, allow_other_N)
% LOAD_STANDARD_PATTERN_SEED 统一读取预生成 ER 强斑图初值。
% 同一 N、K 不因拓扑、alpha、beta、噪声改变而重新训练。
% N=200 缺失时始终报错；只有显式允许且 N~=200 才可生成该尺寸的种子。
% K 改变必须有尺寸匹配的标准文件，N=200 时不会自动生成。
if nargin < 2, allow_other_N = false; end
validateattributes(cfg.N, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.K, {'numeric'}, {'scalar','integer','positive'});
project_dir = fileparts(fileparts(mfilename('fullpath')));
results_dir = fullfile(project_dir, 'results');
seed_file = fullfile(results_dir, sprintf( ...
    'evolution_er_N%d_K%d_p0.030_a0.050_b0.005_s100.0_n0.00_fwd_results.mat', cfg.N, cfg.K));
if ~exist(seed_file, 'file')
    if ~allow_other_N || cfg.N == 200
        error('load_standard_pattern_seed:MissingSeed', ...
            ['缺少预生成标准斑图种子：%s\nN=200 不自动训练；K 改变也必须读取匹配尺寸的文件。' ...
            '\n请先准备该 N、K 的标准 ER 种子，不要用不同尺寸的状态代替。'], seed_file);
    end
    % 不沿用目标模型/拓扑参数：新 N 仍使用项目标准 ER、alpha=.05、beta=.005。
    net_file = fullfile(results_dir, 'topology', 'ER', sprintf('N%d_p0.030.mat',cfg.N));
    if ~exist(net_file,'file')
        error('load_standard_pattern_seed:MissingNetwork','不同 N 的标准 ER 网络文件缺失：%s',net_file);
    end
    data = load(net_file,'nets');
    if ~isfield(data,'nets') || numel(data.nets)<cfg.K
        error('load_standard_pattern_seed:InvalidNetwork','标准 ER 网络不足 %d 层：%s',cfg.K,net_file);
    end
    fprintf('[SEED] N=%d 不同于基准 N=200，仅为新尺寸生成一次标准种子：%s\n',cfg.N,seed_file);
    saved_rng = rng;
    restore_rng = onCleanup(@()rng(saved_rng));
    rng(0,'twister');
    seed_cfg = struct('N',cfg.N,'K',cfg.K,'alpha',0.05,'beta',0.005, ...
        'sigma',100,'noise',0,'noise_seed',0,'T_END',500,'dt',0.001,'steps',2, ...
        'init_perturb',0.1,'detect_convergence',true,'y0',[], ...
        'L_intra',{data.nets(1:cfg.K)},'L_inter',sparse(cfg.K*eye(cfg.K)-ones(cfg.K)));
    [t,Y,generated_cfg] = solve_multiplex(seed_cfg);
    validateSeed(Y,cfg.N,cfg.K,seed_file);
    saved = struct('t',t,'Y',Y,'cfg',generated_cfg);
    temporary_file = [seed_file '.tmp.mat'];
    save(temporary_file,'-struct','saved','-v7');
    [ok,msg] = movefile(temporary_file,seed_file,'f');
    if ~ok, error('load_standard_pattern_seed:SaveFailed','标准种子保存失败：%s',msg); end
    clear restore_rng;
end
data = load(seed_file,'Y');
if ~isfield(data,'Y')
    error('load_standard_pattern_seed:InvalidSeed','标准种子文件缺少 Y：%s',seed_file);
end
validateSeed(data.Y,cfg.N,cfg.K,seed_file);
y_seed = data.Y(end,:)';
fprintf('[SEED] 读取统一标准种子：%s\n',seed_file);
end

function validateSeed(Y,N,K,path)
if ~isnumeric(Y) || ~ismatrix(Y) || isempty(Y) || ~isreal(Y) || size(Y,2)~=2*N*K || ...
        any(~isfinite(Y(end,:)))
    error('load_standard_pattern_seed:InvalidSeed','标准种子必须是有限实数状态，长度 2*N*K=%d：%s',2*N*K,path);
end
NK = N*K;
y = Y(end,:);
A = sqrt(sum((y(1:NK)-5).^2+(y(NK+1:end)-10).^2)/NK);
if A<=0.05
    error('load_standard_pattern_seed:InvalidSeed','标准种子未形成明显斑图 (A=%.6g)：%s',A,path);
end
end
