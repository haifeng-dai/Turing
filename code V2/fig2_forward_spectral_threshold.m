function [sigma_spec, lambda_star] = fig2_forward_spectral_threshold( ...
    topo_type, N, K, p_value, alpha_values, ratio_values)
% 确定性均匀态 (5,10) 的线性谱失稳阈值；不包含乘性噪声。
% J=[10/3,-5;10,-4]。令 x=alpha*lambda，trace 始终为负，
% det=110/3+4*x-sigma*x*(10/3-x)，因此逐模态求 det=0 后取最小值。
if isscalar(alpha_values)
    alpha_values = repmat(alpha_values, size(ratio_values));
elseif isscalar(ratio_values)
    ratio_values = repmat(ratio_values, size(alpha_values));
end
alpha_values = alpha_values(:)';
ratio_values = ratio_values(:)';
if numel(alpha_values) ~= numel(ratio_values) || ...
        any(~isfinite(alpha_values) | alpha_values <= 0) || ...
        any(~isfinite(ratio_values) | ratio_values < 0)
    error('谱计算需要等长的正 alpha 和非负 ratio 参数。');
end
topology_file = fullfile(fileparts(mfilename('fullpath')), 'results', ...
    'topology', upper(topo_type), sprintf('N%d_p%.3f.mat', N, p_value));
S = load(topology_file, 'nets');
if ~iscell(S.nets) || numel(S.nets) < K
    error('网络文件没有足够的层：%s', topology_file);
end
layers = S.nets(1:K);
for k = 1:K
    L = layers{k};
    if ~isequal(size(L), [N N]) || norm(L-L', 'fro') > 1e-8 || ...
            max(abs(full(sum(L, 2)))) > 1e-8
        error('第 %d 层不是所需的对称 N×N 拉普拉斯矩阵。', k);
    end
end
L_intra = blkdiag(layers{:});
L_inter = kron(K*eye(K)-ones(K), speye(N));
sigma_spec = NaN(size(alpha_values));
lambda_star = NaN(size(alpha_values));
ratios = unique(ratio_values);
for r = ratios
    lambda = sort(real(eig(full(L_intra+r*L_inter))));
    if any(lambda < -1e-8)
        error('ratio=%.6g 的网络谱存在明显负特征值。', r);
    end
    lambda = lambda(lambda > 1e-10);
    for i = find(ratio_values == r)
        x = alpha_values(i)*lambda;
        valid = x > 0 & x < 10/3;
        if ~any(valid)
            warning('alpha=%.6g、ratio=%.6g 没有可失稳的非零模态。', ...
                alpha_values(i), r);
            continue;
        end
        thresholds = (110/3+4*x(valid))./(x(valid).*(10/3-x(valid)));
        [sigma_spec(i), idx] = min(thresholds);
        valid_lambda = lambda(valid);
        lambda_star(i) = valid_lambda(idx);
    end
end
end
