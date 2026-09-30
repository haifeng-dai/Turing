function out = ba_lambda_star_spectrum(layers, alpha, ratio, options)
% BA 多层网络在前向确定性谱失稳处的主导谱点；层间是度为 2 的环。
% 不积分、不使用斑图种子。lambda_star 不是最大特征值。
% 大矩阵用 shift-invert eigs；收敛、残差、窗口边界及扩窗稳定性均须通过。
% 默认参数也可由调用者通过 options 覆盖，便于小规模完整谱交叉验证。
if nargin < 4, options = struct; end
defaults = struct('dense_limit', 600, 'dense_verify_limit', 1200, ...
    'eigs_counts', [16 32 64 128 256], 'eigs_tol', 1e-10, ...
    'eigs_maxit', 2000, 'residual_tol', 1e-7, 'comparison_tol', 1e-8);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(options, names{i}), options.(names{i}) = defaults.(names{i}); end
end
K = numel(layers);
assert(K >= 3 && alpha > 0 && isfinite(alpha) && ratio >= 0 && ...
    isfinite(ratio), '需要 K>=3、alpha>0、ratio>=0。');
N = size(layers{1}, 1);
for k = 1:K
    assert(isequal(size(layers{k}), [N N]), '各层矩阵尺寸必须相同。');
end
adj = sparse(1:K, [2:K 1], 1, K, K);
adj = adj + adj';
L_inter = 2*speye(K)-adj;
A = blkdiag(layers{:}) + ratio*kron(L_inter, speye(N));
n = size(A, 1);
lambda_opt = (-110+sqrt(16500))/(12*alpha);
lambda_limit = 10/(3*alpha);
out = struct('lambda_star', NaN, 'sigma_c_spec', NaN, ...
    'lambda_opt', lambda_opt, 'method', '', 'relative_residual', NaN, ...
    'neighborhood_bound', NaN, 'computed_modes', 0, ...
    'dense_verified', false, 'status', 'ok');

if n <= options.dense_limit
    [out.lambda_star, out.sigma_c_spec] = denseMinimum(A, alpha);
    out.method = 'dense';
    out.computed_modes = n;
    out.dense_verified = true;
else
    % 微移位避免恰好落在特征值上造成奇异分解；边界检验以实际移位为中心。
    shift = lambda_opt + 1e-6*max(1, lambda_opt);
    counts = unique(min(options.eigs_counts, n-2));
    counts = counts(counts >= 2);
    previous_sigma = NaN;
    accepted = false;
    last_problem = '候选窗口未通过边界和扩窗检验';
    % 私有随机流只为数值迭代提供起点，不改变拓扑，也不改变全局 rng。
    stream = RandStream('mt19937ar', 'Seed', 7919);
    v0 = randn(stream, n, 1);
    v0 = v0/norm(v0);
    for count = counts
        % 显式矩阵的对称性/实数类型由 eigs 自行识别；issym/isreal 仅用于函数句柄。
        opts = struct('tol', options.eigs_tol, ...
            'maxit', options.eigs_maxit, 'p', min(n, max(2*count+1, 20)), ...
            'v0', v0, 'disp', 0);
        try
            [V, D, flag] = eigs(A, count, shift, opts);
        catch err
            last_problem = err.message;
            previous_sigma = NaN;
            continue;
        end
        if flag ~= 0
            last_problem = sprintf('eigs 未完全收敛，flag=%d', flag);
            previous_sigma = NaN;
            continue;
        end
        lambda = real(diag(D));
        residuals = sqrt(sum(abs(A*V-V*D).^2, 1))' ./ ...
            (max(1, abs(lambda)).*sqrt(sum(abs(V).^2, 1))');
        residual = max(residuals);
        if any(~isfinite(lambda)) || residual > options.residual_tol
            last_problem = sprintf('特征对残差未通过：%.3g', residual);
            previous_sigma = NaN;
            continue;
        end
        [candidate, best] = discreteMinimum(lambda, alpha);
        if ~isfinite(best), previous_sigma = NaN; continue; end

        % sigma_c(lambda) 在 lambda_opt 左侧递减、右侧递增。
        % 对已经收敛的最近邻特征值集合，窗口之外的候选不能优于此下界。
        radius = max(abs(lambda-shift));
        left = shift-radius;
        right = shift+radius;
        outside_bound = Inf;
        if left > 0 && left < lambda_opt
            outside_bound = min(outside_bound, threshold(left, alpha));
        elseif left >= lambda_opt
            outside_bound = -Inf;
        end
        if right < lambda_limit && right > lambda_opt
            outside_bound = min(outside_bound, threshold(right, alpha));
        elseif right <= lambda_opt
            outside_bound = -Inf;
        end
        stable = isfinite(previous_sigma) && ...
            abs(best-previous_sigma) <= options.comparison_tol*max(1, best);
        bounded = outside_bound > best + options.comparison_tol*max(1, best);
        if stable && bounded
            out.lambda_star = candidate;
            out.sigma_c_spec = best;
            out.method = 'sparse_shift_invert';
            out.relative_residual = residual;
            out.neighborhood_bound = outside_bound;
            out.computed_modes = count;
            accepted = true;
            break;
        end
        previous_sigma = best;
    end
    if ~accepted
        if n <= options.dense_verify_limit
            [out.lambda_star, out.sigma_c_spec] = denseMinimum(A, alpha);
            out.method = 'dense_fallback';
            out.computed_modes = n;
            out.dense_verified = true;
        else
            error('ba_lambda_star_spectrum:UnverifiedSpectrum', ...
                '维度 %d 的稀疏谱未能验证：%s。可增大 eigs_counts 或 maxit。', ...
                n, last_problem);
        end
    elseif n <= options.dense_verify_limit
        [~, reference] = denseMinimum(A, alpha);
        if ~isfinite(reference) || ...
                abs(out.sigma_c_spec-reference) > options.comparison_tol*max(1, reference)
            error('ba_lambda_star_spectrum:CrossCheckFailed', ...
                '稀疏谱阈值 %.12g 与完整谱 %.12g 不一致。', out.sigma_c_spec, reference);
        end
        out.dense_verified = true;
    end
end
if ~isfinite(out.lambda_star), out.status = 'no_admissible_mode'; end
end

function [lambda_star, sigma] = denseMinimum(A, alpha)
lambda = eig(full(A));
if any(abs(imag(lambda)) > 1e-8)
    error('ba_lambda_star_spectrum:ComplexSpectrum', '对称网络出现明显复特征值。');
end
[lambda_star, sigma] = discreteMinimum(real(lambda), alpha);
end

function [lambda_star, sigma] = discreteMinimum(lambda, alpha)
if any(lambda < -1e-8)
    error('ba_lambda_star_spectrum:NegativeSpectrum', '拉普拉斯出现明显负特征值。');
end
lambda = sort(lambda(lambda > 1e-10 & alpha*lambda < 10/3));
lambda_star = NaN;
sigma = NaN;
if isempty(lambda), return; end
[sigma, idx] = min(threshold(lambda, alpha));
lambda_star = lambda(idx);
end

function sigma = threshold(lambda, alpha)
x = alpha*lambda;
sigma = (110/3+4*x)./(x.*(10/3-x));
end
