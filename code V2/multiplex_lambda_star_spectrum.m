function out = multiplex_lambda_star_spectrum(layers, L_inter, alpha, ratio, options)
% MULTIPLEX_LAMBDA_STAR_SPECTRUM 计算指定层间拓扑下的主导谱点 lambda_*。
% lambda_* 是使前向确定性谱阈值最小的正非均匀超拉普拉斯特征值，
% 不是最大特征值。layers 和 L_inter 均为拉普拉斯矩阵。
if nargin < 5, options = struct; end
defaults = struct('dense_limit', 600, 'dense_verify_limit', 1200, ...
    'eigs_counts', [16 32 64 128 256], 'eigs_tol', 1e-10, ...
    'eigs_maxit', 2000, 'residual_tol', 1e-7, 'comparison_tol', 1e-8);
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(options, names{i}), options.(names{i}) = defaults.(names{i}); end
end

K = numel(layers);
assert(K >= 2 && alpha > 0 && isfinite(alpha) && ratio >= 0 && ...
    isfinite(ratio) && isequal(size(L_inter), [K K]), ...
    '需要 K>=2、alpha>0、ratio>=0，且层间拉普拉斯尺寸为 K×K。');
N = size(layers{1}, 1);
for k = 1:K
    L = sparse(layers{k});
    assert(isequal(size(L), [N N]) && isreal(L) && ...
        all(isfinite(nonzeros(L))), '第 %d 层拉普拉斯尺寸或数值无效。', k);
    assert(norm(L-L', 1) < 1e-8 && max(abs(sum(L, 2))) < 1e-8, ...
        '第 %d 层矩阵不是对称零行和拉普拉斯。', k);
    layers{k} = L;
end
L_inter = sparse(L_inter);
assert(isreal(L_inter) && norm(L_inter-L_inter', 1) < 1e-10 && ...
    max(abs(sum(L_inter, 2))) < 1e-10, ...
    '层间矩阵必须是对称零行和拉普拉斯。');

A = blkdiag(layers{:}) + ratio*kron(L_inter, speye(N));
A = sparse(A);
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
    shift = lambda_opt + 1e-6*max(1, lambda_opt);
    counts = unique(min(options.eigs_counts, n-2));
    counts = counts(counts >= 2);
    previous_sigma = NaN;
    accepted = false;
    last_problem = '候选谱窗口未通过边界和扩窗检验';
    stream = RandStream('mt19937ar', 'Seed', 7919);
    v0 = randn(stream, n, 1);
    v0 = v0/norm(v0);

    for count = counts
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
        if ~isfinite(best)
            previous_sigma = NaN;
            continue;
        end

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
            error('multiplex_lambda_star_spectrum:UnverifiedSpectrum', ...
                '维度 %d 的稀疏谱未能验证：%s。可增大 eigs_counts 或 maxit。', ...
                n, last_problem);
        end
    elseif n <= options.dense_verify_limit
        [~, reference] = denseMinimum(A, alpha);
        if ~isfinite(reference) || ...
                abs(out.sigma_c_spec-reference) > ...
                options.comparison_tol*max(1, reference)
            error('multiplex_lambda_star_spectrum:CrossCheckFailed', ...
                '稀疏谱阈值 %.12g 与完整谱 %.12g 不一致。', ...
                out.sigma_c_spec, reference);
        end
        out.dense_verified = true;
    end
end
if ~isfinite(out.lambda_star), out.status = 'no_admissible_mode'; end
end

function [lambda_star, sigma] = denseMinimum(A, alpha)
lambda = eig(full(A));
if any(abs(imag(lambda)) > 1e-8)
    error('multiplex_lambda_star_spectrum:ComplexSpectrum', ...
        '对称超拉普拉斯出现明显复特征值。');
end
[lambda_star, sigma] = discreteMinimum(real(lambda), alpha);
end

function [lambda_star, sigma] = discreteMinimum(lambda, alpha)
if any(lambda < -1e-8)
    error('multiplex_lambda_star_spectrum:NegativeSpectrum', ...
        '超拉普拉斯出现明显负特征值。');
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
