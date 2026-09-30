function [sigma_spec, lambda_star] = fig2c_topology_spectrum(layers, L_inter, alpha, ratio)
% 固定网络的确定性线性谱阈值，不含噪声。
% 稳态 (5,10): J=[10/3,-5;10,-4]。
% x=alpha*lambda 时 det=110/3+4*x-sigma*x*(10/3-x)，trace<0。
N = size(layers{1}, 1);
K = numel(layers);
if alpha <= 0 || ratio < 0 || ~isequal(size(L_inter), [K K])
    error('谱计算的 alpha、ratio 或层间矩阵尺寸无效。');
end
L_super = blkdiag(layers{:}) + ratio*kron(L_inter, speye(N));
lambda = sort(real(eig(full(L_super))));
if any(lambda < -1e-8)
    error('超拉普拉斯存在明显负特征值。');
end
lambda = lambda(lambda > 1e-10);
x = alpha*lambda;
valid = x < 10/3;
sigma_spec = NaN;
lambda_star = NaN;
if any(valid)
    thresholds = (110/3+4*x(valid))./(x(valid).*(10/3-x(valid)));
    [sigma_spec, i] = min(thresholds);
    valid_lambda = lambda(valid);
    lambda_star = valid_lambda(i);
end
end
