function cases = fig2c_topology_networks(config)
% 单独使用拓扑种子；生成结束后恢复调用方的随机状态。
previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng));
rng(config.topology_seed, 'twister');
N = config.N;
K = config.K;
if K < 3
    error('混合配置至少需要三层。');
end
types = {'WS', 'BA', 'ER'};
mixed = [types, types(randi(3, 1, K-3))];
mixed = mixed(randperm(K));
tags = {'pure_ws', 'pure_ba', 'mixed'};
labels = {'Pure WS', 'Pure BA', 'WS-BA-ER mixed'};
layer_types = {repmat({'WS'}, 1, K), repmat({'BA'}, 1, K), mixed};
cases = struct('tag', {}, 'label', {}, 'layer_types', {}, ...
    'L', {}, 'mean_degrees', {});
for c = 1:3
    cases(c).tag = tags{c};
    cases(c).label = labels{c};
    cases(c).layer_types = layer_types{c};
    cases(c).L = cell(1, K);
    cases(c).mean_degrees = zeros(1, K);
    for k = 1:K
        switch layer_types{c}{k}
            case 'WS'
                L = gen_ws(N, config.ws_degree, config.ws_rewire);
            case 'BA'
                L = gen_ba(N, config.ba_m);
            case 'ER'
                L = gen_er(N, config.er_p);
        end
        cases(c).L{k} = sparse(L);
        cases(c).mean_degrees(k) = mean(diag(L));
    end
end
end
