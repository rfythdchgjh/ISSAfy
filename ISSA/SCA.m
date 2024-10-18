function [Destination_fitness, Destination_position, Curve_SCA] = SCA(pop, M, lb, ub, dim, fobj)

    % 确保 pop 和 dim 是标量
    assert(isnumeric(pop) && isscalar(pop), 'pop 必须是一个标量');
    assert(isnumeric(dim) && isscalar(dim), 'dim 必须是一个标量');

    %% 初始化
    X = zeros(pop, dim);  % 预先分配空间
    for i = 1:dim
        X(:, i) = lb(i) + rand(pop, 1) .* (ub(i) - lb(i));  % 生成 [lb, ub] 之间的随机数
    end

    % 计算初始种群的适应度
    Objective_values = arrayfun(@(i) fobj(X(i, :)), 1:pop);

    % 初始化目标位置和适应度
    [~, best_idx] = min(Objective_values);
    Destination_position = X(best_idx, :);
    Destination_fitness = Objective_values(best_idx);

    % 初始化记录每一代的最佳适应度
    Curve_SCA = zeros(1, M);
    Curve_SCA(1) = Destination_fitness;

    %% 迭代
    t = 2;  % 从第二次迭代开始
    while t <= M
        % Eq. (3.4)
        a = 2;
        r1 = a - t * ((a) / M);  % r1 线性递减从 a 到 0

        % 位置更新
        for i = 1:pop  % 第 i 个解
            for j = 1:dim  % 第 j 维
                % 更新 r2, r3, 和 r4 用于 Eq. (3.3)
                r2 = (2 * pi) * rand();
                r3 = 2 * rand();
                r4 = rand();

                % Eq. (3.3)
                if r4 < 0.5
                    % Eq. (3.1)
                    X(i, j) = X(i, j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i, j)));
                else
                    % Eq. (3.2)
                    X(i, j) = X(i, j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i, j)));
                end
            end
        end

        % 越界规范
        for i = 1:pop
            % 如果解超出搜索空间，则将其拉回
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = (X(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % 计算目标值
            Objective_values(i) = fobj(X(i, :));

            % 更新目标位置如果找到了更好的解
            if Objective_values(i) < Destination_fitness
                Destination_position = X(i, :);
                Destination_fitness = Objective_values(i);
            end
        end

        Curve_SCA(t) = Destination_fitness;

        t = t + 1;
    end
end