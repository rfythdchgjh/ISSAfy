clc
clear
close all

        pop=50; % 种群数
        M=3000; % 最大迭代次数
        dim = 20; % 可选 2, 10, 20
        
        %%  选择函数
        Function_name=4; % 函数名： 1 - 12
        [lb,ub,dim,fobj] = Get_Functions_cec2022(Function_name,dim);

        %% 调用算法
        Optimal_results = {}; % 保存Optimal results
        index = 1;

        % SCA
        tic
        [Destination_fitness, Destination_position, Curve_SCA] = SCA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "SCA";
        Optimal_results{2, index} = Curve_SCA;
        Optimal_results{3, index} = Destination_fitness;
        Optimal_results{4, index} = Destination_position;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % PSO
        tic
        [Best_score, Best_x, Curve_PSO] = PSO(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "PSO";
        Optimal_results{2, index} = Curve_PSO;
        Optimal_results{3, index} = Best_score;
        Optimal_results{4, index} = Best_x;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % HHO
        tic
        [Rabbit_Energy, Rabbit_Location, Curve_HHO] = HHO(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "HHO";
        Optimal_results{2, index} = Curve_HHO;
        Optimal_results{3, index} = Rabbit_Energy;
        Optimal_results{4, index} = Rabbit_Location;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % TTHHO
        tic
        [TRabbit_Energy, TRabbit_Location, CNVG] = TTHHO(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "TTHHO";
        Optimal_results{2, index} = CNVG;
        Optimal_results{3, index} = TRabbit_Energy;
        Optimal_results{4, index} = TRabbit_Location;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % OOA
        tic
        [Alpha_score, Alpha_pos, Curve_OOA] = OOA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "OOA";
        Optimal_results{2, index} = Curve_OOA;
        Optimal_results{3, index} = Alpha_score;
        Optimal_results{4, index} = Alpha_pos;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % SOOA
        tic
        [SBest_score, SBest_pos, SOOA_curve] = SOOA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "SOOA";
        Optimal_results{2, index} = SOOA_curve;
        Optimal_results{3, index} = SBest_score;
        Optimal_results{4, index} = SBest_pos;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % SSA
        tic
        [fMin, bestX, Curve_SSA] = SSA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "SSA";
        Optimal_results{2, index} = Curve_SSA;
        Optimal_results{3, index} = fMin;
        Optimal_results{4, index} = bestX;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % TSSA
        tic
        [TfMin, TbestX, Curve_TSSA] = TSSA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "TSSA";
        Optimal_results{2, index} = Curve_TSSA;
        Optimal_results{3, index} = TfMin;
        Optimal_results{4, index} = TbestX;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % ASFSSA
        tic
        [AfMin, AbestX, Curve_ASFSSA] = ASFSSA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "ASFSSA";
        Optimal_results{2, index} = Curve_ASFSSA;
        Optimal_results{3, index} = AfMin;
        Optimal_results{4, index} = AbestX;
        Optimal_results{5, index} = toc;
        index = index + 1;

        % ISSA
        tic
        [IfMin, IbestX, Curve_ISSA] = ISSA(pop, M, lb, ub, dim, fobj);
        Optimal_results{1, index} = "ISSA";
        Optimal_results{2, index} = Curve_ISSA;
        Optimal_results{3, index} = IfMin;
        Optimal_results{4, index} = IbestX;
        Optimal_results{5, index} = toc;
        index = index + 1;

          % 初始化一个单元数组来存储达到最优解的迭代次数
        iterations = cell(1, size(Optimal_results, 2));
        
        % 循环遍历每个算法的收敛曲线，找出达到最小值的迭代次数
        for i = 1:size(Optimal_results, 2)
            convergence_curve = Optimal_results{2, i};
            
            % 取整到最接近的整数（向下取整）
            integer_convergence_curve = floor(convergence_curve);
            
            % 找出最小值及其索引
            [min_value, min_idx] = min(integer_convergence_curve);
            
            % 存储达到最小值的迭代次数
            iterations{i} = min_idx + 1; % 迭代次数从1开始计数
        end
        
        % 创建并显示一个包含算法名称和达到最优解的迭代次数的表格
        T = table(Optimal_results(1,:), iterations, 'VariableNames', {'Algorithm', 'Iteration_at_Optimum'});
        
        disp(T);

figure
hold on

% 定义8种不同的颜色用于8种不同的算法
 colors = { '#b0d992', '#99b9e9', '#A89E96', '#E474C4', '#f7df87', '#05B9E2','#54beaa', '#eca680', '#fccccb', '#e3716e'};

% 定义标记类型
marker_types = { 'o', '^', '*','h', 'x', '<','s', '+', 'd', 'v'};
algorithm_labels = {'SCA', 'PSO', 'HHO','TTHHO', 'OOA','SOOA', 'SSA', 'TSSA', 'ASFSSA', 'ISSA'};

% 确保Optimal_results结构正确
assert(size(Optimal_results, 2) == length(colors), 'Number of algorithms does not match the number of colors and markers.');

for i = 1:size(Optimal_results, 2)
    algorithm = Optimal_results{1, i};
    convergence_curve = log10(Optimal_results{2, i});

    % 从颜色列表中获取颜色
    color = colors{i};
    % 从标记类型列表中获取标记
    marker = marker_types{i};

    % 确定标记间隔，例如每20个点标记一次
    marker_spacing = 200;
    marker_indices = 1:marker_spacing:length(convergence_curve);

    % 绘制曲线并指定颜色和标记
    plot(convergence_curve, 'Linewidth', 2, 'Color', color, 'Marker', marker, 'MarkerIndices', marker_indices, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerSize', 5);
end

title(['F' num2str(Function_name), '，Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Fitness' ]);
axis tight
%grid on % 可选：开启网格以增加可读性
box on

% % 显示图例，由于在plot中已经设置了颜色和标记，图例会自动反映这些属性
legend(algorithm_labels, 'Location', 'best');

hold off % 结束hold状态

display(['The best solution obtained by SCA is : ', num2str(Destination_position)]);
display(['The best optimal value of the objective funciton found by SCA is : ', num2str(Destination_fitness)]);

display(['The best solution obtained by PSO is : ', num2str(Best_x)]);
display(['The best optimal value of the objective funciton found by PSO is : ', num2str(Best_score)]);

display(['The best solution obtained by HHO is : ', num2str(Rabbit_Location)]);
display(['The best optimal value of the objective funciton found by HHO is : ', num2str(Rabbit_Energy)]);
display(['The best solution obtained by TTHHO is : ', num2str(TRabbit_Location)]);
display(['The best optimal value of the objective funciton found by TTHHO is : ', num2str(TRabbit_Energy)]);

display(['The best solution obtained by OOA is : ', num2str(Alpha_pos)]);
display(['The best optimal value of the objective funciton found by OOA is : ', num2str(Alpha_score)]);
display(['The best solution obtained by SOOA is : ', num2str(SBest_pos)]);
display(['The best optimal value of the objective funciton found by SOOA is : ', num2str(SBest_score)]);

display(['The best solution obtained by SSA is : ', num2str(bestX)]);
display(['The best optimal value of the objective funciton found by SSA is : ', num2str(fMin)]);
display(['The best solution obtained by TSSA is : ', num2str(TbestX)]);
display(['The best optimal value of the objective funciton found by TSSA is : ', num2str(TfMin)]);
display(['The best solution obtained by ASFSSA is : ', num2str(AbestX)]);
display(['The best optimal value of the objective funciton found by ASFSSA is : ', num2str(AfMin)]);
display(['The best solution obtained by ISSA is : ', num2str(IbestX)]);
display(['The best optimal value of the objective funciton found by ISSA is : ', num2str(IfMin)]);

