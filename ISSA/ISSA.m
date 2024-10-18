
function [IfMin , IbestX,Curve_ISSA ] = ISSA(pop, M,c,d,dim,fobj  )


P_percent = 0.2;    % The population size of producers accounts for "P_percent" percent of the total population size
        if mod(M, 10) == 0

            P_percent = 0.2 + (0.5 - 0.2) * rand(); % 在20%到50%之间随机取值

        end
pNum = round( pop *  P_percent );    % The population size of the producers
        % initialDisturbanceFactor = 1; % 初始扰动因子
        % 
        % decayRate = 0.99; % 衰减率
        

lb= c.*ones( 1,dim );    % Lower limit/bounds/     a vector
ub= d.*ones( 1,dim );    % Upper limit/bounds/     a vector
% Initialization
% ★★改进1：三角形拓扑结构种群初始化★★
X = initializationT(pop, dim, ub, lb);

% ★★改进2：反向学习种群初始化★★
% 注意：这里我们将三角形拓扑结构种群 x_triangle 作为输入传递给反向学习函数
Sol_forward = initialization_for(X, pop, dim, ub, lb);
Sol_backward = initialization_back(Sol_forward, pop, dim, ub, lb);

% 将两个反向学习种群合并为一个大种群
x_all = [Sol_forward; Sol_backward];

% 计算大种群中每个个体的适应度值
 % Calculate fitness values for all individuals in the combined population
    x_all_fitness = zeros(size(x_all, 1), 1);
    for i = 1 : size(x_all, 1)
        x_all_fitness(i) = fobj(x_all(i, :));
    end
    
    % Select the top 'pop' individuals from the combined population based on fitness
    [~, sorted_indexes] = sort(x_all_fitness, 'ascend');
    x = x_all(sorted_indexes(1:pop), :);
    
    % Initialize fitness values for the selected population
    fit = zeros(pop, 1);
    for i = 1 : pop
        fit(i) = fobj(x(i, :));
    end

% 更新全局最优适应度值 pFit 和对应的位置 pX
pFit = fit;
pX = x;
[IfMin, bestI] = min(pFit);
IbestX = x(bestI, :);
%% %发现者位置更新
  % Start updating the solutions.

for t = 1 : M    
  
      
  [ ans, sortIndex ] = sort( pFit );% Sort.
     
  [fmax,B]=max( pFit );
   worse= x(B,:);  

   S=(2*sqrt(rand)-1)*(1+rand)/rand;
   r2=rand(1);
if(r2<0.8)
 r1=rand(1);
 
    for i = 1 : pNum                                                   % Equation (3)
        %% 动态分级
         
          % 确定优势种群、中等种群和劣势种群的大小
          m = round(pop * P_percent); % 优势种群的大小
          n = round((pop - m) / 2);   % 中等种群的大小
          r = pop - m - n;            % 劣势种群的大小
          % 根据排序后的索引，将种群分为三个部分
          advantageIndex = sortIndex(1:m);
          middleIndex = sortIndex(m+1:m+n);
          disadvantageIndex = sortIndex(m+n+1:end);
 for i = 1:m
    % 计算两种权重
    weight1 = exp(1 - (M + t) / (M - t));
    weight2 = ((M - t + 1) / M)^t;

    % 根据第一种权重计算新的解
    newX1 = pX(sortIndex(i), :) .* exp(-(i) / (r1 * M)) .* weight1;
    newX1 = Bounds(newX1, lb, ub);
    newFit1 = fobj(newX1);

    % 根据第二种权重计算新的解
    newX2 = pX(sortIndex(i), :) .* exp(-(i) / (r1 * M)) .* weight2;
    newX2 = Bounds(newX2, lb, ub);
    newFit2 = fobj(newX2);

    % 选择适应度较小的解,贪婪规则
    if newFit1 <= newFit2
        x(advantageIndex(i), :) = newX1;
        fit(advantageIndex(i)) = newFit1;
    elseif newFit1 > newFit2
        x(advantageIndex(i), :) = newX2;
        fit(advantageIndex(i)) = newFit2;
    end


        % % ★★改进2：随机惯性权重策略★★
        % w = ((M - t + 1) / M)^t;
        % x(advantageIndex(i), :) = w* pX( advantageIndex(i), : )*exp(-(i)/(r1*M));
        x(advantageIndex(i), :) = Bounds(x(advantageIndex(i), :), lb, ub);
        fit(advantageIndex(i)) = fobj( x( advantageIndex(i), : ) );  
 end
          for i = 1:n
             % 中等种群的更新策略,这里可以使用常规的更新策略
             % ★★改进3：原来策略★★
             x(middleIndex(i), :) = pX( middleIndex(i), : )*exp(-(i)/(r1*M));
             x(middleIndex(i), :) = Bounds(x(middleIndex(i), :), lb, ub);
             fit(middleIndex(i)) = fobj( x( middleIndex(i), : ) );  
          end
          for i = 1:r
             % 劣势种群的更新策略
             % ★★改进4：%%蝴蝶优化策略
             power_exponent=0.1; % 幂指数
             Fnew=fobj(x( disadvantageIndex(i), : ));
             sensory_modality=0.01; % 感觉因子
             FP=(sensory_modality*(Fnew^power_exponent)); % 每只蝴蝶的香味
             dis = rand * rand * IbestX - pX( disadvantageIndex(i), : );
             x( disadvantageIndex(i), : )=pX( disadvantageIndex(i), : )+dis*FP;
             x(disadvantageIndex(i), :) = Bounds(x(disadvantageIndex(i), :), lb, ub);
             fit(disadvantageIndex(i)) = fobj( x( disadvantageIndex(i), : ) );  
          end
    end
          else
        % ★★改进3：引入柯西变异★★(随机权重优化)
        scale = 0.1; % 设置柯西分布的缩放参数，你可以根据需要调整这个值
        for i = 1 : pNum
            % % 生成柯西分布的随机数
            % cauchyRand = tan(pi * (rand() - 0.5));
            % x(sortIndex(i), :) = pX(sortIndex(i), :) + scale * cauchyRand;

            ori_value = rand(1,dim);
            cauchy_value = tan((ori_value-0.5)*pi);
            x(sortIndex(i), :) = pX(sortIndex(i), :) + pX(sortIndex(i), :).*cauchy_value;

            x(sortIndex(i), :) = Bounds(x(sortIndex(i), :), lb, ub);
            fit(sortIndex(i)) = fobj(x(sortIndex(i), :));
        end
        % for i = 1 : pNum
        % 
        %     x( sortIndex( i ), : ) = pX( sortIndex( i ), : )+randn(1)*ones(1,dim);
        %     x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
        %     fit( sortIndex( i ) ) = fobj( x( sortIndex( i ), : ) );
        % 
        % end

end
    % 更新完优势、中等、劣势种群后，重新计算所有个体的适应度
    for i = 1 : pop
        fit(i) = fobj(x(i, :)); % 重新计算适应度
    end


 [ fMMin, bestII ] = min( fit );      
  bestXX = x( bestII, : );            
        
   % k=2*rand()-1;
   %      c=exp(5*cos(pi*(1-t/M)));
     for i = ( pNum + 1 ) : pop                     % Equation (4)
                    
              A=floor(rand(1,dim)*2)*2-1;
              if( i>(pop/2))
                % 步骤1：变异
                F = 0.5;              % 缩放因子
                n = randperm(pop, 5);
                x( sortIndex( i ), : ) = pX(sortIndex(n(1)), :) + F * (pX(sortIndex(n(2)), :) - pX(sortIndex(n(3)), :)) + F * (pX(sortIndex(n(4)), :) - pX(sortIndex(n(5)), :));
              else
                  x( sortIndex( i ), : )=bestXX+(abs(( pX( sortIndex( i ), : )-bestXX)))*(A'*(A*A')^(-1))*ones(1,dim);
              end
              x( sortIndex( i ), : ) = Bounds( x( sortIndex( i ), : ), lb, ub );
              fit( sortIndex( i ) ) = fobj( x( sortIndex( i ), : ) );
              
     end
    c=randperm(numel(sortIndex));
    b=sortIndex(c(1:20));

for j =  1  : length(b)      % Equation (5)
  % 改进2.1：★★横向交叉策略★★
  if (mod(j,2) == 1)
      for i = 1 : dim
          r = rand();
          c = 2*rand()-1;
          x_criss( sortIndex( b(j) ), i ) = r * pX( sortIndex( b(j) ), i )+(1-r)*pX( sortIndex( b(j+1) ), i )+c*(pX( sortIndex( b(j) ), i )-pX( sortIndex( b(j+1) ), i ));
      end
  else
      for i = 1 : dim
          r = rand();
          c = 2*rand()-1;
          x_criss( sortIndex( b(j) ), i ) = r * pX( sortIndex( b(j-1) ), i )+(1-r)*pX( sortIndex( b(j) ), i )+c*(pX( sortIndex( b(j-1) ), i )-pX( sortIndex( b(j) ), i ));
      end
  end
  
  if fobj( x_criss( sortIndex( b(j) ), : ) ) < fobj( pX( sortIndex( b(j) ), : ) )
      x( sortIndex( b(j) ), : ) = x_criss( sortIndex( b(j) ), : );
  else
      x( sortIndex( b(j) ), : ) = pX( sortIndex( b(j) ), : );
  end
  
  x( sortIndex(b(j) ), : ) = Bounds( x( sortIndex(b(j) ), : ), lb, ub );
  fit( sortIndex( b(j) ) ) = fobj( x( sortIndex( b(j) ), : ) );
  
  % 改进2.2：★★纵向交叉策略★★
  for i = 1 : dim
      d1 = randperm(dim,1);
      d2 = randperm(dim,1);
      r1 = rand();
      x_cross(sortIndex( b(j) ), i) = r1*x( sortIndex(b(j) ), d1 ) + (1-r1)*x( sortIndex(b(j) ), d2 );
  end
  
  if fobj( x_cross( sortIndex( b(j) ), : ) ) < fobj( x( sortIndex( b(j) ), : ) )
      x( sortIndex( b(j) ), : ) = x_cross( sortIndex( b(j) ), : );
  else
      x( sortIndex( b(j) ), : ) = x( sortIndex( b(j) ), : );
  end
  
  x( sortIndex(b(j) ), : ) = Bounds( x( sortIndex(b(j) ), : ), lb, ub );
  fit( sortIndex( b(j) ) ) = fobj( x( sortIndex( b(j) ), : ) );

end
    
    for i = 1 : pop
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end

        if( pFit( i ) < IfMin )
            IfMin= pFit( i );
            IbestX = pX( i, : );


        end
    end

    Curve_ISSA(t)=IfMin;

end


% Application of simple limits/bounds
function s = Bounds( s, Lb, Ub)
% Apply the lower bound vector
temp = s;
I = temp < Lb;
temp(I) = Lb(I);

% Apply the upper bound vector
J = temp > Ub;
temp(J) = Ub(J);
% Update this new move
s = temp;

%---------------------------------------------------------------------------------------------------------------------------
