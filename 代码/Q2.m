clear; clc;
%% 数据初始化
t = 5;
N = 2^t;  % 矩阵阶数
F_N = dftmtx(N);  % 生成DFT矩阵

k = 3;  % 分解产生的矩阵数

q = 3;  % 元素范围
% 随机种子序列
seeds = zeros(1, 2*q+1);
for i = 1:q
    seeds(i) = 2^(i-1);
    seeds(q+i) = -2^(i-1);
end
rand_seeds = zeros(2*q+1);
for i = 1:length(seeds)
    for j = 1:length(seeds)
        rand_seeds(i, j) = seeds(i) + 1i * seeds(j);
    end
end
rand_seeds = reshape(rand_seeds, 1, (2*q+1)^2);

% 遗传算法参数设置
populationSize = 150;         % 种群数量
chromosomeLength = k*N*N+1;  % 染色体长度
max_iter = 500;              % 最大迭代次数
mutationRate = 0.03;         % 变异率
crossoverRate = 0.7;         % 交叉率
tournamentSize =40;          % 锦标赛选择法中每轮锦标赛的参与个体数

%% 遗传算法实现
% 初始化种群
population = zeros(populationSize, k*N*N+1);
for i = 1:populationSize
    for j = 1:k*N*N
        perm = rand_seeds(randperm(numel(rand_seeds)));
        population(i, j) = perm(1);
    end
    population(i,k*N*N+1)=rand();
end

% 存储每次迭代的最优适应度值
bestFitnessHistory = zeros(max_iter, 1);

for iter = 1:max_iter
    % 计算适应度
    fitness = zeros(populationSize, 1);
    for i = 1:populationSize
        fitness(i) = norm_F(population(i, :), F_N, k, N);
    end

    % 记录当前迭代的最优适应度值以及最优解
    [bestFitness,index] = min(fitness);
    bestSolution = population(index,:);
    bestFitnessHistory(iter) = bestFitness;
    %fprintf('Iteration: %d, Best Fitness: %.4f\n', iter, bestFitness);

    % 选择
    % 锦标赛选择法
    tournamentPopulation = zeros(populationSize, chromosomeLength);
    for i = 1:populationSize
        tournamentIndices = randperm(populationSize, tournamentSize);
        tournament = population(tournamentIndices, :);
        tournamentFitness(i) = min(fitness(tournamentIndices));
        tournamentPopulation(i, :) = tournament(find(fitness(tournamentIndices) == tournamentFitness(i), 1, 'first'), :);
    end

    % 交叉
    for i = 1:2:populationSize-1
        r = rand();
        if r < crossoverRate
            % 使用单点交叉法进行交叉
            crossPoint = randi([1, chromosomeLength-1]);
            tournamentPopulation([i, i+1], crossPoint:end-1) = tournamentPopulation([i+1, i], crossPoint:end-1);
        end
    end

    % 变异
    for i = 1:populationSize
        for j = 1:chromosomeLength-1
            r = rand();
            if r < mutationRate
                % 小概率进行变异
                perm = rand_seeds(randperm(numel(rand_seeds)));
                tournamentPopulation(i, j) = perm(1);
            end
        end
    end

    % 更新种群
    population = tournamentPopulation;
end

%% 绘制适应度随迭代次数变化的图表
figure;
plot(1:max_iter, bestFitnessHistory, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Fitness');
title('Fitness Evolution');
grid on;

%% 解码染色体
A=cell(1,k+1);
for i=1:k
    A{i}=reshape(bestSolution(N^2*(i-1)+1:N^2*i),N,N);
end
A{k+1}=bestSolution(end);

%% 计算硬件复杂度
Ak=A{1};l=0;
for i=2:k
    l=l+computeComplexity(Ak,A{i});
    Ak=Ak*A{i};
end
if k==1 C=0;
else
    C= q*l;
end
disp(['硬件复杂度为：' num2str(C)]);

%% 适应度函数和其他函数
%适应度-目标函数
function f = norm_F(x,F_N,k,N)
    A=cell(1,k);
    for i=1:k
        A{i}=reshape(x(N^2*(i-1)+1:N^2*i),N,N);
    end
    Ak=eye(N);
    for i=1:k
        Ak=Ak*A{i};
    end
    f=1/N*norm(F_N-x(end)*Ak,'fro');
end
