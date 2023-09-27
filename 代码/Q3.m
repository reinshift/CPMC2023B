clear; clc;
%% 数据初始化
t = 1;
N = 2^t;  % 矩阵阶数
F_N = dftmtx(N);  % 生成DFT矩阵

k = 4;  % 分解产生的矩阵数

q = 3;  % 元素范围
% 随机种子序列
seeds = zeros(1, 2*q+1);
for i = 1:q
    seeds(i) = 2^(i-1);
    seeds(q+i) = -2^(i-1);
end

% 遗传算法参数设置
populationSize = 200;         % 种群数量
chromosomeLength = k*N*N+1;  % 染色体长度
max_iter = 250;              % 最大迭代次数
mutationRate = 0.1;         % 变异率
crossoverRate = 0.85;         % 交叉率
tournamentSize =5;          % 锦标赛选择法中每轮锦标赛的参与个体数

%% 遗传算法实现
% 初始化种群
population = zeros(populationSize, chromosomeLength);
A=cell(1,k+1);
for i = 1:populationSize
    for j=1:k
        A{j}=Matrix_generate(N,seeds);
        population(i,(j-1)*N*N+1:j*N*N)=reshape(A{j},1,N*N);
    end
    population(i,end)=rand();
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
        tournamentIndices = randperm(populationSize, tournamentSize);%从所有个体中抽取出参与锦标赛的
        tournament = population(tournamentIndices, :);%出列
        tournamentFitness(i) = min(fitness(tournamentIndices));%让竞赛者竞争产生赢家
        %保留这个赢家的数据
        tournamentPopulation(i, :) = tournament(find(fitness(tournamentIndices) == tournamentFitness(i), 1, 'first'), :);
    end

    % 交叉
    for i = 1:2:populationSize-1
        r = rand();
        if r <= crossoverRate
            % 使用单点交叉法进行交叉，同时不破坏约束一
            if k<=2
                crossPoint = randi([1,N-1])*N;
            else
                crossPoint = min(randi([2,k])*N*N+1,randi([2,k-1])*N*N+1+randi([1,N-1])*N);
            end
            tournamentPopulation([i, i+1], crossPoint:end) = tournamentPopulation([i+1, i], crossPoint:end);
        end
    end

    % 变异
    for i = 1:populationSize
        for j = 1:chromosomeLength-1
            r = rand();%骰点
            if r < mutationRate
                r=rand();%骰点
                % 对染色体进行变异
                if r <=1 && r>0.1   %策略一
                    % 生成一个随机值替换当前基因
                    if population(i,j)~=0
                        population(i,j)=seeds(randperm(size(seeds,2),1))+1i*seeds(randperm(size(seeds,2),1));
                    end
                else                %策略二
                    % 敲除当前基因
                    population(i, j) = 0;
                end
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

%生成随机矩阵
function M=Matrix_generate(N,x)
    M=zeros(N);
    for i=1:N
        x1=zeros(1,N);
        x1(1)=x(randperm(size(x,2),1))+1i*x(randperm(size(x,2),1));
        x1(2)=x(randperm(size(x,2),1))+1i*x(randperm(size(x,2),1));
        x1=x1(randperm(size(x1,2),size(x1,2)));
        M(i,:)=x1;
    end
end
