clear;clc;
%% 数据初始化
% k=2;
% n=4;
% x=[1.09,2.41,3.321,4.55,5.767];
% V=Vander_G(n,k,x);
q=16;
t=5;
t_pre=t;%由于下面将t作为了循环变量，故这里再储存一个t
k=0;
N=2^t;
n=N-1;
omega = exp(-2j * pi / N);  % 复数单位根
x=zeros(n);
for i=0:n
    x(i+1)=omega^i;
end
V=Vander_G(n,k,x);

l=cell(n,1);%L子分解
u=cell(n,1);%U子分解
for i=1:n
    l{i}=zeros(n+1);
end
for i=1:n
    u{i}=zeros(n+1);
end

%% 1-带宽 LU分解实现
L=eye(n+1);
for m=1:n
    for i=1:n+1
        for j=1:n+1
            if i==j
                l{m}(i,j)=1;
            elseif i==j+1 && i>n-m+1
                l{m}(i,j)=(x(i)/x(j))^k;
                for t=0:m-n+i-3
                    l{m}(i,j)=l{m}(i,j)*(x(i)-x(i-1-t))/(x(i-1)-x(i-2-t));
                end
            end
        end
    end
    l{m}=Threshold(l{m});
end

U=eye(n+1);
for m=1:n-1
    for i=1:n+1
        for j=1:n+1
            if i==j && i<=n-m+1
                u{m}(i,j)=1;
            elseif i==j && i>n-m+1
                u{m}(i,j)=x(i)-x(n-m+1);
            elseif i==j-1 && i>=n-m+1
                u{m}(i,j)=x(m-n+i);
                for t=1:m-n+i-1
                    u{m}(i,j)=u{m}(i,j)*(x(i)-x(i-t))/(x(i+1)-x(i+1-t));
                end
            end

        end
    end
    u{m}=Threshold(u{m});
end

for i=1:n+1
    for j=1:n+1
        if i==j && i==1
            u{n}(i,j)=x(1)^k;
        elseif i==j && i~=1
            u{n}(i,j)=x(i)^k*(x(i)-x(1));
        elseif i==j-1
            u{n}(i,j)=x(i)^(k+1);
            for t=1:i-1
                u{n}(i,j)=u{n}(i,j)*(x(i)-x(i-t))/(x(i+1)-x(i+1-t));
            end
        end
    end
end

%% 计算复乘次数
complexity = 0;

% 计算 l{1} * l{2} * ... * l{n} 的复杂度
L=l{1};
for i = 2:n    
    complexity = complexity + computeComplexity(L,l{i});
    L=L*l{i};
end

% 计算 u{n} * u{n-1} * ... * u{1} 的复杂度
U=u{n};
for i = n-1:-1:1
    complexity = complexity + computeComplexity(U,u{i});
    U=U*u{i};
end

complexity=complexity+computeComplexity(L,U);

% 最终的复杂度为复杂度计数器乘以常数 q
C = complexity * q;
disp('L*U=');disp(L*U);
disp(['复杂度：' num2str(complexity)]);

%% 计算最小误差
options=optimoptions('fmincon', 'Algorithm', 'sqp', 'TolCon', 1e-6);
[best_beta,min_RMSE]=fmincon(@(x)norm_F(x,V,L*U),1,[],[],[],[],[],[],[],options);
disp(['best_beta=' num2str(best_beta)]);
disp(['min_RMSE=' num2str(min_RMSE)]);
disp(['C=' num2str(C)]);

% %保存当前t的复杂度、最小误差，不用时可注释
% filename = 'C_results.txt'; % 文件名
% fileID = fopen(filename, 'a'); % 打开文件以追加写入模式
% fprintf(fileID, 't=%d, RMSE=%.6f, C=%f\n', t_pre, min_RMSE, C);
% fclose(fileID);

%% 函数定义
%生成范德蒙德矩阵
function V=Vander_G(n,k,x)
    V=zeros(n+1);
    for i=1:n+1
        for j=1:n+1
            V(i,j)=x(i)^(k+j-1);
        end
    end
end

%优化目标函数
function f = norm_F(x,F_N,LU)
    beta=x;
    f=norm(beta*F_N-LU,'fro');
end

%阈值函数
function f=Threshold(M)
    for i = 1:numel(M)
        if abs(M(i))<1e-8
               M(i)=0;
        end
    end
    f=M;
end

