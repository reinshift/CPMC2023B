%DFT函数定义
function F = dftmtx(N)
    F = zeros(N, N);
    omega = exp(-2j * pi / N);  % 复数单位根

    for n = 0:N-1
        for k = 0:N-1
            F(n+1, k+1) = omega^(n*k);
        end
    end

    F = 1 / sqrt(N) * F;  % 归一化
end
