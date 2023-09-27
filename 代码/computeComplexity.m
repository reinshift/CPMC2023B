%硬件复杂度计算
function C=computeComplexity(M1,M2)
    n=size(M1,1);C=0;
    for i=1:n
        for j=1:n
            for k=1:n
                if isInclude(M1(i,k),M2(k,j))
                    C=C+1;
                end
            end
        end
    end
end

function flag=isInclude(a,b)
    if any([a,b].^2==1) || any([a,b]==1+1i) || any([a,b]==1-1i) ||...
           any([a,b]==-1+1i) || any([a,b]==-1-1i) || any([a,b]==0)
        flag=0;
    else
        flag=1;
    end
end