function [point, a, b] = deCasteljau(u, pts, i, j)
    if i == 1
        point = pts(:, j);
    else
        %point = u * deCasteljau(u, pts, i - 1, j) + (1 - u) * deCasteljau(u, pts, i - 1, j - 1);
        a = deCasteljau(u, pts, i - 1, j);
        b = deCasteljau(u, pts, i - 1, j - 1);
        point = u * a + (1 - u) * b;
    end
end

