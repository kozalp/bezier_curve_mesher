clear all;

% Init control polygon
figure;
axis([0 1 0 1]);
[x, y] = getpts();
x = x';
y = y';
canManipulatePts = false;

while (true)
    clf;
    % Plot control polygon
    plot(x, y, 'b-o');
    hold on;

    % Allocate Memory for curve
    stepSize = 0.01; % hundreds pts + 1
    u = 0: stepSize: 1;
    du = zeros(length(u),1)
    numOfU = length(u);
    c = zeros(2, numOfU);

    % Iterate over curve and apply deCasteljau
    numOfPts = length(x);
    pts = [x; y];
    for i = 1: numOfU
        ui = u(i);
        [c(:, i),da ,db] = deCasteljau(ui, pts, numOfPts, numOfPts)
        du(i) = (db(2)-da(2)) / (db(1)-da(1));
    end
    dL = zeros(1,length(c));
    for i = 1: length(c)-1
        x1 = c(1,i);
        x2 = c(1,i+1);
        y1 = c(2,i);
        y2 = c(2,i+1);
        dL(i) = ((x2-x1)^2 + (y2-y1)^2)^0.5;
    end
    
    L =sum(dL);

    % Plot curve
    axis([0 1 0 1]);
    plot(c(1, :), c(2, :), '-r');
    canManipulatePts = true;

    % Manipulate points
    if (canManipulatePts)
        pts = reposition(pts);
        x = pts(1, :);
        y = pts(2, :);
    end
end

