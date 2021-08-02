%=========================================================================================================================================================================
% This code meshes a given Bézier curve using a similar approach as in the
% paper below:
%
% https://iopscience.iop.org/article/10.1088/1748-3190/ababb0/pdf?casa_token=9cPbFTHHpkoAAAAA:WwSZzRXZmBxv-6IVI3gMkWqwbTxL88IrEMDuAE6hvp8gQVWhVkhitQ9hukRQIEheCZSwH_klCg
% 
% Source for derivative of a Bezier curve:
% https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
% 
% by Kemal OZALP
%=========================================================================================================================================================================
%% 1. INITIALIZE VARIABLES & PARAMETERS

% 1.A. Bezier Curve and mesh variables

weight_points = pts; % Weight/Control points of the Bézier curve. Output of curve.m
length_curve = L; % Length of the Bézier curve. Output of curve.m
                  % (YOU MIGHT HAVE TO CHANGE THIS WHEN IMPLEMENTING IT INTO TOBLERON)
                  
spatial_resolution = 20; % Number of pieces that you want to separate the curve into.
mesh_step_size = 1 / spatial_resolution; % 
t = 0: mesh_step_size: 1; % Dimensionless Bézier curve coordinates, 0=< t_i =< 1
num_t = length(t); % Number of points on the Bézier curve                  
dL = zeros(1, num_t - 1); % Array that stores length of each piece of line 
                          % making up the Bézier curve
                          
num_weight_points = length(pts); % Number of weight points
bezier_curve_points = zeros(2, num_t); % Matrix to store the coordinates of the Bézier curve points
derivative_t = zeros(1, num_t); % Matrix to store derivatives of points on the Bézier curve

ds = length_curve / spatial_resolution; % Desired length of each each mesh piece

% 1.B. Cost function and loop variables
                         
dJ = zeros(1, num_t); % Derivative of the cost function 
alpha = 5*10^(-5); % Learning rate
temp = 1; % 
iter = 10^5; % number of iterations

%% 2. CALCULATE BÉZIER CURVE POINTS & LENGTH OF EACH LINE PIECE
% Instead of using an error parameter as a threshold, I used number of
% iterations to construct the loop

while temp < iter
    
% Iterate over curve and apply deCasteljau

    for i = 1: num_t
        ti = t(i);
        [bezier_curve_points(:, i),a ,b] = deCasteljau(ti, weight_points,...
            num_weight_points, num_weight_points);
        
        % Here, a & b are the end point coordinates of the last line
        % created to approximate the point on the Bézier curve via de
        % Casteljau algorithm. For more information please see the link:
        % https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html
        derivative_t(i) = (b(2)-a(2)) / (b(1)-a(1));
    end

% Calculate the distance between each point on the Bézier curve


    for i = 1: length(bezier_curve_points) - 1
        x1 = bezier_curve_points(1,i);
        x2 = bezier_curve_points(1,i+1);
        y1 = bezier_curve_points(2,i);
        y2 = bezier_curve_points(2,i+1);
        dL(i) = ((x2-x1)^2 + (y2-y1)^2)^0.5;
    end

%% 3. CALCULATE COST FUNCTION (J(t)) & UPDATE BÉZIER CURVE POINTS


    % Cost function

    err = ((dL - ds) / ds).^2;
    J = sum(err) / (2 * (num_weight_points));

    % Calculate new points on the Bézier curve
    % I didn't use the exact equations for the derivative of the cost function 
    % (J(t)). Equations can be found in this paper on page 6:
    % https://iopscience.iop.org/article/10.1088/1748-3190/ababb0/pdf?casa_token=9cPbFTHHpkoAAAAA:WwSZzRXZmBxv-6IVI3gMkWqwbTxL88IrEMDuAE6hvp8gQVWhVkhitQ9hukRQIEheCZSwH_klCg
    % The reason is that I though the start and the end points for
    % dimensionless Bezier curve coordinates (t(i)) shouldn't change. So, I 
    % set dJ(1) to zero and didn't set any rules for dJ(n). This is where
    % things might be going wrong.
    
    for j = 1: num_t

        if j == 1 || j == num_t
            dJ(j) = 0; 
            %dJ(j) = -1/(num_weight_points-1)*ds^2 * (dL(j)-ds)/dL(j) * dL(j)*derivative_t(j);
        elseif j == length(bezier_curve_points) - 1
            dJ(j) = 1/(num_weight_points-1)*ds^2 * (dL(j-1)-ds)/dL(j-1) * dL(j-1)*derivative_t(j);
        else
            dJ(j) = 1/(num_weight_points-1)*ds^2 * ((dL(j-1)-ds)/dL(j-1) * dL(j-1)*derivative_t(j)) - ((dL(j)-ds)/dL(j) * dL(j)*derivative_t(j));
        end

    end

    t = t - alpha * dJ;
    
    temp = temp + 1;
end

%%4. CALCULATE FINAL BÉZIER CURVE POINTS & PLOT CURVE AND MESH

for i = 1: num_t
        ti = t(i);
        [bezier_curve_points(:, i),~ ,~] = deCasteljau(ti, weight_points,...
            num_weight_points, num_weight_points);
end

 % Plot curve
figure;
axis([0 1 0 1]);
hold on 
plot(c(1, :), c(2, :), '-r');
plot(bezier_curve_points(1, :), bezier_curve_points(2, :), 'xb');
hold off