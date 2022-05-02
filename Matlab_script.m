% FEM code for solving clamped thin plate under tension using Q4 element

close all; clear; clc

% Mesh info
coords = [0, 0; 5, 0; 10, 0; 0, 1; 5, 1; 10, 1; 0, 2; 5, 2; 10, 2];
conn = [1, 2, 5, 4; 2, 3, 6, 5; 4, 5, 8, 7; 5, 6, 9, 8];
nbc_sides = [3, 6; 6, 9];
dbc_nodes = [1, 4, 7];

num_nodes = size(coords, 1);
num_elems = size(conn, 1);

% Gauss points and weights
num_qp = 2;
q_points = [-1/sqrt(3), 1/sqrt(3)];
weights = [1, 1];

% Material properties
E = 10e3;
nu = 0.3;
h = 0.1; % thickness of the plate

% Plane stress
D = E / (1 - nu^2) * [1,  nu,           0;
                      nu, 1,            0;
                      0,  0, (1 - nu) / 2];

% Calculation of K
K = zeros(2 * num_nodes);

for elem = 1 : num_elems
    nodes_elem = conn(elem, :);
    coords_elem = coords(nodes_elem, :);
    dofs = reshape([2 * nodes_elem - 1; 2 * nodes_elem], 1, 2 * size(nodes_elem, 2));
    for s = 1 : num_qp
        for t = 1 : num_qp
            [B, J] = calculateQ4_B_J(coords_elem, q_points, s, t);
            Ke = h * B' * D * B * det(J) * weights(s) * weights(t);
            K(dofs, dofs) = K(dofs, dofs) + Ke;
        end
    end
end

% Calculation of F: Neumann BC
F = zeros(2 * num_nodes, 1);
% use of Edge2 element
for edge = 1 : size(nbc_sides, 1)
    nodes_edge = nbc_sides(edge, :);
    coords_edge = coords(nodes_edge, :);
    J = norm(diff(coords_edge)) / 2;
    for qp = 1 : num_qp
        shape = [(1 - q_points(qp)) / 2, 0, (1 + q_points(qp)) / 2, 0;
                  0, (1 - q_points(qp)) / 2, 0, (1 + q_points(qp)) / 2]; % based on shape functions of Edge2 element
        y = [(1 - q_points(qp)) / 2, (1 + q_points(qp)) / 2] * coords_edge(:, 2);
        T = [100 * (2 - y) * y; 0]; % applied traction
        Fe = h * shape' * T * J * weights(qp);
        dofs = reshape([2 * nodes_edge - 1; 2 * nodes_edge], 1, 2 * size(nodes_edge, 2));
        F(dofs) = F(dofs) + Fe;
    end
end

% Dirichlet BC
u_applied = 0;

dofs = reshape([2 * dbc_nodes - 1; 2 * dbc_nodes], 1, 2 * size(dbc_nodes, 2));
for i = 1 : size(dbc_nodes, 2)
    F = F - u_applied * K(:, dofs(i));
end
K(dofs, :) = 0;
K(:, dofs) = 0;
K(dofs, dofs) = eye(numel(dofs));
F(dofs) = u_applied;


% Solve
U = K \ F;
disps = [U(1 : 2 : end), U(2 : 2 : end)];

% Calculate strains and stresses at Gauss points within a element
% need to figure out how to visualize the straines and stresses

element = struct('strain', zeros(3, 4), 'stress', zeros(3, 4));
strain_qp = zeros(3, 4);
for elem = 1 : num_elems
    nodes_elem = conn(elem, :);
    coords_elem = coords(nodes_elem, :);
    disp_elem(1:2:7) = disps(nodes_elem, 1);
    disp_elem(2:2:8) = disps(nodes_elem, 2);
    for s = 1 : num_qp
        for t = 1 : num_qp
            [B, J] = calculateQ4_B_J(coords_elem, q_points, s, t);
            strain_qp(:, 2 * (s - 1) + t) = B * disp_elem';
        end
    end
    element(elem).strain = strain_qp;
    element(elem).stress = D * strain_qp;
end

% potential visualization issue due to linear interpolation used in patch
figure(1)
patch('Faces', conn, 'Vertices', coords, 'FaceVertexCData', disps(:, 1), 'FaceColor', 'interp')
title('Displement x')
axis equal
colormap jet
colorbar

function [gradient, jacobian] = calculateQ4_B_J(coords_elem, q_points, s, t)
gradient = zeros(3, 8);
for i = 1 : 4
    diff_st = 1 / 4 * [-1 + q_points(t), 1 - q_points(t), 1 + q_points(t), -1 - q_points(t);
                       -1 + q_points(s), -1 - q_points(s), 1 + q_points(s), 1 - q_points(s)];
    jacobian = diff_st * coords_elem;
    diff_xy = inv(jacobian) * diff_st(:, i);
    B_k = [diff_xy(1),        0;
                 0,          diff_xy(2);
           diff_xy(2),  diff_xy(1)];
    gradient(:, 2 * i - 1 : 2 * i) = B_k;
end
end
