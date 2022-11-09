function global_f = Input_2D_Force(nodes,elems, b, h)

%  Determine the number of nodes
node_size = size(nodes);
num_of_nodes = node_size(1);

%  Determine the number of elements
elem_size = size(elems);
num_of_elem = elem_size(1);

% Create space for vertical force
force_Y = zeros(num_of_nodes,1);
force_Y(3) = -52205*9.81/2;

% Create space for Moment
moment = zeros(num_of_nodes,1);
global_f = zeros(num_of_nodes*3,1);

for j = 1:num_of_elem

    %   Obtain global node numbers
    node_1 = elems(j,1);
    node_2 = elems(j,2);

    node_1_coord = nodes(node_1,:);
    node_2_coord = nodes(node_2,:);

    % Obtain the length of element
    L = ((node_2_coord(1)-node_1_coord(1))^2 + (node_2_coord(2)-node_1_coord(2))^2)^0.5;

    % Material Properties
    density = get_material_prop('p');
    area = b * h;

    % Calculate half element weight
    weight = area * L * density * 9.81 / 2;
    avg_X = (node_1_coord(1) + node_2_coord(1))/2;

    force_Y(node_1) = force_Y(node_1) - weight;
    force_Y(node_2) = force_Y(node_2) - weight;

    moment(node_1) = moment(node_1) + weight * (avg_X - node_1_coord(1));
    moment(node_2) = moment(node_2) + weight * (avg_X - node_2_coord(1));
end

for i=1:num_of_nodes
    global_f(3*i - 1) = force_Y(i);
    global_f(3*i) = moment(i);
end
