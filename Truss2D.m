nodes = load('node_2.txt');
elems = load('element_2.txt');

%  Split the elements recursively
[nodes,elems] = elementsplit(nodes,elems,0);
num_of_elem = size(elems,1);

dimensions = load('dimensions.txt');
dimensions_size = size(dimensions,1);

%max stress for each set of dimensions
stress_array = zeros(num_of_elem,dimensions_size);

for k = 1:dimensions_size
%  Determine the number of nodes
node_size = size(nodes);
num_of_nodes = node_size(1);

%  Allocate space for [K] matrix and [F] vector
K = zeros(3*num_of_nodes, 3*num_of_nodes);

%  Determine the number of elements
elem_size = size(elems);
num_of_elem = elem_size(1);

%input dimensions b & h
b = dimensions(k,1); % element breadth
h = dimensions(k,2); % element height

%Assembly
for j = 1:num_of_elem
    
    %   Obtain global node numbers
    node_1 = elems(j,1);
    node_2 = elems(j,2);


    %transformation
    theta = 0;
    if (nodes(node_2,1) - nodes(node_1,1)< 0)
        theta = theta+pi;
    end
    theta = theta+ atan((nodes(node_2,2) - nodes(node_1,2))/ (nodes(node_2,1) - nodes(node_1,1)));
    A2 = zeros(6,6);

    %theta = angle from x axis
    A2([1,2],[1,2]) = Transform2D( theta,0 );
    A2(3,3) = 1;
    A2([4,5],[4,5]) = Transform2D( theta,0 );
    A2(6,6) = 1;
    
    %	Flocal = A2Inv * K * A2 * ulocal
    K_bar = Sub_bar_stiffness(nodes(node_1,:), nodes(node_2,:), b, h);
    K_bar_global = A2 * Input_1D_k(K_bar) * A2';
    
    K_beam = Sub_beam_stiffness(nodes(node_1,:), nodes(node_2,:), b, h);
    K_beam_global = A2 * Input_1D_k_beam(K_beam) * A2';
    K_elem_global = K_bar_global + K_beam_global;
   
    %divide K_bar_global into 4 matrices (6x6)
    first_first = K_elem_global(1:3,1:3);
    first_second = K_elem_global(1:3,4:6);
    second_first = K_elem_global(4:6,1:3);
    second_second = K_elem_global(4:6,4:6);
    
    elem_global_dof_1 = (3*node_1 - 2 : 3*node_1);
    elem_global_dof_2 = (3*node_2 - 2 : 3*node_2);
    
    K( elem_global_dof_1, elem_global_dof_1 ) = K( elem_global_dof_1, elem_global_dof_1 ) + first_first;
    K( elem_global_dof_1, elem_global_dof_2 ) = K( elem_global_dof_1, elem_global_dof_2 ) + first_second;
    K( elem_global_dof_2, elem_global_dof_1 ) = K( elem_global_dof_2, elem_global_dof_1 ) + second_first;
    K( elem_global_dof_2, elem_global_dof_2 ) = K( elem_global_dof_2, elem_global_dof_2 ) + second_second;
    
end

% Load force and moment 
F = Input_2D_Force(nodes,elems, b, h);

%  Apply boundary condition
disp_data = load('boundary condition.txt');
disp_nodes = disp_data(:,1);
disp_dof = disp_data(:,2);
disp_values = disp_data(:,3);
num_of_disp = length(disp_values);

for j = 1:num_of_disp
    % Determine global dof
    if (disp_dof(j) == 1)
        global_dof = 3*disp_nodes(j) - 2; % X direction
    elseif (disp_dof(j) == 2) 
        global_dof = 3*disp_nodes(j) - 1;  % Y direction
    else
        global_dof = 3*disp_nodes(j);
    end   
    F(global_dof) = disp_values(j);
    K(global_dof,:) = zeros(1,3*num_of_nodes);
    K(global_dof,global_dof) = 1; 
end

%  Solve for the displacements
u = K\F;
%input material properties
elem_stress = zeros(num_of_elem,1);

for i = 1:num_of_elem
    %determine which nodes for each element
    node_1 = elems(i,1);
    node_2 = elems(i,2);

    %extract material properties
    E = get_material_prop('E');
    h = h;

    %extracting undeformed node coordinates
    X1 = nodes(node_1, 1);
    Y1 = nodes(node_1, 2);
    X2 = nodes(node_2, 1);
    Y2 = nodes(node_2, 2);
    L_undeformed = sqrt((Y2-Y1)^2 + (X2-X1)^2);

    %calculating undeformed rotation of element
    r_elem = atan2((Y2-Y1) , (X2-X1));
    %calculating undeformed tangent of element
    t_elem = (Y2-Y1) / (X2-X1);
    
    %extracting nodal displacements X, Y, rotation
    D1 = u(3*node_1-2 : 3*node_1);
    D2 = u(3*node_2-2 : 3*node_2);

    %calculating deformed node coordinates
    x1 = X1 + D1(1);
    y1 = Y1 + D1(2);
    x2 = X2 + D2(1);
    y2 = Y2 + D2(2);
    L_deformed = sqrt((y2-y1)^2 + (x2-x1)^2);

    %calculating deformed tangent at each node
    t1 = tan(r_elem + D1(3));
    t2 = tan(r_elem + D2(3));

    %calculating stress
    elem_bending_stress = calc_bending_stress(x1, y1, t1, x2, y2, t2, E, h);
    elem_axial_stress = E * abs(L_deformed - L_undeformed)/L_undeformed;
    stress_array(i,k) = (elem_axial_stress + elem_bending_stress);

end

end

%array of areas based on input dimensions
area = zeros(1,dimensions_size);
for i=1:dimensions_size
    area(i) = dimensions(i,1)*dimensions(i,2);
end

%number of rows and columns of stress array
[r,c] = size(stress_array);

%array variables
stress = zeros(1,c);
min_b = zeros(r,1);
min_h = zeros(r,1);
cross_section = zeros(r,1);

for j=1:r
    stress = stress_array(j,:);

    x = area(:);
    y = stress(:);
    
    %obtain curvefit for element
    f = fit(x,y,'power1');
    
    %plot(f)
    %title('Element' + string(j));
    %figure;
    
    %equation coefficient
    a=f.a;
    b=f.b;
    
    syms x y z 
    assume(x > 0)
    M = [a*x^b == 117.5e6;];
    [xsol] = vpasolve(M, [x], 'random',true);
    
    %populate variables with minimum dimensions
    cross_section(j) = xsol;
    min_h(j) = sqrt(xsol/1.5);
    min_b(j) = 1.5*min_h(j);
end
