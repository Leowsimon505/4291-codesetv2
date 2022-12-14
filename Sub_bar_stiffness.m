function Ke = Sub_bar_stiffness(node_1_coord, node_2_coord, b, h)
%   Returns a 2 by 2 element stiffness matrix
%   Inputs: global x-coordinates of node 1 and node 2
%   Requires Gauss quadrature data
%   Requires E & A function for bar

  % Load Gauss data
  Gauss_data = load('Gauss_04.txt');
  Gauss_point = Gauss_data(:,1);
  Gauss_weight = Gauss_data(:,2);

  % Obtain the length of element
  L = ((node_2_coord(1)-node_1_coord(1))^2 + (node_2_coord(2)-node_1_coord(2))^2)^0.5;

  % Obtain E & A of element
  E = get_material_prop('E');
  A = b * h;

  % For Gauss quadrature sum
  sum = 0;
  for j=1:length(Gauss_point)
    % Calculate the global coordinate of Gauss point
    X = (node_2_coord(1) + node_1_coord(1))/2 + (node_2_coord(1) - node_1_coord(1))/2 * Gauss_point(j);
    Y = (node_2_coord(2) + node_1_coord(2))/2 + (node_2_coord(2) - node_1_coord(2))/2 * Gauss_point(j);
    sum = sum + Gauss_weight(j) * E * A/L^2;
  end
  sum = sum/2*L;

  Ke = sum* [1, -1; -1, 1];
  
end