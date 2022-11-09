function val = get_material_prop(prop_type)
% returns material properties based on property type and element type

material_prop = zeros(2,1);

% properties for hot rolled S235 structural steel
material_prop(1,1) = 210e9;                  % Young's Modulus                
material_prop(2,1) = 7850;                   % density

switch prop_type
    case 'E'
        prop_type = 1;
    case 'p'
        prop_type = 2;
end

val = material_prop(prop_type);
