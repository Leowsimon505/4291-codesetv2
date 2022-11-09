function A = Transform2D(theta,dir)
%Global to Local, dir = 1
if (dir == 1)
    A = [cos(theta), sin(theta); -sin(theta) , cos(theta)];
 %Local to Global, dir = 0
elseif (dir == 0)
    A = [cos(theta), -sin(theta); sin(theta) , cos(theta)];
   
end