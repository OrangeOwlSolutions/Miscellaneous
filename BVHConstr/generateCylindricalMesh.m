% --- Generate a trinagular mesh of a cylindrical shape. The cylinder's
% axis is aligned to the z axis.
clc;
close all;
clear all;

lambda = 0.125;                                                         % --- Wavelength

% radius  = 0.25;                                                       % --- Radius of the cylinder
radius  = 10 * lambda;                                                  % --- Radius of the cylinder
height  = 15 * lambda;                                                  % --- Height of the cylinder
% height  = 2 * lambda;                                                 % --- Height of the cylinder

Offsetx = 30 * lambda;
Offsety = 0  * lambda;
Offsetz = 0  * lambda;

discretizationLevel = 2;                                               % --- Discretization level

numVerticalPoints   = height * discretizationLevel / lambda;            % --- Number of vertical points
% numVerticalPoints   = 360;                                            % --- Number of vertical points
numHorizontalPoints = ceil(2 * pi * radius * discretizationLevel / lambda);   
                                                                        % --- Number of horizontal points
% numHorizontalPoints = 50;                                             % --- Number of horizontal points

if (mod(numVerticalPoints,   2) ~= 0) numVerticalPoints   = numVerticalPoints   + 1; end
if (mod(numHorizontalPoints, 2) ~= 0) numHorizontalPoints = numHorizontalPoints + 1; end

numVerticalPoints   = numVerticalPoints   + 1;
numHorizontalPoints = numHorizontalPoints + 1;

% --- Total number of vertices and faces
numVertices = numVerticalPoints * numHorizontalPoints + 2 * (numVerticalPoints + 1);
numFaces    = 2 * (numVerticalPoints - 1) * (numHorizontalPoints - 1) + 2 * (numVerticalPoints - 1);

vertices = zeros(3, numVertices);                                       % --- Mesh vertices
normals  = zeros(3, numVertices);                                       % --- Mesh normals
X1       = zeros(3, numVertices);                                       % --- Principal directions at the vertices
k1       = zeros(1, numVertices);                                       % --- First principal curvatures
k2       = zeros(1, numVertices);                                       % --- Second principal curvatures
faces    = zeros(3, numFaces);                                          % --- List of the mesh faces. Each faceCount (triangle) is identified by the indices of 
                                                                        %     the three vertices.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATERAL CYLINDER SURFACE                                                         % 
% TRIANGLE VERTICES, PRINCIPAL CURVATURES AND PRINCIPAL DIRECTIONS AT THE VERTICES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vertexCount = 1;
for row = 0 : numHorizontalPoints - 1
    for col = 0 : numVerticalPoints - 1        
        angle = 2* col * pi / (numVerticalPoints - 1);
        x = radius * cos(angle);
        y = radius * sin(angle);
        z = row * height / (numHorizontalPoints - 1) - height / 2.0;            
        vertices(:, vertexCount) = [x; y; z];
        normals(:, vertexCount) = [cos(angle); sin(angle); 0];
        X1(:, vertexCount) = [0; 0; 1];
        k1(vertexCount) = 0;
        k2(vertexCount) = 1 / radius;        
        vertexCount = vertexCount + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOTTOM CAPS                                                                      % 
% TRIANGLE VERTICES, PRINCIPAL CURVATURES AND PRINCIPAL DIRECTIONS AT THE VERTICES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bottomIndex = vertexCount - 1;
vertices(:, vertexCount) = [0; 0; -height / 2.0];
normals (:, vertexCount) = [0; 0; -1];
X1(:, vertexCount) = [1; 0; 0];
k1(vertexCount) = 0;
k2(vertexCount) = 0;
vertexCount = vertexCount + 1;
for col = 0: numVerticalPoints - 1        
    angle = 2 * col *pi / (numVerticalPoints - 1);
    x = radius * cos(angle);
    y = radius * sin(angle);
    z = -height / 2.0;            
    vertices(:, vertexCount) = [x; y; z];
    normals (:, vertexCount) = [0; 0; -1];
    X1(:, vertexCount) = [1; 0; 0];
    k1(vertexCount) = 0;
    k2(vertexCount) = 0;
    vertexCount = vertexCount + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOP CAPS                                                                         % 
% TRIANGLE VERTICES, PRINCIPAL CURVATURES AND PRINCIPAL DIRECTIONS AT THE VERTICES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topIndex = vertexCount - 1;
vertices(:, vertexCount) = [0; 0; height / 2.0];
normals (:, vertexCount) = [0; 0; 1];
X1(:, vertexCount) = [1; 0; 0];
k1(vertexCount) = 0;
k2(vertexCount) = 0;
vertexCount = vertexCount + 1;
for col = 0 : numVerticalPoints - 1        
    angle = 2 * col *pi / (numVerticalPoints - 1);
    x = radius * cos(angle);
    y = radius * sin(angle);
    z = height / 2.0;            
    vertices(:, vertexCount) = [x; y; z];
    normals (:, vertexCount) = [0; 0; 1];
    X1(:, vertexCount) = [1; 0; 0];
    k1(vertexCount) = 0;
    k2(vertexCount) = 0;
    vertexCount = vertexCount + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATERAL CYLINDER SURFACE %                                                          
% FACES                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faceCount = 1;
for row = 0 : numHorizontalPoints - 2
    for col = 0 : numVerticalPoints - 2
        fx = col + row * numVerticalPoints;
        fy = col + (row + 1) * numVerticalPoints;
        fz = (col + 1) + (row + 1) * numVerticalPoints;
        faces(:, faceCount) = [fx; fy; fz];
        faceCount = faceCount + 1;
        
        fx = col + row * numVerticalPoints;
        fy = (col + 1) + (row + 1) * numVerticalPoints;
        fz = (col + 1) + row * numVerticalPoints;
        faces(:, faceCount) = [fx; fy; fz];
        faceCount = faceCount + 1;
    end
end

%%%%%%%%%%%%%%%
% BOTTOM CAPS % 
% FACES       %
%%%%%%%%%%%%%%%
for col = 0 : numVerticalPoints - 2
    fx = bottomIndex;
    fy = bottomIndex + 1 + col;
    fz = bottomIndex + 2 + col;
    faces(:, faceCount) = [fx; fy; fz];
    faceCount = faceCount + 1;  
end

%%%%%%%%%%%%
% TOP CAPS % 
% FACES    %
%%%%%%%%%%%%
for col = 0 : numVerticalPoints - 2
    fx = topIndex;
    fy = topIndex + 1 + col;
    fz = topIndex + 2 + col;
    faces(:, faceCount) = [fx; fy; fz];
    faceCount = faceCount + 1;  
end

%%%%%%%%%%%%%%%%
% SAVE TO FILE % 
%%%%%%%%%%%%%%%%
fileID = fopen('cylinder_MIMO.ply', 'w');
fprintf(fileID, 'ply\n');
fprintf(fileID, 'format ascii 1.0\n');
fprintf(fileID, 'comment Created by matlab\n');
fprintf(fileID, 'element vertexCount %d \n', numVertices);
fprintf(fileID, 'property float x\n');
fprintf(fileID, 'property float y\n');
fprintf(fileID, 'property float z\n');
fprintf(fileID, 'property float nx\n');
fprintf(fileID, 'property float ny\n');
fprintf(fileID, 'property float nz\n');
fprintf(fileID, 'property float x1x\n');
fprintf(fileID, 'property float x1y\n');
fprintf(fileID, 'property float x1z\n');
fprintf(fileID, 'property float k1\n');
fprintf(fileID, 'property float k2\n');
fprintf(fileID, 'element faceCount %d\n', numFaces);
fprintf(fileID, 'property list uchar uint vertex_indices\n');
fprintf(fileID, 'end_header\n');

for k = 1 : numVertices
%     fprintf(fileID, '%.15f %.15f %.15f  %.15f %.15f %.15f   %.15f %.15f %.15f   %.15f %.15f\n', vertices(1,k), vertices(2,k), vertices(3,k), normals(1,k), normals(2,k), normals(3,k), X1(1,k), X1(2,k), X1(3,k), k1(k), k2(k));
    fprintf(fileID, '%.15f %.15f %.15f \n', vertices(1,k), vertices(2,k), vertices(3,k));
end

for k = 1 : numFaces
    fprintf(fileID, '%d %d %d %d \n', 3, faces(1,k), faces(2,k), faces(3,k));
end

fclose(fileID);

%%%%%%%%%%%%%
% PLOT MESH % 
%%%%%%%%%%%%%
Tri = (faces.' + 1);
X = vertices(1, :);
Y = vertices(2, :);
Z = vertices(3, :);
figure
hold on
trimesh(Tri, X, Y, Z);

% --- Draw normals
for k= 1 : numVertices
    plot3([vertices(1, k), vertices(1, k) + radius / 5 * normals(1, k)], [vertices(2, k), vertices(2, k) + radius / 5 * normals(2, k)], [vertices(3, k), vertices(3, k) + radius / 5 * normals(3,k)])
end

% --- Draw principal directions
for k = 1 : numVertices
    plot3([vertices(1, k), vertices(1, k) + height / (1.5 * numHorizontalPoints) * X1(1, k)], [vertices(2, k), vertices(2, k) + height / (1.5 * numHorizontalPoints) * X1(2, k)], [vertices(3, k), vertices(3, k) + height / (1.5 * numHorizontalPoints) * X1(3, k)],'k')
end
