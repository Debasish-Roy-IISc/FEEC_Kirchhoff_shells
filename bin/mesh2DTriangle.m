classdef mesh2DTriangle < handle
    %Creates a mesh on the rectangular domain [0, l] x [0, w].
    properties
        nodes
        elements
        noEl
        noNd
    end
    methods
        function obj = mesh2DTriangle(l,w,nx,ny)
            x=linspace(l(1),l(2),nx);
            y=linspace(w(1),w(2),ny);
            [X,Y] = meshgrid(x,y);
            obj.nodes=[X(:) Y(:)];
            obj.elements=delaunay(obj.nodes(:,1),obj.nodes(:,2));
            obj.noEl=size(obj.elements,1);
            obj.noNd=size(obj.nodes,1);
        end
    end
end