classdef EchoImageType < ImageType
%ECHOIMAGETYPE   class to handle nD echo images
%    This class provides some convenient features for ultrasound images (particularly 3D US), like finding the transducer position from pyramidal frustums and storing B-Mode and colour Doppler data.
%
%   See also IMAGETYPE, VECTORIMAGETYPE.

%   Written by Alberto Gomez 2011
%   King's College London
%   OpenSource code under the BSD-2 license
%   This software is distributed with no warranties.
%   This software was designed for research purposes and not for clinical use.
    
    
    properties(GetAccess = 'public', SetAccess = 'public')
        beamSource = [0 0 0]';
        velocityBounds = [-100 100];
        dataVelocity = [];
    end
    
    methods(Access = public)
        %constructor
        function obj = EchoImageType(size,origin,spacing,orientation)
            if (nargin==1)
                % argument "size" is another images
                obj = EchoImageType(size.size,size.origin,size.spacing,size.orientation);
            elseif (nargin==2)
                % two arguments: first gray image, then colour
                obj = EchoImageType(size.size,size.origin,size.spacing,size.orientation);
                obj.data = size.data;
                obj.dataVelocity = origin.data;
                
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                obj.data = zeros(obj.size(:)');
                obj.dataVelocity = zeros(obj.size(:)');
                obj.paddingValue = zeros(numel(size),1);
                obj.beamSource = [];
                obj.velocityBounds = [-100 100];
                
                obj.D = orientation;
                
                for i=1:numel(obj.size)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);    
                end
            end
        end
        
        function [out_mesh direction angles radiuses]= GetFrustumMesh(obj,varargin)
            % generates a mesh which is the frustum. Also returns the axis
            % of the frustum (i.e. the axis of the mesh)
            % resolution of the spherical caps is fixed to 8 in each
            % direction
            N=8;
            dgb=false;
            colorbox=false; % this flag says wether the frustum is computed from the bmode or from the color
            MAX_CHUNK_SIZE = 50;
            for i=1:size(varargin,2)
                if (strcmp(varargin{i},'debug'))
                    dbg=true;
                elseif (strcmp(varargin{i},'N'))
                    N = varargin{i+1};
                    i=i+1;
                elseif (strcmp(varargin{i},'colorbox'))
                    colorbox=true;
                end
                
            end
            
            
            
            if ~numel(obj.beamSource)
                obj.GetBeamSource();
            end
            
            notOrientedImage = ImageType(obj);
            if colorbox
                notOrientedImage.data = obj.dataVelocity;
            else
                notOrientedImage.data = obj.data;
            end
            refdata = obj.data;
            notOrientedImage.orientation = eye(3);
            notOrientedImage.D = notOrientedImage.orientation;
            notOrientedImage.D(:,1)=notOrientedImage.D(:,1)*notOrientedImage.spacing(1);
            notOrientedImage.D(:,2)=notOrientedImage.D(:,2)*notOrientedImage.spacing(2);
            notOrientedImage.D(:,3)=notOrientedImage.D(:,3)*notOrientedImage.spacing(3);
            notOrientedImage.origin = -1*(obj.size-1).*obj.spacing/2;
            
            % find center of mass of the image
            nonzerosindices = find(abs(notOrientedImage.data)>0);
            [x, y, z] = ind2sub(notOrientedImage.size', nonzerosindices );
            positions = obj.GetPosition([x(:) y(:) z(:)]');
            meanpos = mean(positions,2);
            
            meanIndex = obj.GetContinuousIndex(meanpos );
            
            slicex=round(meanIndex(1));
            slicey=round(meanIndex(2));
            
            % extract the two central slices
            slice1 = squeeze(refdata(:,slicey,:));
            [i1, angle1, angle2] = obj.findSliceSource(slice1,[1 1],false);
            % correct for the spacing
            angle1 = atan(tan(angle1*pi/180)*obj.spacing(1)/obj.spacing(3));
            angle2 = atan(tan(angle2*pi/180)*obj.spacing(1)/obj.spacing(3));
            
            
            slice2 = squeeze(refdata (slicex,:,:));
            [ i2, angle3, angle4 ] = obj.findSliceSource(slice2,[1 1],false);
            % correct for the spacing
            angle3 = atan(tan(angle3*pi/180)*obj.spacing(2)/obj.spacing(3));
            angle4 = atan(tan(angle4*pi/180)*obj.spacing(2)/obj.spacing(3));
            
            
            index=[i1(2) i2(2) (i1(1)+i2(1))/2]';
            
            
            nonOrientedSource1 = [ notOrientedImage.GetPosition(index); 1];
            nonOrientedSource = nonOrientedSource1-[notOrientedImage.origin(:); 0 ];
            M= [obj.orientation obj.origin(:); 0 0 0 1];
            orientedSource = M * nonOrientedSource;
            new_beamSource= orientedSource(1:3,1);
            
            % Get the maximum radius
            nonzeropositions = positions-new_beamSource*ones(1,numel(nonzerosindices));
            clear positions;
            
            
            
            radiuses = sqrt(nonzeropositions(1,:).^2+nonzeropositions(2,:).^2+nonzeropositions(3,:).^2);
            
            max_radius = max(radiuses );
            min_radius = min(radiuses );
            
            
            anglex_span = (angle2-angle1)/N;
            angley_span = (angle4-angle3)/N;
            
            angles_x = angle1:anglex_span:angle2;
            angles_y = angle3:angley_span:angle4;
            sp_counter = 0;
            for ax = angles_x
                for ay = angles_y
                    sp_counter = sp_counter+1;
                    %surface_points(sp_counter,:) = [sin(ax) cos(ax)*sin(ay) cos(ax)*cos(ay) ];
                    surface_points(sp_counter,:) = [sin(ax)*cos(ay) cos(ax)*sin(ay) sqrt(1-cos(ax).^2*sin(ay).^2-cos(ay).^2*sin(ax).^2) ];
                end
            end
            
            % topology of the surfaces
            
            points_surfaces_2D = surface_points(:,1:2);
            DT = DelaunayTri(points_surfaces_2D);
            triangles_surface = DT.Triangulation;
            
            points = [ min_radius*surface_points; max_radius*surface_points];
            
            %  indices of points in the sides of the surfacein  counter-clockwise order from the beam source
            indices_of_surfaces_sides_down = [ 1:numel(angles_y):((numel(angles_x)-1)*numel(angles_y)+1) ...
                ((numel(angles_x)-1)*numel(angles_y)+1):1:(numel(angles_x)*numel(angles_y)) ...
                (numel(angles_x)*numel(angles_y)):-numel(angles_y):numel(angles_x) ...
                numel(angles_x):-1:1];
            indices_of_surfaces_sides_down(numel(angles_y))=[];
            indices_of_surfaces_sides_down(numel(angles_y)+numel(angles_x)-1)=[];
            indices_of_surfaces_sides_down(2*numel(angles_y)+numel(angles_x)-2)=[];
            % indices_of_surfaces_sides_down(end)=[];
            
            indices_of_surfaces_sides_up = indices_of_surfaces_sides_down+numel(angles_x)*numel(angles_y);
            
            triangles = [triangles_surface; triangles_surface+size(surface_points,1)];
            
            for i=1:(numel(indices_of_surfaces_sides_down)-1)
                triangles = [triangles
                    indices_of_surfaces_sides_down(i) indices_of_surfaces_sides_down(i+1) indices_of_surfaces_sides_up(i)
                    indices_of_surfaces_sides_down(i+1) indices_of_surfaces_sides_up(i+1) indices_of_surfaces_sides_up(i) ];
            end
            
            
            out_mesh = MeshType(size(points,1),size(triangles,1));
            out_mesh.triangles = triangles;
            % the points have to be transformed to the appropiate orientation and origin
            
            
            
            points = [points' ; ones(1,size(points,1))]- [notOrientedImage.origin(:)-nonOrientedSource1(1:3); 0 ]*ones(1,size(points,1))	;
            orientedpoints = M *points ;
            
            out_mesh.points = orientedpoints(1:3,:)';
            
            direction = obj.orientation;
            angles = [angle1 angle2 angle3 angle4]*180/pi;
            radiuses = [min_radius max_radius];
            
            
        end
        
        
        function angles= GetFrustumAngles(obj)
            % note: the angles given are with respect to the non oriente
            % version!
            notOrientedImage = ImageType(obj);
            notOrientedImage.data = obj.data;
            notOrientedImage.orientation = eye(3);
            notOrientedImage.D = notOrientedImage.orientation;
            notOrientedImage.D(:,1)=notOrientedImage.D(:,1)*notOrientedImage.spacing(1);
            notOrientedImage.D(:,2)=notOrientedImage.D(:,2)*notOrientedImage.spacing(2);
            notOrientedImage.D(:,3)=notOrientedImage.D(:,3)*notOrientedImage.spacing(3);
            notOrientedImage.origin = -1*(obj.size-1).*obj.spacing/2;
            
            % find center of mass of the image
            [x y z]=ndgrid(1:obj.size(1),1:obj.size(2), 1:obj.size(3));
            positions = obj.GetPosition([x(:) y(:) z(:)]');
            
            meanx = sum(positions(1,:).*obj.data(:)')/sum(obj.data(:));
            meany = sum(positions(2,:).*obj.data(:)')/sum(obj.data(:));
            meanz = sum(positions(3,:).*obj.data(:)')/sum(obj.data(:));
            
            meanIndex = obj.GetContinuousIndex([meanx meany meanz]');
            
            slicex=round(meanIndex(1));
            slicey=round(meanIndex(2));
            
            % extract the two central slices
            slice1 = squeeze(notOrientedImage.data(:,slicey,:));
            [i1 angle1 angle2] = obj.findSliceSource(slice1,[ 1 1]);
            
            % convert i1 to continuous point
            
            slice2 = squeeze(notOrientedImage.data(slicex,:,:));
            [ i2 angle3 angle4 ] = obj.findSliceSource(slice2,[1 1]);
            
            % rotate the angles with orientation
            angles = [angle1 angle2 angle3 angle4];
            
        end
        
        function  bs = GetBeamSource(obj,varargin)
            
            dbg = false;
            if numel(varargin)==1;
                dbg=varargin{1};
            end
            
            % get the beam source location, using the hough transform
            notOrientedImage = ImageType(obj);
            notOrientedImage.data = obj.data;
            notOrientedImage.orientation = eye(3);
            notOrientedImage.D = notOrientedImage.orientation;
            notOrientedImage.D(:,1)=notOrientedImage.D(:,1)*notOrientedImage.spacing(1);
            notOrientedImage.D(:,2)=notOrientedImage.D(:,2)*notOrientedImage.spacing(2);
            notOrientedImage.D(:,3)=notOrientedImage.D(:,3)*notOrientedImage.spacing(3);
            notOrientedImage.origin = -1*(obj.size-1).*obj.spacing/2;
            
            % find center of mass of the image
            [x y z]=ndgrid(1:obj.size(1),1:obj.size(2), 1:obj.size(3));
            positions = obj.GetPosition([x(:) y(:) z(:)]');
            if ~sum(abs(obj.data(:)))
                obj.beamSource= [NaN NaN NaN]';
                disp('WARNING: image is all black');
                return;
            end
            meanx = sum(positions(1,:).*abs(obj.data(:))')/sum(abs(obj.data(:)));
            meany = sum(positions(2,:).*abs(obj.data(:))')/sum(abs(obj.data(:)));
            meanz = sum(positions(3,:).*abs(obj.data(:))')/sum(abs(obj.data(:)));
            
            meanIndex = obj.GetContinuousIndex([meanx meany meanz]');
            
            slicex=round(meanIndex(1));
            slicey=round(meanIndex(2));
            
            % extract the two central slices
            slice1 = squeeze(notOrientedImage.data(:,slicey,:));
            
            i1 = obj.findSliceSource(slice1,[1 1],dbg);
            
            % convert i1 to continuous point
            
            slice2 = squeeze(notOrientedImage.data(slicex,:,:));
            i2 = obj.findSliceSource(slice2,[1 1],dbg);
            % convert i1 to continuous point
            
            index=[i1(2) i2(2) (i1(1)+i2(1))/2]';
            
            nonOrientedSource = [ notOrientedImage.GetPosition(index); 1];
            nonOrientedSource = nonOrientedSource-[notOrientedImage.origin(:); 0 ];
            orientedSource = [obj.orientation obj.origin(:); 0 0 0 1] * nonOrientedSource;
            
            obj.beamSource= orientedSource(1:3,1);
            bs = orientedSource(1:3,1);
        end
    end
    methods(Access = private)
        function [ point angle1 angle2] = findSliceSource(obj,slice,spacing,dbg)
            % retrns the slice source and the angle of the in-plane frustum
            n=40;
            slice(abs(slice)>=1)=100;
            slice(abs(slice)<1)=0;
            BW = edge(slice,'canny');
            [H,T,R] = hough(BW,'RhoResolution',0.6,'Theta',[-89:0.8:-35 35:0.8:89]);
            % find the two highes points
            P = houghpeaks(H,2,'Threshold',0.25*max(H(:)),'NHoodSize',ceil(size(H)/40)*2+1);
            lines = houghlines(BW,T,R,P,'FillGap',50,'MinLength',4);
            
            % find two theta
            spacing = spacing(:)';
            lines(1).point1 = lines(1).point1.*spacing;
            u1 = [ lines(1).point1(1)-lines(1).point2(1) lines(1).point1(2)-lines(1).point2(2)];
            u1 = u1/norm(u1);
            u2 = [ lines(2).point1(1)-lines(2).point2(1) lines(2).point1(2)-lines(2).point2(2)];
            u2 = u2/norm(u2);
            
            vec = [0 1];
            a1 = acos(u1*vec')*180/pi;
            a2 = acos(u2*vec')*180/pi;
            
            angle1 = min([a1 a2])-90;
            angle2 = max([a1 a2])-90;
            
            p1 = lines(1).point1;
            p2 = lines(2).point1;
            
            A = [u1(1) -u2(1); u1(2) -u2(2)];
            b = [p2(1)-p1(1); p2(2)-p1(2)];
            
            x = A\b;
            % find intersection of both lines
            point = p1 + u1*x(1);
            
            % for testing
            if dbg
                figure;
                imagesc(BW);
                hold on;
                line([lines(1).point1(1); lines(1).point2(1)],[lines(1).point1(2); lines(1).point2(2)],'Color',[1 0 0],'LineWidth',4);
                line([lines(2).point1(1); lines(2).point2(1)],[lines(2).point1(2); lines(2).point2(2)],'Color',[0 1 0],'LineWidth',4);
                hold off;
                axis equal;
            end
            
            
            
        end
        
    end
end
