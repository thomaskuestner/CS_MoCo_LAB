classdef VectorImageType < ImageType
%VECTORIMAGETYPE   class to handle nD vector images (i.e. each voxel has 3 components)
%    This class is designed to store 3 arrays with the 3 components of each voxel.
%
%   See also IMAGETYPE, VECTORIMAGETYPE.

%   Written by Alberto Gomez 2011
%   King's College London
%   OpenSource code under the BSD-2 license
%   This software is distributed with no warranties.
%   This software was designed for research purposes and not for clinical use.    
    properties(GetAccess = 'public', SetAccess = 'public')
        datax=[];
        datay=[];
        dataz=[];
    end
    
    methods(Access = public)
        %constructor
        function obj = VectorImageType(size,origin,spacing,orientation)
            if (nargin==1)
                % argument "size" is another images
                obj = VectorImageType(size.size,size.origin,size.spacing,size.orientation);
            elseif(nargin > 0)
                obj.size = size(:);
                obj.spacing =  spacing(:);
                obj.origin = origin(:);
                obj.orientation = orientation;
                s = obj.size;
                obj.data = zeros(s(:)');
                obj.datax = zeros(s(:)');
                obj.datay = zeros(s(:)');
                obj.dataz = zeros(s(:)');
                obj.paddingValue = s*0;
                obj.D = orientation;
                for i=1:numel(obj.size)
                    obj.D(:,i)=obj.D(:,i)*obj.spacing(i);
                end
            end
        end
        
        
        function P = GetValue(obj, pos, varargin)
            % get the pixel value (vector) at a non index position
            % can use different interpolation schemes:
            % im.GetValue(pos)  returns the value using nearest neighrbor
            % interpolation
            %   pos = world coordinates of the position where the value is
            %   desired
            % im.GetValue(pos, mode) uses the following interpolation
            %   mode = 'NN'     nearest neighbor
            %   mode = 'linear'    (tri) linear interpolation
            %   mode = 'spline'    (cubic) b-spline interpolation
            % im.GetValue(pos, mode,field) returns a scalar with the field
            % that can be:
            %   field = 'data'
            %   field = 'datax'
            %   field = 'datay'
            %   field = 'dataz'
            
            P=obj.paddingValue;
            mode = 'NN'; % NN
            field = '';
            if  (size(varargin,2)>0)
                mode = varargin{1};
            end
            if  (size(varargin,2)>1)
                field= varargin{2};
            end
            
            index = obj.GetContinuousIndex(pos);
            
            round_index = round(index);
            
            % find the indexes inside the range
            c_in = zeros(1,size(index,2));
            ndims = numel(obj.size);
            for i=1:ndims
                c_in = c_in | round_index(i,:)<1 | round_index(i,:)>obj.size(i);
            end
            
            
            if (strcmp(mode,'NN'))
                round_index(:,c_in) = round_index(:,c_in)*0+1;
                
                in_1D = sub2ind(obj.size',round_index(1,:),round_index(2,:),round_index(3,:));
                if numel(field) &&   any(strcmp(properties(obj), field))
                    P = obj.(field)(in_1D);
                else
                    P = [obj.datax(in_1D);obj.datay(in_1D);obj.dataz(in_1D)];
                end
                
            elseif (strcmp(mode,'linear'))
                
                index(:,c_in) = index(:,c_in)*0+1;
                P = obj.evaluate_at_point_linear(index,field);
            elseif (strcmp(mode,'spline'))
                index(:,c_in) = index(:,c_in)*0+1;
                P = obj.evaluate_at_point_spline(index,field);
            end
            
        end
        
          function out = extractFrame(obj,nframe)
            if numel(obj.size)~=4
                disp('WARNING: input image is not 4D. There might be problems')
            end
            out = VectorImageType(obj.size(1:end-1),obj.origin(1:end-1),obj.spacing(1:end-1),obj.orientation(1:end-1,1:end-1));
            out.data = obj.data(:,:,:,nframe);
            out.datax = obj.datax(:,:,:,nframe);
            out.datay = obj.datay(:,:,:,nframe);
            out.dataz = obj.dataz(:,:,:,nframe);
            
        end
    end
    methods(Access = private)
        
        function  value = evaluate_at_point_linear(obj,continuous_index,varargin)
            % continuous_index is K x N
            field = '';
            if  (size(varargin,2)>0)
                field = varargin{1};
            end
            
            
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index(' num2str(i) ',:)'];
            end
            
            if numel(field) &&  any(strcmp(properties(obj), field))
                value = eval(['interpn(obj.data' str ',''linear'')']);
            else
                valuex = eval(['interpn(obj.datax' str ',''linear'')']);
                valuey = eval(['interpn(obj.datay' str ',''linear'')']);
                valuez = eval(['interpn(obj.dataz' str ',''linear'')']);
                
                
                value = [valuex; valuey; valuez];
            end
            
            
        end
        
        
        
        function  value = evaluate_at_point_spline(obj,continuous_index,varargin)
            % continuous_index is 3 x N
            
            field = '';
            if  (size(varargin,2)>0)
                field = varargin{1};
            end
            
            
            str = '';
            for i = 1:ndims(obj.data)
                str = [str ', continuous_index(' num2str(i) ',:)'];
            end
            
            if numel(field) &&  any(strcmp(properties(obj), field))
                value = eval(['interpn(obj.data' str ',''cubic'')']);
            else
                valuex = eval(['interpn(obj.datax' str ',''cubic'')']);
                valuey = eval(['interpn(obj.datay' str ',''cubic'')']);
                valuez = eval(['interpn(obj.dataz' str ',''cubic'')']);
                
                
                value = [valuex; valuey; valuez];
            end
            
            
        end
        
    end
end
