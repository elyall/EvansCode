classdef ROIdataset < handle
    
    properties (Constant = true)
        DataSetID = now;                % define dataset's unique ID
    end %constant properties
        
    properties
        ImagingFile
        FrameSize
        Date
        Depth
        FileLocation = nan(1,3)         % frame #, depth, channel
        Vertices
        Mask
        NeuropilMask
        RawData
        RawNeuropil
        ParsedDataSets
        Labels
        UserData
    end %main properties
    
%     properties (SetObservable = true)
% 
%     end %observable properties
    
    properties (Dependent = true)
        Pixels
        Centroid
    end %dependent properties
    
    methods
        
        % Constructor
        function obj = ROIdataset(ImagingFile, varargin)
            
            % Define imaging file dataset is associated with
            if nargin > 0
                obj.ImagingFile = ImagingFile;
            end
            
            % Define associated information
            for index = 1:2:numel(varargin)
                try
                    obj.(varargin{index}) = varargin{index+1};
                end
            end
            
        end %constructor
        
        % Set file
        function set.ImagingFile(obj,val)
            obj.ImagingFile = val;  % define file ROIdata is associated with
            
            % Define descriptors of data set
            Config = load2PConfig(val);
            obj.FrameSize = [Config.Height, Config.Width];
            obj.FileLocation(Config.size(3:end)==1) = 1;
        end
        
        % Set Vertices
        function set.Vertices(obj,val)
            obj.Vertices = val;
%             obj.Pixels(poly2mask(obj.Vertices(:,1), obj.Vertices(:,2), obj.FrameSize(1), obj.FrameSize(2)));
        end
        
%         % Set Pixels
%         function set.Pixels(obj,val)
%             obj.Pixels = val;
%             temp = regionprops(obj.Pixels, 'centroid');
%             obj.Centroid(temp.Centroid);
%         end
%         
%         % Set Centroid
%         function set.Centroid(obj,val)
%             obj.Centroid = val;
%         end
        
        % Get Pixels
        function val = get.Pixels(obj)
            val = poly2mask(obj.Vertices(:,1), obj.Vertices(:,2), obj.FrameSize(1), obj.FrameSize(2));
        end
        
        % Get Centroid
        function val = get.Centroid(obj)
            temp = regionprops(obj.Pixels, 'centroid');
            val = temp.Centroid;
        end
                
    end %methods
end %classdef