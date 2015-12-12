classdef ROI < handle
    
    properties (Hidden = true)
        roiID                           % ROI's unique identifier
    end % hidden properties
    
    properties
        Mouse                           % Mouse ROI is associated with
        DataSets                        % ROIdataset objects that contain the various datasets the ROI is found in
        UserData                        % user editable field
    end %main properties
    
    methods
        
        % Constructor
        function obj = ROI(MouseID, varargin)
            obj.roiID = now;                    % set ROI's unique ID

            % Define mouse ROI is tied to
            if nargin > 0
                obj.Mouse = MouseID;
            end
            
            % Add datasets to ROI
            if nargin > 1
                for index = 1:numel(varargin)
                    obj.DataSets(index) = ROIdataset(varargin{index});
                end
            else
                obj.DataSets = ROIdataset;
                obj.DataSets(1) = [];
            end
            
        end %constructor
        
        % Add dataset
        function obj = addDataSet(obj, ImagingFile, varargin)
            obj.DataSets(end+1) = ROIdataset(ImagingFile, varargin);
        end
        
        % Remove dataset
        function obj = removeDataSet(obj, index)
            if index < 1 || mod(index,1)~=0
                error('Value has to be a positive integer');
            elseif index > length(obj.DataSets)
                error('ROI only has %d dataset(s)', length(obj.DataSets));
            end
            obj.DataSets(index) = [];
        end
        
    end %methods
end %classdef


    