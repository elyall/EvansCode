classdef ROI < handle
    
    properties (Constant = true)
        roiID = now;                    % set ROI's unique ID
    end %constant properties
    
    properties (GetAccess = private)
        numDataSets = 0;
    end %private properties
    
    properties
        Mouse                           % Mouse ROI is associated with
        DataSets                        % ROIdataset objects that contain the various datasets the ROI is found in
        UserData                        % user editable field
    end %main properties
    
    methods
        
        % Constructor
        function obj = ROI(MouseID, varargin)
            
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
            end
            
        end %constructor
        
        % Add dataset
        function obj = addDataSet(obj, ImagingFile, varargin)
            obj.DataSets(obj.numDataSets+1) = ROIdataset(ImagingFile, varargin);
            obj.numDataSets = obj.numDataSets + 1;
        end
        
        % Remove dataset
        function obj = removeDataSet(obj, index)
            obj.DataSets(index) = [];
            obj.numDataSets = obj.numDataSets - 1;
        end
        
    end %methods
end %classdef


    