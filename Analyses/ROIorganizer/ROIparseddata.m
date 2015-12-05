classdef ROIparseddata < handle
    properties
        ParseDataID
        Data            % (Ntrials x Nframes)
        Neuropil        % (Ntrials x Nframes)
        dFoF            % (Ntrials x Nframes)
        NumFrames       % (Ntrials x 3) [before, stimulus, after]
        StimulusInfo    % (Ntrials x Nframes) logical 1 when stimulus present, 0 when stimulus not present
        StimID          % (Ntrials x 1)
        TrialID         % (Ntrials x 1)
        UserData        % user editable field
    end
    
    methods
        
        % Constructor
        function obj = ROIparseddata(ROIdataset, Experiment, NumFrames)
            obj.ParseDataID = now; % define parseddata's unique ID
            
            % Parse data
            
            
        end %constructor
        
    end %methods
end %classdef