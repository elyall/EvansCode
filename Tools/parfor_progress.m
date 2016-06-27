function [Fn,percent] = parfor_progress(varargin)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   FILE = PARFOR_PROGRESS(N) initializes the progress monitor for a set of
%   N upcoming calculations.
%
%   PARFOR_PROGRESS(FILE) updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(FILE, 0) deletes parfor_progress.txt and finalizes
%   progress bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      F = parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress(F);
%      end
%      parfor_progress(F,0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/
% Edited by Evan Lyall

narginchk(0,2);

percent = 0;
w = 50; % Width of progress bar

if nargin == 0
    
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
    
elseif nargin == 1
    
    if isnumeric(varargin{1}) && varargin{1} ~= 0 % create and open file
        
        if ~exist('parfor_progress.txt','file')
            Fn = 'parfor_progress.txt';
        else
            index = 1;
            while exist(sprintf('parfor_progress%d.txt',index),'file')
                index = index + 1;
            end
            Fn = sprintf('parfor_progress%d.txt',index);
        end
        
        f = fopen(Fn, 'w');
        if f<0
            error('Do you have write permissions for %s?', pwd);
        end
        fprintf(f, '%d\n', varargin{1}); % Save N at the top of progress.txt
        fclose(f);
        
        if nargout <= 1
            disp(['  0%[>', repmat(' ', 1, w), ']']);
        end
        
    elseif isnumeric(varargin{1}) && varargin{1} == 0 % delete file
        
        delete('parfor_progress.txt');
        percent = 100;
        
        if nargout == 0
            disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
        end
        
    elseif ischar(varargin{1}) % increment file
        
        if ~exist(varargin{1}, 'file')
            error('%s not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.', varargin{1});
        end
        
        f = fopen(varargin{1}, 'a');
        fprintf(f, '1\n');
        fclose(f);
        
        f = fopen(varargin{1}, 'r');
        progress = fscanf(f, '%d');
        fclose(f);
        percent = (length(progress)-1)/progress(1)*100;
        
        if nargout <= 1
            perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
            disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
        end
        
    end
    
elseif nargin == 2 % delete file

    delete(varargin{1});
    percent = 100;
    
    if nargout <= 1
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
    
end
