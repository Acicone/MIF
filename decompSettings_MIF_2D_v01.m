function options = decompSettings_MIF_2D_v01(varargin)
% 
% decompSettings_MIF_2D_v01 Constructs option structure for the algorithms
%
% EXAMPLES
%
% OPTIONS = DECOMPSETTINGS_MIF_2D_V01 with no input arguments returns
%   setting structure with default values
%
% OPTIONS = DECOMPSETTINGS_MIF_2D_V01('NAME1',VALUE1,'NAME2',VALUE2,...) creates a
%   solution options structure OPTIONS in which the named properties have
%   the specified values.  Any unspecified properties have default values.
%   It is sufficient to type only the leading characters that uniquely
%   identify the property.  Case is ignored for property names.
%
% OPTIONS = DECOMPSETTINGS_MIF_2D_V01(OLDOPTIONS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTIONS.
%
%
% DECOMPSETTINGS_V5 PROPERTIES :
%
% GENERAL
% 
%   saveEnd          - (0) If value is 1 save outputs in "name file"_decomp_MIF_2D_vXX.mat (only at termination)
%   saveIntermediate - (0) If value is 1 save outputs in "name file"_inter_decomp_MIF_2D_vXX.mat
%                      every inner step
%
%   verbose          - (1) 0 the method is silent, 1 normal, >1 loud, 
%
%   maxTime          - (Inf) Approximate max time in seconds one allows the
%                      algorithm to run. At the end of each iteration the  
%                      algorithm checks the elapsed time, if it has run for
%                      more than maxTime or maxTime - (time of last iteration)
%                      it stops iterating.
%                      Termination due to reaching of maxTime is signaled
%                      with a msg and in info.status. 
%       
%   plots            - (0) the algorithm does not produce plots,
%                      1 it produces them 
%   
%   saveplots        - (0) no saving
%                      
%
% SPECIFIC PARAMETERS
%
%  MIF.delta       (0.001) Stopping criterion
%  MIF.ExtPoints   (3)     Number of extrema allowed in the remainder
%  MIF.NIMFs       (1)     Number of IMFs we want to produce, not counting
%                         the remainder
%  MIF.Xi          (1.6)   Parameter we use to tune the mask length
%  MIF.extensionType 
%                 ('c') - constant
%                  'p'  - periodical
%                  'r'  - reflection
%  MIF.MaxInner    (200)   Max number of inner steps
%
% ------------------------------------------------------
% EXAMPLE
%          
%   >> options = decompSettings_MIF_2D_v01('MIF.delta',0.08,'MIF.NIMFs',5,'plots',1) 
%   >> IMF = Decomp_MIF_2D_v10(x,options)
%              
%  Executes algorithm MIF with delta = 0.08, stop after we have at most 5 IMFs and a trend = 5, it produces plots.                            
% ------------------------------------------------------      
%
% See also DECOMP_MIF_2D_V10
%
%
%  Ref: A. Cicone, H. Zhou. 'Multidimensional Iterative Filtering method 
%      for the decomposition of high-dimensional non-stationary signals'.
%      Preprint ArXiv http://arxiv.org/abs/1507.07173
% 

% (Ripped from sdpsettings.m by Johan Lufberg)


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    help decompSettings_MIF_2D_v01
    return;
end


Names = {
    % General   
    'saveEnd'
    'saveIntermediate'
    'verbose'
    
    'maxTime'
    'plots'
    'saveplots'
     % 'saveinIt'
     % 'logfile'
     'algorithm'
     
    % MIF
    'MIF.delta'
    'MIF.ExtPoints'
    'MIF.NIMFs'
    'MIF.Xi'
    'MIF.extensionType'
    'MIF.MaxInner'
};

obsoletenames ={ % Use when options have become obsolete
};

[m,n] = size(Names);
names = lower(Names);

if (nargin>0) && isstruct(varargin{1})
    options = varargin{1};
    paramstart = 2;
else
    paramstart = 1;
    
    % General 
    %options.saveinIt = 0;
    options.saveIntermediate = 0;
    options.saveEnd = 0;
    options.verbose = 1;
    %options.logfile = 1;
    options.maxTime = Inf;
    options.plots = 0.0;
    options.saveplots = 0; 
           
    % MIF
    options.MIF.delta = 0.001;
    options.MIF.ExtPoints=3;
    options.MIF.NIMFs=1;
    options.MIF.Xi=1.6;
    options.MIF.extensionType='c';
    options.MIF.MaxInner= 200;
end

i = paramstart;
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;       % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string property name.', i));
        end
        
        lowArg = lower(arg);
        
        j_old = strmatch(lowArg,obsoletenames);
        if ~isempty(j_old)
            % For compability... No need yet
        end
        
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized property name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous property name ''%s'' ', arg);
                msg = [msg '(' deblank(Names{j(1)})];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names{k})];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;    % we expect a value next
    else
        eval(['options.' Names{j} '= arg;']);
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for property ''%s''.', arg));
end

end

