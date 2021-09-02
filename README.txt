Bayes_Op_estim_2021-09-01.m requires the following Matlab toolboxes:
	Statistics and Machine Learning Toolbox
	The Signal Processing Toolbox
	DSP Toolbox
	Data Acquisition Toolbox
Additionally, for studies in Stieve et al., modifications were made to Bayesian.m file to keep the length-scale at a constant value:
Within the function: function GP = iFitrgpRobust(X, Y, SigmaLowerBound, varargin); starting on line 3596 in the 2016-2017 version of BayesianOptimization.m, add:
 KernelParameters_ind = find(cellfun(@(x)ischar(x)&&strcmp(x,'KernelParameters'),varargin));
        ConstantKernelParameters_ind = find(cellfun(@(x)ischar(x)&&strcmp(x,'ConstantKernelParameters'),varargin));
        constant_length_scale = 1.2; % hard coded length scale value
        if isempty(KernelParameters_ind)
            error('Didnt find KernelParameters')
        end
        varargin{KernelParameters_ind+1}(1:end-1) = constant_length_scale; % hard coded length scale for covariance
        constant_KP_replacement = false(size(varargin{KernelParameters_ind+1}));
        constant_KP_replacement(1:end-1) = true; % keep all but the last constant.
        if isempty(ConstantKernelParameters_ind)
            warning('Didnt find ConstantKernelParameters, tacking onto the end')
            varargin{length(varargin)+1} = 'ConstantKernelParameters';
            varargin{length(varargin)+1} = constant_KP_replacement;
        else
            varargin{ConstantKernelParameters_ind+1} = constant_KP_replacement;
        end


Any questions, concerns, or comments can be directed to ekrookma@umn.edu; attn: Bethany
