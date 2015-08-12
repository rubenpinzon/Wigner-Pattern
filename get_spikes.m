function varargout = get_spikes(varargin)

clusters = varargin{1};
spikes = varargin{2};
laps = varargin{3};

exVars = 0;
if nargin>3; exVars = nargin - 3; end

for j=1:max(clusters) %for each neuron  
   spk_per_neuron{j}  = spikes(clusters==j);
   for v = 1 : exVars
      eval(['aux' num2str(v) '{j} = varargin{v+3}(clusters==j);']); 
   end
   for k=1:length(laps)-1
        index = spk_per_neuron{j}>=laps(k) & spk_per_neuron{j}<laps(k+1);
        spk_per_lap{k,j} = spk_per_neuron{j}(index);
        for v = 1 : exVars
            eval(['aux_lap' num2str(v) '{k,j} = aux' num2str(v) '{j}(index);']); 
        end
   end
end

varargout{1} = spk_per_neuron;
varargout{2} = spk_per_lap;

for v = 1 : exVars
   varargout{2 + v} = eval(['aux' num2str(v)]);
   varargout{2 + exVars + v} = eval(['aux_lap' num2str(v)]);

end
