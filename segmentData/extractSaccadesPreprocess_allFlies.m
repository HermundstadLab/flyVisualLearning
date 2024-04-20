function rsout = extractSaccadesPreprocess_allFlies(rs,arena,alignment)
% EXTRACTSACCADESPREPROCESS_ALLFLIES is a wrapper script that loops over 
% flies in the input data structure, and returns a new structure with
% preprocessed data for segmentation 
% 
% INPUTS:
%   rs: data structure containing the raw heading trajectories
%   arena: a string that specifies which flight arena is used. Can take 
%       values 'upright' or '2P'
%   alignment: a string that specifies how to align the behavioral data. 
%       Can take values 'safe', 'UpBar', 'DnBar', or 'none'. Default is
%       'safe'
%
% OUTPUT:
%   rsout: a data structure indexed by fly, containing preprocessed data for
%       segmenting into saccades and fixations
%
% See also: EXTRACTSACCADESPREPROCESS


if nargin<3
    alignment = 'safe';
end

rsout = struct();
for j=1:numel(rs.fly)
    for k=1:numel(rs.fly(j).trial)
        rsout = extractSaccadesPreprocess(rs,rsout,j,k,arena,alignment);
    end
end

end