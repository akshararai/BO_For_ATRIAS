% Funtion to be used by covWarp covariance function.
%
% When init_* are specified the function creates a persistent variable the_map
% and for each row i of init_grid constructs the key and assigns:
% the_map(mkey) = init_data(i,:)
% the_map persists between calls to locomap(), but is accessible only to 
% locomap() function (so not a global variable).
%
% If only inx is specified, for each row of inx, the function looks it up in the_map.
%
% If inx and newx are specified, for each row of inx, the function constructs 
% the key and updates: the_map(mkey) = newx(i,:)
%
% the_map: map of controller to kernel value
% the_map_dim : dimension of map output - gpml needs this. Usually 1 or 2
% the_out_diff : % should the output be the difference of 2 dims. true in
% hwadjust = 2, false in all other cases
% outx is the feature transform for inx
function [outx] = locomap(inx, newx, init_grid, init_data, out_diff)
  persistent the_map; 
  persistent the_map_dim;
  persistent the_out_diff; 
  outx = NaN;
  if nargin<5, out_diff = false; end
  if nargin >= 4 %initialize map
    fprintf('locomap(): loading %d entries\n', size(init_data, 1))
    assert(size(init_data, 1) == size(init_grid,1));
    the_map = containers.Map;
    the_map_dim = size(init_data,2);
    the_out_diff = false;
    if out_diff  % use k^{v_0}_{DoG_{adjust}}
        assert(size(init_data,2) == 2);  % 2nd dim contains \bar{g}_{new}
        the_map_dim = 1;
        the_out_diff = out_diff;
    end
    % generate the map from data
    for i = 1:size(init_grid,1)
        mkey = sprintf('%0.4g', init_grid(i,:));
        the_map(mkey) = init_data(i,:);
    end
  elseif nargin == 2 %update map from hardware data
      for i = 1:size(inx,1)
        mkey = sprintf('%0.4g', inx(i,:));
        the_map(mkey) = newx(i,:);
      end
  elseif nargin == 1 %query map for feature transform
      outx = nan([size(inx,1),the_map_dim]);
      for i = 1:size(inx,1)
        mkey = sprintf('%0.4g', inx(i,:));
        if ~the_map.isKey(mkey)
            fprintf('locomap: %s not found\n', mkey)
            disp(inx)
            assert(the_map.isKey(mkey));
        end
        if the_out_diff  % \phi(x)=\phi_{sim}(x)-\bar{g}_{new}
            out = the_map(mkey);
            if (length(out) < 2), disp(out); end
            outx(i,:) = out(1)-out(2);
        else
            outx(i,:) = the_map(mkey);
        end
      end
  else  % nargin == 3
      fprintf('locomap: invalid number of input arguments\n')
      assert(false)
  end
end
