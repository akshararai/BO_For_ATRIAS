% Funtion to be used by meanMmap function.
%
% When init_* are specified the function creates a persistent variable the_mmap
% and for each row i of init_grid constructs the key and assigns:
% the_mmap(mkey) = init_data(i,:)
% the_mmap persists between calls to mmap(), but is accessible only to 
% mmap() function (so not a global variable).
%
% If only inx is specified, for each row of inx, the function constructs 
% the key and looks it up in the_mmap.
%
% If inx and newx are specified, for each row of inx, the function constructs 
% the key and updates: the_mmap(mkey) = newx(i,:)
function [outx] = mmap(inx, newx, init_grid, init_data)
  persistent the_mmap;
  persistent the_mmap_dim;
  outx = NaN;
  if nargin >= 4
    fprintf('mmap(): loading %d entries\n', size(init_data, 1))
    assert(size(init_data, 1) == size(init_grid,1));
    the_mmap = containers.Map;
    the_mmap_dim = size(init_data,2);
    for i = 1:size(init_grid,1)
        mkey = sprintf('%0.4g', init_grid(i,:));
        the_mmap(mkey) = min(init_data(i,:),120);
    end
  elseif nargin == 2
      for i = 1:size(inx,1)
        mkey = sprintf('%0.4g', inx(i,:));
        the_mmap(mkey) = newx(i,:);
      end
  elseif nargin == 1
      outx = nan([size(inx,1),the_mmap_dim]);
      for i = 1:size(inx,1)
        mkey = sprintf('%0.4g', inx(i,:));
        if ~the_mmap.isKey(mkey)
            fprintf('mmap: %s not found\n', mkey)
            disp(inx)
            assert(the_mmap.isKey(mkey));
        end
        outx(i,:) = the_mmap(mkey);
      end
  else  % nargin == 3
      fprintf('mmap: invalid number of input arguments\n')
      assert(false)
  end
end
