function v1 = box_boundary(tree,varargin)
  p = inputParser;
  depth = tree.depth; n1 = tree.count(depth); 
  
  defaultInd = (1:n1)';
  defaultLimit = 1;
  addRequired(p,'tree');
  addParamValue(p,'Indices',defaultInd,@isvector);
  addParamValue(p,'Limit',defaultLimit,@isnumeric);
  
  parse(p, tree, varargin{:});
  ind = p.Results.Indices;
  limit = p.Results.Limit;
  
  b = tree.boxes(depth); d = tree.dim; b = b(:,ind); r = b(d+1:2*d,1); 
  r0 = tree.boxes(0); c = r0(1:d,1); r0 = r0(d+1:2*d); tree = Tree(c,r0);
  tree.insert(b(1:d,:),depth);
  n = size(b,2); v = zeros(n,1);
  for i = 1:n
      x = b(1:d,i);
      if limit, if any(abs(x-c)>=r0-2*r), continue; end; end
      if length(tree.search_box(x,2*r))<3^d
          v(i) = 1;
      end
  end
  v1 = zeros(n1,1);
  v1(ind) = v;
end