function t = anova2D(y,group)

if (size(y,1)==1), y=y(:); end

group = group(:);
ng    = length(group);

% Calling sequence with named arguments
model    = 1; %1: 'linear'
sstype   = 2;
varnames = cellstr([repmat('X',ng,1) strjust(num2str((1:ng)'),'left')]);
display  = 'off';
alpha    = 0.05;

okargs   = {'random'     'continuous' 'nested'};
defaults = {false(ng,1)  []           []};
[randomvar continuous vnested] = internal.stats.parseArgs(okargs,defaults);

continuous = ismember(1:ng,continuous);
randomvar  = randomvar(:);
% Check optional arguments
% variable names
termlist = eye(ng);

% STEP 1:  Remove NaN's and prepare grouping variables.
% Also, make sure all groups are rows, and create group index and name arrays.
% Find NaNs among response and group arrays
n = size(y,1);
nanrow = any(isnan(y),2);
ng = length(group);
for j=1:ng
   gj = group{j};
   if (size(gj,1) == 1), gj = gj(:); end
   if (size(gj,1) ~= n)
      error(message('stats:anovan:GroupVarSizeMismatch', j, n));
   end
   if (ischar(gj)), gj = cellstr(gj); end
   if ~isvector(gj)
       error(message('stats:anovan:GroupNotVectorOrCharArray'));
   end

   group{j} = gj;
   if (isnumeric(gj))
      nanrow = (nanrow | isnan(gj));
   elseif isa(gj,'categorical')
      nanrow = (nanrow | isundefined(gj));
   else
      nanrow = (nanrow | strcmp(gj,''));
   end
end

% Remove rows with NaN anywhere
y(nanrow,:) = [];
n = size(y,1);
ng = length(group);
dfvar = zeros(ng,1);
allgrps = zeros(n, ng);
allnames = cell(ng,1);

% Get arrays describing the groups
for j=1:ng
   gj = group{j};
   gj(nanrow,:) = [];
   group{j} = gj;
   if continuous(j)
      dfvar(j) = 1;
      allgrps(:,j) = gj;
      allnames{j} = {''};
   else
      [gij,gnj] = grp2idx(gj);
      nlevels = size(gnj,1);
      dfvar(j) = nlevels - 1;
      allgrps(:,j) = gij;
      allnames{j} = gnj;
   end
end

% The df and allnames information does not yet reflect nesting.  We will
% fix up allnames for nested variables later, and not use df for them.

varinfo.df = dfvar;
varinfo.allnames = allnames;

n = size(y,1);

% STEP 2:  Create dummy variable arrays for each grouping variable.
% STEP 2a: For type 3 ss, create constraints for each grouping variable.
varinfo.varnames = varnames;
ng       = size(allgrps,2);
vdum     = cell(ng,1);      % dummy terms for variables
vconstr  = cell(ng,1);   % constraints
vlevels  = ones(ng,1);   % levels corresponding to each dummy term
vmeans   = zeros(ng,1);   % variable means, used for continuous factors
allnames = varinfo.allnames;

for j=1:ng
   if continuous(j)
      % For continuous variables, use the variable itself w/o constraints
      vdum{j} = allgrps(:,j);
      vconstr{j} = ones(0,1);
      vmeans(j) = mean(allgrps(:,j));
   elseif isempty(vnested) || ~any(vnested(j,:))
      % Create dummy variable arrays for each grouping variable
      % using a sum-to-zero constraint
      vdum{j} = idummy(allgrps(:,j), 3);
      nlevgrp  = size(vdum{j},2);
      vconstr{j} = ones(1,nlevgrp);  % sum-to-zero constraint for categorical
      vlevels(j) = nlevgrp;
   else
      % Create dummy variable arrays for a nested variable
      nesternums = find(vnested(j,:));
      allvars = allgrps(:,[nesternums, j]);
      [ugrps,~,grpnum] = unique(allvars,'rows');
      vdum{j} = idummy(grpnum,3);
      vconstr{j} = []; % to be computed later
      vlevels(j) = size(vdum{j},2);
      
      % Fix up names for nested variables
      levelnames = varinfo.allnames{j};
      namesj = cell(vlevels(j),1);
      nester1names = allnames{nesternums(1)};
      for rnum = 1:vlevels(j)
          nesterlist = nester1names{ugrps(rnum,1)};
          for k = 2:length(nesternums)
              nesterknames = allnames{nesternums(k)};
              nesterlist = sprintf('%s,%s',nesterlist,...
                                   nesterknames{ugrps(rnum,k)});
          end
          namesj{rnum} = sprintf('%s(%s)',levelnames{ugrps(rnum,end)},...
                                          nesterlist);
      end
      varinfo.allnames{j} = namesj;
   end
end

varinfo.vdum = vdum;
varinfo.vconstr = vconstr;
varinfo.vlevels = vlevels;
varinfo.vmeans = vmeans;

% STEP 3:  Create dummy variable arrays for each term in the model.
nterms = size(termlist,1);
fulltermlist = showtermnesting(termlist,vnested);
tnested = gettermnesting(fulltermlist,vnested,continuous);

[~,sindex] = sortrows(fulltermlist);
terminfo = makedummyterms(sindex, termlist, varinfo, ...
                          continuous, tnested);

% STEP 4:  Create the full design matrix
[dmat,cmat,termname] = makedesignmat(terminfo,n);

% STEP 5:  Fit the full model and compute the residual sum of squares
mu  =  mean(y,1);
y   = bsxfun(@minus,y,mu); % two passes to improve accuracy
mu2 = mean(y,1);
mu  = mu + mu2;
y   = bsxfun(@minus,y,mu2);
sst = sum(y.^2);

[ssx,dfx,dmat2,y] = dofit(dmat, y, cmat, sst, mu);
   
sse = sst - ssx;
dfe = n - dfx;

% STEP 6:  Determine which models to compare for testing each term
ssw  = -ones(nterms, size(y,2));      % sum of squares with this term
sswo = ssw;                   % sum of squares without this term
dfw  = ssw(:,1);                   % residual d.f. with this term
dfwo = ssw(:,1);                   % residual d.f. without this term

modw  = termsnotcontained(fulltermlist,continuous);
k     = (sum(modw,2) == nterms);
TnotC = modw;

modw  = logical(modw);
modwo = logical(modw - eye(nterms));   % get model without this term

% STEP 7:  Fit each model, get its residual SS and d.f.
ssw(k,:) = repmat(ssx,[length(k) 1]);        % for full model we already know the results
dfw(k) = dfx;

% Fit each model separately
dfboth = [dfw; dfwo];

% Consider interactions before their components for type 3 ss, so
% examine the terms in decreasing order of the number factors in the term
sindices = [(1:size(termlist,1)), (1:size(termlist,1))]';

% Fit all required subsets of the full model
fitsubmodels;   % nested function

% STEP 8:  Compute the sum of squares attributed to each term
ssterm = max(0, ssw-sswo);
dfterm = dfw-dfwo;
ssterm(dfterm==0,:) = 0;    % make this exact

% STEP 9:  Compute the mean square for each term
mse = bsxfun(@times,sse,(dfe>0) / max(1, dfe));

t = struct;
t.mse = mse;
t.sst = sst;
t.ssterm = ssterm;
t.dfterm = dfterm;

% ----- nested function
function fitsubmodels
   for j=length(sindices):-1:1
      % Find the next model index to fit
      k = sindices(j);
      
      % Look in unsorted arrays to see if we have already fit this model
      if j>nterms
         k0 = k+nterms;
      else
         k0 = k;
      end
      if dfboth(k0)~=-1
         continue
      end
      
      % Find the model with this index
      if (j > nterms)
         thismod = modwo(k, :);
      else
         thismod = modw(k, :);
      end
      
      % Get the design matrix for this model
      keepterms = find(thismod);
      clist = ismember(termname, [0 keepterms]);
      X = dmat2(:,clist);
      C = cmat(:,clist);

      % Fit this model
      [ssx0,dfx0] = dofit(X, y, C);
      
      % Use these results for each term that requires them
      mod0 = repmat(thismod, nterms, 1);
      k = find(all(modw == mod0, 2));
      ssw(k,:) = repmat(ssx0,[length(k) 1]);
      dfw(k) = dfx0;
      dfboth(k) = 0;
      k = find(all(modwo == mod0, 2));
      sswo(k,:) = repmat(ssx0,[length(k) 1]);
      dfwo(k) = dfx0;
      dfboth(nterms+k) = 0;
   end
end

end % of main function


% --------------------------
function m = termsnotcontained(terms,continuous)
%TERMSNOTCONTAINED Creates a logical matrix indicating term containment
%   m(i,j)==1  iff  t(i) is not contained by t(j)

% For starters, omit continuous variables from this test
contnum = find(continuous);
t = terms;
t(:,contnum) = 0;

% The main test:  no overlap between terms included in t(i) and not in t(j)
m = (t*~t') > 0;

% Now consider continuous variables (X)
nterms = size(terms,1);
for j=1:length(contnum)
   v = terms(:,contnum(j));
   vv = repmat(v,1,nterms);
   % A cannot be contained in B if they don't match on X
   m = m | vv~=vv';
end

% Set diagonals to 1 because we want proper containment
m(1:(nterms+1):end) = 1;
end

% --------------------------
function [ssx,dfx,dmat2,y2,stats] = dofit(dmat,y,cmat,sst,mu)
%DOFIT Do constrained least squares fit and reduce data for subsequent fits

% Find the null space of the constraints matrix
[Qc,Rc,~] = qr(cmat');
pc = Rrank(Rc);
Qc0 = Qc(:,pc+1:end);

% Do qr decomposition on design matrix projected to null space
Dproj = dmat*Qc0;
[Qd,Rd,Ed] = qr(Dproj,0);
Dproj = []; % no longer needed but potentially big
dfx = Rrank(Rd);
Qd = Qd(:,1:dfx);
Rd = Rd(1:dfx,1:dfx);

% Fit y to design matrix in null space
y2 = Qd' * y;            % rotate y into that space
zq = Rd \ y2;            % fit rotated y to projected design matrix

z = zeros(length(Ed),size(y,2)); % coefficient vector extend to full size ...
z(Ed(1:dfx),1:size(y,2)) = zq;       % ... and in correct order for null space

b = Qc0 * z;             % coefficients back in full space
ssx = sum(y2.^2);        % sum of squares explained by fit

% Return reduced design matrix if requested
if nargout>2
   dmat2 = Qd' * dmat;
end

% Prepare for multiple comparisons if requested
if (nargout >= 3)
   stats.source = 'anovan';
   
   % Get residuals
   % yhat = Dproj * z;
   yhat = Qd * Rd * zq;
   stats.resid = y - yhat;
   
   % Calculate coefficients, then adjust for previously removed mean
   t = b;
   if ~isempty(t)
      t(1,:) = t(1,:) + mu;            % undo mean adjustment before storing
   end
   stats.coeffs = t;

   RR = zeros(dfx,length(Ed));
   RR(:,Ed(1:dfx)) = Rd;
   stats.Rtr = RR';
   [~,rr,ee] = qr(dmat,0);
   pp = Rrank(rr);
   rowbasis = zeros(pp,size(rr,2));
   rowbasis(:,ee) = rr(1:pp,:);
   stats.rowbasis = rowbasis;
   stats.dfe = length(y) - dfx;
   stats.mse = (sst - ssx) / max(1, stats.dfe);
   stats.nullproject = Qc0;
end
end

% ---------------------
function p = Rrank(R,tol)
%RRANK Compute rank of R factor after qr decomposition
if (min(size(R))==1)
   d = abs(R(1,1));
else
   d = abs(diag(R));
end
if nargin<2
   tol = 100 * eps(class(d)) * max(size(R));
end
p = sum(d > tol*max(d));
end

% --------------------------------
function terminfo = makedummyterms(sindex,termlist,varinfo,continuous,tnested)

ncols = 1;
[nterms,nfactors] = size(termlist);

termdum = cell(nterms, 1);      % dummy variable = design matrix cols
termconstr = cell(nterms,1);    % constraints to make each term well defined
levelcodes = cell(nterms, 1);   % codes for levels of each M row
tnames = cell(nterms, 1);       % name for entire term, e.g. A*B
dfterm0 = zeros(nterms, 1);     % nominal d.f. for each term
termvars = cell(nterms, 1);     % list of vars in each term
termlength = zeros(size(sindex));% length of each term (number of columns)

auxtermlist = zeros(0,nfactors);  % may need to create temp terms for nesting
auxtermdum = cell(0,1);           % may need their dummy variables

isnested = false;

% For each term,
for j=1:nterms
   % Get dummy columns, var list, term name, etc
   sj = sindex(j);
   tm = termlist(sj,:);
   if ~isempty(tnested)
      isnested = any(tnested(sj,:));
   end
   [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,...
                                         termlist,termdum,termconstr,tnames,dfterm0);

   % Store this term's information
   k = size(tdum, 2);
   termlength(sindex(j),1) = k;
   ncols = ncols + k;
   termdum{sj} = tdum;
   termvars{sj} = vars;
   levelcodes{sj} = fliplr(fullfact(varinfo.vlevels(vars(end:-1:1))));
   tnames{sj,1} = tn;
   dfterm0(sj) = df0;

   % For a nested term, figure out the constraint now
   if isnested
      % First get dummy variables for this term and its nesters
      if any(tm(continuous))
         % Ignore any continuous variable contributions to this term
         tmcat = tm;                % this term, cat predictors only
         tmcat(continuous) = 0;
         ntmcat = termlist(tnested(sj,:),:);
         ntmcat(:,continuous) = 0;  % nesting terms, cat predictors only
         [Ydum,auxtermlist,auxtermdum] =findtermdum(tmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
         Ydum = Ydum{1};
         [Xdum,auxtermlist,auxtermdum] =findtermdum(ntmcat,termlist,termdum,...
                              auxtermlist,auxtermdum,varinfo);
      else
         % No continuous variables, so work with the whole term
         Ydum = tdum;
         Xdum = termdum(tnested(sj,:)>0);
      end

      % Constrain the part of this term in its nesters to be 0
      tconstr = gettermconstr(Ydum, Xdum);
   end
   termconstr{sj} = tconstr;
end
tnames{length(tnames)+1,1} = getString(message('stats:anovan:Error'));

% Package up term info into a structure
terminfo.termdum = termdum;
terminfo.termconstr = termconstr;
terminfo.levelcodes = levelcodes;
terminfo.tnames = tnames;
terminfo.dfterm0 = dfterm0;
terminfo.termvars = termvars;
terminfo.termlength = termlength;
end

% -----------------------------
function [dmat,cmat,termname] = makedesignmat(terminfo,n)

termlength = terminfo.termlength;
ncols = sum(termlength);

nconstr = sum(cellfun('size',terminfo.termconstr,1));
dmat = ones(n, ncols+1);      % to hold design matrix
cmat = zeros(nconstr,ncols);  % to hold constraints matrix
cbase = 0;                    % base from which to fill in cmat
termname = zeros(ncols,1);
termstart = cumsum([2; termlength(1:end-1)]);
termend = termstart + termlength - 1;
for j=1:length(termlength)
   clist = termstart(j):termend(j);
   dmat(:, clist) = terminfo.termdum{j};
   C = terminfo.termconstr{j};
   nC = size(C,1);
   cmat(cbase+1:cbase+nC,clist) = C;
   termname(clist) = j;
   cbase = cbase + nC;
end
end


% -----------------------------
function tnested = gettermnesting(fulltermlist,vnested,continuous)
% Create matrix with (i,j) indicating if term i is nested in term j

% No nested variables implies no nested terms
if isempty(vnested)
    tnested = [];
    return
end

% Work with categorical and continuous terms separately
ctermlist = fulltermlist(:,continuous);
fulltermlist(:,continuous) = 0;

nterms = size(fulltermlist,1);
tnested = zeros(nterms);

% If A(B), then any term containing A may be nested in one containing B
[nestee,nester] = find(vnested);
for j=1:length(nester)
    hasnestee = (fulltermlist(:,nestee(j))>0);
    hasnester = (fulltermlist(:,nester(j))>0);
    tnested(hasnestee,hasnester) = 1;
end

% But terms are not nested within themselves
tnested(1:nterms+1:end) = 0;

% And there is no nesting relationship if the continuous vars don't match
for j=1:size(ctermlist,2)
    c = repmat(ctermlist(:,j),1,nterms);
    tnested = tnested & (c == c');
end

% And the full representation for a nesting term must be a subset of the
% full representation for a nested term
[nestee,nester] = find(tnested);
for j=1:length(nester)
    if any(fulltermlist(nestee(j),:) < fulltermlist(nester(j),:))
        tnested(nestee(j),nester(j)) = 0;
    end
end

end

% ----------------------------
function constr = gettermconstr(tdum,nesterdum)

% Create X matrix containing design matrix columns of nester terms
if iscell(nesterdum)
    Xnester = cat(2,nesterdum{:});
else
    Xnester = nesterdum;
end
nnested = size(Xnester,2);

% Find unique combinations of nested and nesting terms
Xboth = [Xnester tdum];
Uboth = unique(Xboth,'rows');

% Find constraints to force the part of the nested terms in the column
% space of the nesting terms to be zero
[Q,R,~] = qr(Uboth(:,1:nnested));
rnk = Rrank(R);
Q = Q(:,1:rnk);
constr = Q*(Q'*Uboth(:,nnested+1:end));

% Remove redundant constraints
[~,R,E] = qr(constr',0);
rnk = Rrank(R);
constr = constr(E(1:rnk),:);
end

% ----------------------------
function fulltermlist = showtermnesting(termlist,vnested)

fulltermlist = termlist;

[nestee,nester] = find(vnested);
for j=1:length(nester)
    t = fulltermlist(:,nestee(j))>0;
    fulltermlist(t,nester(j)) = 1;
end
end


% ----------------------------
function [tdum,vars,tn,df0,tconstr] = maketerm(tm,isnested,varinfo,j,sindex,termlist,...
                                             termdum,termconstr,tnames,dfterm0)
% Make term info such as dummy vars, name, constraints

% Loop over elements of the term
df0 = 1;
tdum = [];         % empty term so far
tconstr = 1;       % empty constraints so far
tn = '';           % blank name so far
vars = find(tm);   % list of variables making up terms
pwrs = tm(vars);   % and their exponents
tm = (tm>0);       % and a boolean mask for them
for varidx = 1:length(vars)
   % Process each variable participating in this term
   varnum = vars(varidx);          % variable name
   thispwr = pwrs(1);              % power of this variable
   tm(varnum) = 0;                 % term without this variable
   pwrs(1) = [];                   % powers without this variable
   df0 = df0 * varinfo.df(varnum); % d.f. so far

   % Combine its dummy variable with the part computed so far
   G = varinfo.vdum{varnum};       % dummy vars for this grouping var
   thisname = varinfo.varnames{varnum};
   if thispwr>1
       G = G .^ thispwr;
       thisname = sprintf('%s^%d',thisname,thispwr);
   end
   tdum = termcross(G,tdum);   % combine G into term dummy vars

   % Construct the term name and constraints matrix
   if nargout>1
      if (isempty(tn))
         tn = thisname;
         if ~isnested
            tconstr = varinfo.vconstr{varnum};
         end
      else
         tn = [tn '*' thisname];
         if ~isnested
            tconstr = [kron(varinfo.vconstr{varnum},eye(size(tconstr,2)));
                       kron(eye(length(varinfo.vconstr{varnum})),tconstr)];
         end
      end
   end

   if varidx<length(vars) && j>1
      % If the rest of this term is computed, take advantage of that
      prevterms = termlist(sindex(1:j-1),:);
      oldtm = find(all(prevterms(:,tm) == repmat(pwrs,size(prevterms,1),1),2));
      oldtm = oldtm(~any(prevterms(oldtm,~tm),2));
      if ~isempty(oldtm)
         k = sindex(oldtm(1));
         tdum = termcross(termdum{k}, tdum);
         if nargout>1
            if ~isnested
               oconstr = termconstr{k};
               tconstr = [kron(tconstr,              eye(size(oconstr,2)));
                          kron(eye(size(tconstr,2)), oconstr)];
            end
            tn = [tn '*' tnames{k}];
            df0 = df0 * dfterm0(k);
         end
         break;
      end
   end
end
if (isempty(tn)), tn = getString(message('stats:anovan:Constant')); end
end

% ----------------------------
function [dum,auxtermlist,auxtermdum] =findtermdum(tms,termlist,termdum,...
                   auxtermlist,auxtermdum,varinfo)

 
nterms = size(tms,1);
dum = cell(nterms,1);
for j=1:nterms
    tm = tms(j,:);

    % Try to find this term among the terms in the model
    oldtm = find(all(termlist == repmat(tm,size(termlist,1),1),2));
    if ~isempty(oldtm)
        dum{j} = termdum{oldtm(1)};
        continue
    end
        
    % Try to find this term among the auxiliary terms
    oldtm = find(all(auxtermlist == repmat(tm,size(auxtermlist,1),1),2));
    if ~isempty(oldtm)
        dum{j} = auxtermdum{oldtm(1)};
        continue
    end
    
    % Create an auxiliary term
    dum{j} = maketerm(tm,[],varinfo,0);

    auxtermdum{end+1} = dum{j};
    auxtermlist(end+1,:) = tm;
end
end

function d = idummy(x, method)
% function d = idummy(x, method)

n = length(x);
g = max(x);
ncols = g - (method ~= 3);
d = repmat(0, n, ncols);

if (g > 1 || method==3)
   % Fill in -1 for the first level
   if (method == 1)
      d((x == 1),:) = -1;
   end
   
   % Fill in 1 in the appropriate column for other levels
   m3 = (method == 3);
   for j=(2-m3):g
      d((x == j),j-1+m3) = 1;
   end
end
end

function ab = termcross(a,b)

if (isempty(a)), ab = b; return, end
if (isempty(b)), ab = a; return, end

na = size(a,2);
nb = size(b,2);
acols = repmat((1:na), 1, nb);
bcols = reshape(repmat((1:nb), na, 1), 1, na*nb);
ab = a(:,acols) .* b(:,bcols);
end

function p = fpval(x,df1,df2)
xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);
end