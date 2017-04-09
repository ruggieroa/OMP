function [sig residue_hist] = my_OMP(N,x,D,K)

% my OMP algorithm, modeled after the methods of:
% http://www.mathworks.com/help/wavelet/ug/matching-pursuit-algorithms.html
%
%
%   x - the measurement taken -> x = D*s
%   s - the signal we wish to recover
%   D - the demodulation matrix (or measurement matrix)
%   K - sparsity

% set the tolerance level to escape the loop
tol = 5e-11;

% define matrix functions
Df = @(x) D*x;
Dt = @(x) D'*x;

% initialize residues
r = x;
Dr = Dt(r);

nAt = size(Dr,1); % num atoms
sAt = size(r,1); % size of atoms

residue_hist = [];

if K > nAt
  error("K is bigger than the number of atoms!");
endif

% norm of the measurement
max_L2 = norm(x);

% intialize everything else to 0, i.e. allocate it
%x = zeros(nAt,1);
unitVec = zeros(N,1);
index_set = zeros(length(x),1);
index_sort = zeros(length(x),1);
A_T = zeros(sAt, length(x));
A_T_not = zeros(sAt, length(x));
s = [];
s = zeros(N, 1);

for j = 1:K
  % find new index and atom
  [dummy,temp] = max(abs(Dr));
  
  index_set(j) = temp;
  index_sort(1:j) = sort(index_set(1:j));
  
  
  % octave matrix is ROW by COL
  
  unitVec(temp)     = 1;
  new_atom                = Df( unitVec );
  unitVec(temp)     = 0;
  %new_atom = D(:,temp);
  

  
  if (1 == 2)
  %
  % slow mode, the '\' operator is not time efficient
  % 
    s_T = A_T_not(:, 1:j)\x;
    
    s(index_set(1:j)) = s_T;
    
    r = x - A_T_not(:, 1:j)*s_T;
  else
  %%
  %% fast mode
  %%
  %%
  for i=1:j
    % MGS to make new_atom orthogonal in A_T
    new_atom = new_atom - (A_T(:,i)'*new_atom)*A_T(:,i);
  endfor
  
  new_atom = new_atom/norm(new_atom); % normalization (else norm(A_T) tends to infinity)
  A_T(:,j) = new_atom;  % add the new atom to the orthogonal projection operator
  
  % Least Squares solve
  % A_T is already orthogonal so A_T' = A_T^-1 (faster than inverting)
  % else the proper way is s_T = A_T(:,1:j)\x;
  s_T = A_T(:,1:j)'*x;
  s(index_set(1:j),1) = s_T;
  
  % update the residue
  r = x - A_T(:,1:j)*s_T;
  %%
  %%
  %% end fast mode
  %%
  endif
  
  %disp(['j=' num2str(j) ' and size = ' num2str(size(index_set))]);
  
  normR = norm(r);
  
  residue_hist(j) = normR;
  
  %disp(['Current norm(r) = ' num2str(normR)]);
  fflush(stdout);
  
  if (normR < tol)
    disp(["Final norm = " num2str(normR)]);
    disp(['K = ' num2str(K)])
    disp(['# of Iterations = ' num2str(i)])
    fflush(stdout);
    break;
  endif
  
  Dr = Dt(r);
  
endfor



%disp(['The original norm of the signal was  ' num2str(max_L2)]);

sig = s/norm(s);
