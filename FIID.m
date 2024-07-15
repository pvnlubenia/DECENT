% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    FIID                                                                   %
%                                                                           %
%                                                                           %
% OUTPUT: Returns the finest nontrivial incidence independent decomposition %
%    of a chemical reaction network (CRN), if it exists. If no such         %
%    decomposition exists, a message appears saying so. It also shows a     %
%    table summarizing the network numbers of the entire network and of     %
%    the subnetworks. The sum of the subnetworks' network numbers, which    %
%    may provide additional insights, also appear at the end of the table.  %
%    The output variables 'model', 'I_a', 'G', and 'P' allow the user to    %
%    view the following, respectively:                                      %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%       - Incidence matrix of the network                                   %
%       - Undirected graph of I_a                                           %
%       - Partitions representing the decomposition of the reactions        %
%       - Table of network numbers                                          %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% References:                                                               %
%    [1] Arceo C, Jose E, Lao A, Mendoza E (2017) Reactant subspaces and    %
%           kinetics of chemical reaction networks. J Math Chem             %
%           56(5):395–422. https://doi.org/10.1007/s10910-017-0809-x        %
%    [2] Hernandez B, Amistas D, De la Cruz R, Fontanil L, de los Reyes V   %
%           A, Mendoza E (2022) Independent, incidence independent and      %
%           weakly reversible decompositions of chemical reaction networks. %
%           MATCH Commun Math Comput Chem, 87(2):367-396.                   %
%           https://doi.org/10.46793/match.87-2.367H                        %
%    [3] Hernandez B, De la Cruz R (2021) Independent decompositions of     %
%           chemical reaction networks. Bull Math Biol 83(76):1–23.         %
%           https://doi.org/10.1007/s11538-021-00906-3                      %
%    [4] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical        %
%           reaction network theory. Bioinform 25(21):2853–2854.            %
%           https://doi.org/10.1093/bioinformatics/btp513                   %
%                                                                           %
% Created: 14 July 2024                                                     %
% Last Modified: 15 July 2024                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [model, I_a, G, P, results] = FIID(model)
    
%
% Step 1: Determine the finest independent decomposition of the network
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);

% Count the number of species
m = numel(model.species);

% Initialize the matrix of reactant complexes
reactant_complexes = [ ];

% Initialize the matrix of product complexes
product_complexes = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
  
    % Initialize the vector for the reaction's reactant complex
    reactant_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end);
    
    % If the reaction is reversible
    if model.reaction(i).reversible
      
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complexes(:, end+1) = product_complexes(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complexes(:, end+1) = reactant_complexes(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Count the total number of reactions
r = size(N, 2);

% Get just the unique complexes
% index(i) is the index in all_complex of the reactant complex in reaction i
[all_complex, ~, index] = unique([reactant_complexes product_complexes]', 'rows');

% Construct the matrix of complexes
all_complex = all_complex';

% Count the number of complexes
n = size(all_complex, 2);

% Initialize matrix (complexes x total reactions) for the reacts_in relation
% This is the incidence matrix I_a
reacts_in = zeros(n, r);

% Fill out the entries of the matrices
for i = 1:r
    
    % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
    reacts_in(index(i), i) = -1;
    reacts_in(index(i+r), i) = 1;
end

% Construct the incidence matrix
I_a = reacts_in;

% Get the transpose of N: Each row now represents the reaction vector a reaction
R = I_a';

% Write R in reduced row echelon form: the transpose of R is used so 'basis_reaction_num' will give the pivot rows of R
%    - 'basis_reaction_num' gives the row numbers of R which form a basis for the rowspace of R
[~, basis_reaction_num] = rref(R');

% Form the basis for the rowspace of R
basis = R(basis_reaction_num, :);

% Construct the vertex set of undirected graph G
% Initialize an undirected graph G
G = graph();

% Add vertices to G: these are the reaction vectors that form a basis for the rowspace of R
for i = 1:numel(basis_reaction_num)
    
    % Use the reaction number as label for each vertex
    G = addnode(G, strcat('R', num2str(basis_reaction_num(i))));
end

% Write the nonbasis reaction vectors as a linear combination of the basis vectors
% Initialize matrix of linear combinations
linear_combo = zeros(r, numel(basis_reaction_num));

% Do this for the nonbasis reactions vectors
for i = 1:r
    if ~ismember(i, basis_reaction_num)
      
      % This gives the coefficients of the linear combinations
      % The basis vectors will have a row of zeros
      linear_combo(i, :) = basis'\R(i, :)';
    end
end

% Round off values to nearest whole number to avoid round off errors
linear_combo = round(linear_combo);

% Construct the edge set of undirected graph G
% Get the reactions that are linear combinations of at least 2 basis reactions
% These are the reactions where we'll get the edges
get_edges = find(sum(abs(linear_combo), 2) > 1);
    
% Initialize an array for sets of vertices that will form the edges
vertex_set = { };
 
% Identify which vertices form edges in each reaction: get those with non-zero coefficients in the linear combinations
for i = 1:numel(get_edges)
    vertex_set{i} = find(linear_combo(get_edges(i), :) ~= 0);
end

% Initialize the edge set
edges = [ ];

% Get all possible combinations (not permutations) of the reactions involved in the linear combinations
for i = 1:numel(vertex_set)
    edges = [edges; nchoosek(vertex_set{i}, 2)];
end

% Get just the unique edges
edges = unique(edges, 'rows');

% Add these edges to graph G
for i = 1:size(edges, 1)
    G = addedge(G, strcat('R', num2str(basis_reaction_num(edges(i, 1)))), strcat('R', num2str(basis_reaction_num(edges(i, 2)))));
end

% Check if G is connected, i.e., has only one connected component
% Determine to which component each vertex belongs to
component_numbers = conncomp(G);

% Determine the number of connected components of G: this is the number of partitions R will be decomposed to
num_components = max(component_numbers);

% For the case of only one connected component
if num_components == 1
    P = [ ];
    disp([model.id ' has no nontrivial independent decomposition.']);
    
    % 'return' exits the function; we don't need to continue the code
    % If we wanted to just get out of the loop, we use 'break'
    return
end

% If G is NOT connected, form the partitions of R
% Initialize the list of partitions
P = cell(1, num_components);

% Basis vectors: assign them first into their respective partition based on their component number
for i = 1:numel(component_numbers)
    P{component_numbers(i)}(end+1) = basis_reaction_num(i);
end

% Nonbasis vectors: they go to the same partition as the basis vectors that form their linear combination
for i = 1:numel(P)
    for j = 1:numel(P{i})
        
        % Get the column number representing the basis vectors in 'linear_combo'
        col = find(basis_reaction_num == P{i}(j));
        
        % Check which reactions used a particular basis vector and assign them to their respective partition
        P{i} = [P{i} find(linear_combo(:, col) ~= 0)'];
    end
end

% Get only unique elements in each partition
for i = 1:numel(P)
    P{i} = unique(P{i});
end

% Check if all the reactions are in the partitions
%    - If not all reactions are partitions, usually it's because of the computation of linear combinations
% If some reactions are missing, then redo the end of Step a up to Step d (labeled here but not above)
if length(cell2mat(P)) ~= size(R, 1)

    % Step a end: Do not round off the coefficients of the linear combinations
    linear_combination = linear_combo;

    % Step b
    get_edges = find(sum(abs(linear_combination), 2) > 1);
    vertex_set = { };
    for i = 1:numel(get_edges)
        vertex_set{i} = find(linear_combination(get_edges(i), :) ~= 0);
    end
    edges = [ ];
    for i = 1:numel(vertex_set)
        edges = [edges; nchoosek(vertex_set{i}, 2)];
    end
    edges = unique(edges, 'rows');
    for i = 1:size(edges, 1)
        G = addedge(G, strcat('R', num2str(basis_reaction_num(edges(i, 1)))), strcat('R', num2str(basis_reaction_num(edges(i, 2)))));
    end
    
    % Step c
    component_numbers = conncomp(G);
    num_components = max(component_numbers);
    if num_components == 1
        P = [ ];
        disp([model.id ' has no nontrivial independent decomposition.']);
        return
    end
    
    % Step d
    P = cell(1, num_components);
    for i = 1:numel(component_numbers)
        P{component_numbers(i)}(end+1) = basis_reaction_num(i);
    end
    for i = 1:numel(P)
        for j = 1:numel(P{i})
            col = find(basis_reaction_num == P{i}(j));
            P{i} = [P{i} find(linear_combination(:, col) ~= 0)'];
        end
    end
    for i = 1:numel(P)
        P{i} = unique(P{i});
    end
end

% Display the independent decomposition
% Use 'fprintf' instead of 'disp' to interpret '\n' as 'newline'
fprintf('\nFinest incidence independent decomposition of %s:\n\n', model.id)
for i = 1:numel(P)
    subnetwork = sprintf('R%d, ', P{i});
    subnetwork(end-1:end) = [ ]; % To clean the trailing comma and space at the end of the list
    fprintf('N%d: %s \n', i, subnetwork);
end
fprintf('\n\n')



%
% Step 2: Generate the network numbers of the network
%

% Generate the list of network numbers for the whole network
[model, characteristics, notation, network] = networkNumbers(model);

% Prepare a list of separators between the numbers for the entire network and the subnetworks
separator = {'|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'; '|'};

% Convert the list to a table
headers = {'Characteristics', 'Notation', 'Network', '|'};
results = cell2table([characteristics, notation, network, separator], 'VariableNames', headers);
results = convertvars(results, results.Properties.VariableNames, 'categorical');

% Create a vector of reaction numbers for the total number of reactions
reac_num = [ ];
for i = 1:numel(model.reaction)
    if model.reaction(i).reversible == 0
        reac_num(end+1) = i;
    else
        reac_num(end+1) = i;
        reac_num(end+1) = i;
    end
end



%
% Step 3: Generate the network numbers of the subnetworks
%

% Initialize list of network numbers
sums = { };

% Go through each subnetwork
for i = 1:numel(P)

    % Get the reactions for the subnetwork
    model_P(i).id = [model.id ' - Subnetwork ' num2str(i)];
    model_P(i).species = { };
    reac_P = unique(reac_num(P{i}));
    for j = 1:numel(reac_P)
        model_P(i).reaction(j) = model.reaction(reac_P(j));
    end

    % Generate the network numbers for the subnetwork
    [~, ~, ~, network_] = networkNumbers(model_P(i));

    % Add the network numbers to the list
    sums(:, end+1) = network_;

    % Append subnetwork network numbers to the original table of network numbers
    header_ = {['N', num2str(i)]};
    results = addvars(results, network_, 'NewVariableNames', header_);
end

% Get the sum of each network number of the subnetworks
sums = sum(cell2mat(sums), 2);

% Append this sum to the summary table
header_ = {'Sum'};
results = addvars(results, sums, 'NewVariableNames', header_);

% Fix the appearance of the table
results = convertvars(results, results.Properties.VariableNames, 'string');
results = convertvars(results, results.Properties.VariableNames, 'categorical');

% Display the table
fprintf('Network numbers:\n\n')
disp(results)

end










% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                 % %
% % The following function is used in the algorithm % %
% %                                                 % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% networkNumbers                                                      %
%                                                                     %
%    - Purpose: To determine the network numbers of chemical reaction %
%          network                                                    %
%    - Input                                                          %
%         - model: complete or incomplete structure                   %
%    - Outputs                                                        %
%         - model: complete structure                                 %
%         - characteristics: list of network numbers                  %
%         - notation: list of variables of the network numbers        %
%         - network: list of values of the network numbers            %
%    - Used in                                                        %
%         - FID_check (Step 2)                                        %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model, characteristics, notation, network] = networkNumbers(model)
    
%
% Create a list of all species indicated in the reactions
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);



%
% Species
%

% Count the number of species
m = numel(model.species);



%
% Complexes
%

% Initialize the matrix of reactant complexes
reactant_complexes = [ ];

% Initialize the matrix of product complexes
product_complexes = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
    
    % Initialize the vector for the reaction's reactant complex
    reactant_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complexes(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end); % N %
    
    % If the reaction is reversible
    if model.reaction(i).reversible
        
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complexes(:, end+1) = product_complexes(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complexes(:, end+1) = reactant_complexes(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Get just the unique complexes
% index(i) is the index in all_complex of the reactant complex in reaction i
[all_complex, ~, index] = unique([reactant_complexes product_complexes]', 'rows');

% Construct the matrix of complexes
all_complex = all_complex';

% Count the number of complexes
n = size(all_complex, 2);



%
% Reactant complexes
%

% Get just the unique reactant complexes
reactant_complexes_unique = unique([reactant_complexes]', 'rows')';

% Count the number of unique reactant complexes
n_r = size(reactant_complexes_unique, 2);



%
% Reversible, irreversible, and total reactions
%

% Count the number of reversible, irreversible, and total reactions
r_rev = 0;
r_irrev = 0;
for i = 1:numel(model.reaction)
    if (model.reaction(i).reversible == 1)
        r_rev = r_rev + 1;
     else
        r_irrev = r_irrev + 1;
    end
end
r = r_irrev + 2*r_rev;



%
% Linkage classes
%

% Initialize a matrix (complexes x complexes) for the reacts_to relation
% This is for testing reversibility of the network
reacts_to = false(n, n);

% Initialize matrix (complexes x total reactions) for the reacts_in relation
% This is the incidence matrix I_a
reacts_in = zeros(n, r);

% Fill out the entries of the matrices
for i = 1:r
    
    % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
    reacts_to(index(i), index(i + r)) = true;
    
    % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
    reacts_in(index(i), i) = -1;
    reacts_in(index(i+r), i) = 1;
end

% Linkage classes
% Count number of connected components of an undirected graph
linkage_class = conncomp(graph(reacts_to | reacts_to'));

% Count the number of linkage classes
l = max(linkage_class);



%
% Strong linkage classes
%

% Check if the network is reversible
is_reversible = isequal(reacts_to, reacts_to');

% Strong linkage classes
% Count number of strongly connected components of an directed graph
if is_reversible
    strong_linkage_class = linkage_class;
else
    % Count number of connected components of a directed graph
    strong_linkage_class = conncomp(digraph(reacts_to));
end

% Count the number of strong linkage classes
sl = max(strong_linkage_class);



%
% Terminal linkage classes
%

% Count the number of terminal strong linkage classes
% Initialize the count
t = 0;
is_nontrivial_terminal_slc = false(sl, 1);
is_terminal_complex = false(n, 1);
for i = 1:sl
    
    % Locate the indexes in Y of the complexes present in strong-linkage class i
    complexes_i = find(strong_linkage_class == i);
    
    % Locate the indexes in Y of the complexes which in some reactions are products of complexes_i
    products_of_complexes_i = find(any(reacts_to(complexes_i, :), 1));
    
    % Products_of_complexes_i is a subset of complexes_i, so the strong-linkage class i is terminal
    if all(ismember(products_of_complexes_i, complexes_i))
        t = t + 1;
        is_terminal_complex(complexes_i) = true;
        if numel(complexes_i) > 1
            is_nontrivial_terminal_slc(i) = true;
        end
    end
end



%
% Rank
%

% Get the rank of the reaction network
% S = Im N
% dim S = dim (Im N) = rank(N)
% Note: We talk of "dimension of a linear transformation" and "rank of a matrix"
s = rank(N);



%
% Reactant rank
%

% Construct the incidence matrix
% We can decompose this into I_a = I_a^+ - I_a^-
I_a = reacts_in;

% Construct I_a^-: for reactant complexes only
I_a_minus = I_a;
I_a_minus(I_a_minus > 0) = 0;
I_a_minus(I_a_minus < 0) = 1;

% Construct N^-: "reactant subspace" matrix
N_minus = all_complex*I_a_minus;

% Get the reactant rank
% R = Im N^-
% dim R = dim (Im N^-) = rank(N^-)
q = rank(N_minus);


%
% Deficiency
%

% Compute the deficiency of the reaction network
delta = n - l - s;



%
% Reactant deficiency
%

% Compute the reactant deficiency
delta_p = n_r - q;



%
% Table columns
%

% List of characteristics
characteristics = {'Species'; 'Complexes'; 'Reactant complexes'; 'Reversible reactions'; 'Irreversible reactions'; 'Reactions'; 'Linkage classes'; 'Strong linkage classes'; 'Terminal linkage classes'; 'Rank'; 'Reactant rank'; 'Deficiency'; 'Reactant deficiency'};

% List of notations
notation = {'m'; 'n'; 'n_r'; 'r_rev'; 'r_irrev'; 'r'; 'l'; 'sl'; 't'; 's'; 'q'; 'delta'; 'delta_p'};

% List of network numbers
network = {m; n; n_r; r_rev; r_irrev; r; l; sl; t; s; q; delta; delta_p};

end
