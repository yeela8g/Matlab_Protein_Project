%main script for managing the program running.

%c
clear;
close all;

%load pdb file
pdb_data = select_and_read_pdb();

%validate heterogen atom names is a row array, instead of 2D matrix.
pdb_data.HeterogenName = arrayfun(@toRow, pdb_data.HeterogenName);

%D - main program loop - user can keep choosing ligand and chain as pleased.
while true

%get from user the desired chain and ligand
[chain_id, ligand_id] = select_chain_ligand(pdb_data);

%get the protein id code
idcode = pdb_data.Header.idCode;

%graph title
text = ['file ' idcode ' chain ' chain_id ' hetero ' ligand_id];

%create struct of user choices.
user_info = struct('idcode', [], 'ligand_id', [], 'chain_id', [], 'foreign_atoms', [], ...
    'protein_atoms', [], 'start_index_atoms', [], 'end_index_atoms', [],'resNames', []);

%E) update the user info struct
user_info.ligand_id = ligand_id;
user_info.chain_id = chain_id;
user_info.idcode = idcode;

% Get only the Ligand's Atoms that are as the same of the chain_id (or all ligands)
foreign_atoms = retrieve_data_for_chain(pdb_data.Model.HeterogenAtom, chain_id);

% Getting Atoms of the specific LIgand (first_heterogen_name.hetID)
user_info.foreign_atoms = retrieve_atoms(foreign_atoms, ligand_id);

% Get only the protein's Atoms the are as the same of the chain_id
user_info.protein_atoms = retrieve_data_for_chain(pdb_data.Model.Atom, chain_id);

%Get the start and end atoms of each amino acid + update in user info struct
[start_index_atoms, end_index_atoms] = start_end_atoms(user_info.protein_atoms);
user_info.start_index_atoms = start_index_atoms;
user_info.end_index_atoms = end_index_atoms;

%Get the start atoms on the chosen chain
binary = zeros(1, length([user_info.protein_atoms.AtomSerNo]));
i = 1;
for x = [user_info.protein_atoms.AtomSerNo]
    a = x == [user_info.protein_atoms(user_info.start_index_atoms).AtomSerNo];%check if x is an serNum of amino acid's start atom
    if any(a)
        binary(i) = 1; %mark this atom as a start atom
    end
    i = i + 1;
end
protein_atoms_subset = user_info.protein_atoms(logical(binary)); %convery binary from double to logic

%save the names of all amino acid
user_info.resNames = strtrim({protein_atoms_subset.resName});

%get backbone atoms coordinates
backbone_coordinates = retrieve_backbone_coordinates(user_info);
spatial_presentation = figure; %open new graphic window

% Use the coordinates variable to plot a 3D scatter graph 
scatter3(backbone_coordinates(:,1), backbone_coordinates(:,2), backbone_coordinates(:,3),"cyan","filled",".");
line(backbone_coordinates(:,1), backbone_coordinates(:,2), backbone_coordinates(:,3))
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
axis;
title('3D Scatter Plot of Atoms');
hold on; % Keep the plot window open

%get the coordinates of the foreign atoms
foreign_coordinates = extract_coordinates(user_info.foreign_atoms);
%Add to the spinal graph the atoms of the foreign material in unconnected black markers
scatter3(foreign_coordinates(:, 1), foreign_coordinates(:, 2), foreign_coordinates(:, 3), 'black', 'MarkerEdgeColor', 'black');

dis_metrix = calc_distances(user_info);
minmax_dis = struct("min_distance",[min(dis_metrix(:))], "max_distance",[max(dis_metrix(:))]);


% calculate the edges of the bins
edges = linspace(minmax_dis.min_distance, minmax_dis.max_distance, 49);
% plot the histogram in a new window and save the window ID in a variable
histograma  = figure;
histogram(dis_metrix(:), edges);

% add axis labels and a title
xlabel('Distance');
ylabel('counts');
title(text);
heatm = figure('visible','off');
minmax_dis2 = struct("min_distance",[min(dis_metrix(:))], "max_distance",[max(dis_metrix(:))]);
while true
  figure(histograma);
  minmax_dis2 = receive_user_input(minmax_dis2);
  figure(heatm);
  %clf;
  graphical_display_of_the_mutual_distances(user_info, dis_metrix, minmax_dis2)

  %ask the user for exiting or repeating the loop
  answer1 = questdlg("would you like to change the minimum or the maximum distance?");
  switch answer1
      case 'Yes'
          clf
      case 'No'
          break;
      case 'Cancel'
          break;
      case ' '
          break;
  end
end


%finds the amino acids that close enough to the ligand
amino_acids = find_closest_aminoAcids(user_info, dis_metrix, minmax_dis2.min_distance);
figure(spatial_presentation);
% Add a title to the graph
title(text);
%struct of colors
colors = cool(length(amino_acids));
% Cell array for curve titles
curveTitles = {};
% Cell for ID-chain backbone characters
curveTitles{1} = {'backbone ' user_info.chain_id};

% Cell for foreign material code
curveTitles{2} = {user_info.ligand_id};

t = 1;
%show the atoms on the figure
for i = amino_acids
    %create a list of the indexes between the first and the last atom in the aa
    integer_aa = user_info.start_index_atoms(i):user_info.end_index_atoms(i);
    % extract all the atoms of the aa
    amino_acid_atoms = user_info.protein_atoms(integer_aa);
    amino_coordinates = extract_coordinates(amino_acid_atoms);
    %takes a color
    color = colors(t,:);
    %show the aa on the spatial presentation 
    scatter3(amino_coordinates(:,1), amino_coordinates(:,2), amino_coordinates(:,3), 'filled', 'Marker', 'o', 'MarkerFaceColor', color);
    aminoAcidIndex = amino_acid_atoms(1); %to find the name
    aminoAcidName = aminoAcidIndex.resName; %get the name
    curveTitles{t+2} = {aminoAcidName};
    t = t+1;
end

curveTitles = curveTitles(~cellfun(@isempty,curveTitles));
legendLabels = [curveTitles{:}];
legend(legendLabels);


%ask the user for exiting or repeating the loop
answer = questdlg("would you like to choose another chain and ligand from this file?");
switch answer
    case 'Yes'
        close all;
    case 'No'
        break;
    case 'Cancel'
        break;
    case ' '
        break;
end

end %of main script loop
