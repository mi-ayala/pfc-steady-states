%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------EXTRACT COEFFICIENT FROM ACCUMULATED----------------------------
%%%%%%%%%%%%%% Saves the coefficients and pfc_g of a given state in an accumulated file
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Sep 18 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accumulated_SS_file = 'CoefficientFiles/all_SS_Table5.mat';
extract_state_id = 1;
save_to_file = 'CoefficientFiles/table_5/state_row1.mat';

load(accumulated_SS_file);
A = list_SS{extract_state_id}.A
save(save_to_file, 'pfc_g', 'A')
