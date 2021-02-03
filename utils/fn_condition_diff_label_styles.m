function [labels, names, diff_pairs, colors, line_styles] = fn_condition_diff_label_styles(factor_name)
%% Converts a set of conditions into difference pairs (labels, plotting colors/styles)

switch factor_name
    %----------------------------------------------------------------------
    % Target Time Contrasts:
    %   NOTE: RPE feature matching is prioritized over outcome features
    %   since those features better explain the neural data. Outcome
    %   feature mismatcvhes are noted in side comments.
    case 'Neg-Pos'
        error('Subtract negative from positive instead!');
    case 'RewP'
        % Valence contrast: negative - positive
        labels   = {{'AllPos','AllNeg'}}; % Positive: 'HdWn','EzWn','HdSu'; Negative: 'EzLs','HdLs','EzSu'
        names    = {'Pos-Neg'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = zeros([numel(labels) 2]);
        for pair_ix = 1:numel(labels)
            diff_pairs(pair_ix,1) = find(strcmp(cond_lab,labels{pair_ix}(1)));
            diff_pairs(pair_ix,2) = find(strcmp(cond_lab,labels{pair_ix}(2)));
        end
        colors = {[0 0 0]};
        line_styles = {'-'};
    case 'Pos-Neg'
        % Valence contrast: negative - positive
        labels   = {{'HdWn','EzLs'},... % lo prob, lg mag (outcome mag = lg)
                    {'EzWn','HdLs'},... % hi prob, sm mag (outcome mag = lg) 
                    {'HdSu','EzSu'}};   % lo prob, sm mag (outcome mag = sm)
        names    = {'HW-EL: loPrb,lgMag',...
                    'EW-HL: hiPrb,lgMag',...
                    'HS-ES: loPrb,smMag'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = zeros([numel(labels) 2]);
        for pair_ix = 1:numel(labels)
            diff_pairs(pair_ix,1) = find(strcmp(cond_lab,labels{pair_ix}(1)));
            diff_pairs(pair_ix,2) = find(strcmp(cond_lab,labels{pair_ix}(2)));
        end
        colors = {[0 0 0], [0.2 0.2 0.2], [0.4 0.4 0.4]}; %Blue/Red: {[31 120 180]./256, [227 26 28]./256};
        line_styles = {'-', '-.', ':'};
    case 'Large-Small'
        % Magnitude contrast: large - small
        labels   = {{'EzLs','EzSu'},... % neg val, lo prob (outcome val = neg - neu)
                    {'HdWn','HdSu'}};   % pos val, lo prob (outcome val = neu - pos)
        names    = {'EL-ES: negVal,loPrb',...
                    'HW-HS: posVal,loPrb'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = zeros([numel(labels) 2]);
        for pair_ix = 1:numel(labels)
            diff_pairs(pair_ix,1) = find(strcmp(cond_lab,labels{pair_ix}(1)));
            diff_pairs(pair_ix,2) = find(strcmp(cond_lab,labels{pair_ix}(2)));
        end
        colors = {[0 0 0], [0.3 0.3 0.3]};
        line_styles = {'-', ':'};
    case 'Unlik-Lik'
        % Probability contrast: unlikely - likely
        labels   = {{'EzSu','HdLs'},... % neg val, sm mag (outcome mag = lo - hi)
                    {'HdSu','EzWn'}};   % pos val, sm mag (outcome mag = lo - hi)
        names    = {'ES-HL: negVal,loMag',...
                    'HS-EW: posVal,loMag'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
        diff_pairs = zeros([numel(labels) 2]);
        for pair_ix = 1:numel(labels)
            diff_pairs(pair_ix,1) = find(strcmp(cond_lab,labels{pair_ix}(1)));
            diff_pairs(pair_ix,2) = find(strcmp(cond_lab,labels{pair_ix}(2)));
        end
        colors = {[0 0 0], [0.3 0.3 0.3]};
        line_styles = {'-', ':'};
    
    %----------------------------------------------------------------------
    % Oddball Contrasts:
    case 'TarStd'
        labels  = {'Tar-Std'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles('Odd');
        diff_pairs= {[find(strcmp(cond_lab,'Tar')) find(strcmp(cond_lab,'Std'))]};
        colors = {[31 120 180]./256};
        line_styles = {'-'};
    case 'OddStd'
        labels  = {'Odd-Std'};
        [cond_lab,~,~,~,~] = fn_condition_label_styles('Odd');
        diff_pairs= {[find(strcmp(cond_lab,'Odd')) find(strcmp(cond_lab,'Std'))]};
        colors = {[31 120 180]./256};
        line_styles = {'-'};
    otherwise
        error(['Difference wave not defined for factor ' factor_name]);
end

%% OLD Contrasts:
%     case 'DifOutUE'
%         labels   = {'HW-EW', 'EL-HL'};
%         [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
%         diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
%             [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]};
%         colors = {[31 120 180]./256, [227 26 28]./256};
%         line_styles = {'-', '-'};
%     case 'DifOutdO'
%         labels   = {'EW-HW', 'EL-HL'};
%         [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
%         diff_pairs = {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'HdWn'))]...
%             [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]};
%         colors = {[31 120 180]./256, [227 26 28]./256};
%         line_styles = {'-', '-'};
%     case 'DifOutWL'
%         labels   = {'EW-EL', 'HW-HL'};
%         [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
%         diff_pairs = {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'EzLs'))]...
%             [find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'HdLs'))]};
%         colors = {[31 120 180]./256, [227 26 28]./256};
%         line_styles = {'-', '-'};
%     case 'DifOutS'
%         labels   = {'HW-EW', 'EL-HL', 'Su'};
%         [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
%         diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
%             [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]...
%             [find(strcmp(cond_lab,'Su')) 0]};
%         colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
%         line_styles = {'-', '-', '-'};
%     case 'DifOutSx'
%         labels   = {'HW-EW', 'EL-HL', 'Su-EWHL'};
%         [cond_lab,~,~,~,~] = fn_condition_label_styles(factor_name);
%         diff_pairs = {[find(strcmp(cond_lab,'HdWn')) find(strcmp(cond_lab,'EzWn'))]...
%             [find(strcmp(cond_lab,'EzLs')) find(strcmp(cond_lab,'HdLs'))]...
%             {find(strcmp(cond_lab,'Su')) {[find(strcmp(cond_lab,'EzWn')) find(strcmp(cond_lab,'HdLs'))]}}};
%         colors = {[31 120 180]./256, [227 26 28]./256, [51 160 44]./256};
%         line_styles = {'-', '-', '-'};

end