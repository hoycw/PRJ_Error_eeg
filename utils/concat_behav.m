function bhv_mod = concat_behav(bhv, bhv_oddball);
fn_bhv = fieldnames(bhv)';
fn_bhv_oddball = fieldnames(bhv_oddball)';
for index = 1:(length(fn_bhv))
    a = find(strcmp(fn_bhv{index}, fn_bhv_oddball));
    if a
        bhv_field = bhv.(fn_bhv{index});
        bhv_field_oddball = bhv_oddball.(fn_bhv_oddball{a});
        bhv_mod.(fn_bhv{index}) = vertcat(bhv_field_oddball, bhv_field);
    end
end