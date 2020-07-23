function ps_rename(ps_new)
    % rename a today's date .ps file to a new input name.
ps = ['spm_' datestr(datetime('now'),'yyyymmmdd') '.ps'];
ps_root = extractBefore(ps_new,'.ps');
if exist(ps,'file')
    if exist(ps_new,'file')
        ps_num = dir([ps_root '*']);
        ps_new = [ps_root '_run' num2str(ps_num) '.ps'];
    end
    movefile(ps,ps_new)
else
    fprintf(2,'No spm .ps file found in working directory\n')
end