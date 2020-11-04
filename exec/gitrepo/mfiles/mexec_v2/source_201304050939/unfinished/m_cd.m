% not function
% script to change directory
% assume target directory is in variable 'dd'

mstar_list = DIRLIST(:,1);
target_list = DIRLIST(:,2);

k = strmatch(dd,mstar_list);
if isempty(k)
    disp('target not found')
    return
end
target = target_list{k};

fulldir = [DIR_ROOT '/' target];

cd(fulldir);