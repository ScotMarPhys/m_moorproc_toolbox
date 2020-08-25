% routine just to list the available data report tools
% run without any arguments
function data_report_tools
    directory=which('data_report_tools');
    index=strfind(directory,'/');
    directory=directory(1:index(end));
    ls([directory '*.m'])
end