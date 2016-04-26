function file_list=get_file_list(directory)
%file_list=get_file_list(directory)
%returns a filelist cell{}, if no directory is given, '.' is used.
if nargin<1 
  directory='.';
end

f_struct=dir(directory);
file_list={};
j=1;
for i=1:length(f_struct)
  if ~strcmp(f_struct(i).name, '.') & ~strcmp(f_struct(i).name, ...
					      '..') & ~strcmp(f_struct(i).name, '.svn')
    file_list{j}=f_struct(i).name;
    j=j+1;
  end %if
end %for i
end
