function varargout = new_dir( varargin )
%%new_dir extends MATLAB's dir, includes filtering options and executables
%the_dir = new_dir( path ) acts just like dir, but two extra fields are
%returned, which are "path", for the path of the file, and "number", which
%is the last number in the filename, and is useful for sequentially
%numbered data files (see examples below).
%[the_dir, number] = new_dir( path, search_string) returns also the numbers
%in the filenames, and search_string uses regexp rules for searching the
%filenames
%new_dir( path, search_numbers ) locates the files with numbers inside the
%specified search_numbers
%new_dir( path, ..., -function_name ) executes function_name for each
%file in the directory, examples are -display, -delete, -load, -open
%new_dir( path, ..., \-search_string ) performs search where '-' is
%included (i.e. not executed as a function)
%
%Options are 
%-first, which is used in conjunction with functions, and
%restricts the function to only operate on the first file.
%-dir, which is used to retrieve only directories
%~-dir, which is used to retrieve non-directories
%
%examples:
%let's generate a directory with files for experimenting
%new_dir( '-generate' ); new_dir test_new_dir
%to display files containing `data1' in the filename with extension '.mat',
%and having suffix numbers 12:15 in their filename
%new_dir( 'test_new_dir', '\.mat', 'data1', 12:15 );
%to load file w/number 22
%new_dir( 'test_new_dir', 22, '-load' ); i
%to delete file w/number 22
%new_dir( 'test_new_dir', 22, '-delete' ); new_dir test_new_dir \.mat
%to display files with extension .m:
%the_dir = new_dir( 'test_new_dir', '\.m$', '-display' );
%to not display files with extension .mat:
%[the_dir, number] = new_dir( 'test_new_dir/*.mat' );
%to display files with extension .mat and load the first one into the workspace:
%new_dir( 'test_new_dir/*.mat', '-load', '-first' ), i
%to find files with the fragment '.*ata11.*' in either name or path(for unix-based systems):
%[the_dir, number] = new_dir( 'test_new_dir/*.mat', '-find', '.*ata11.*' );
%to locate and open files with extension .m
%new_dir( 'test_new_dir/*.m', '-open' ); 
%  or with regexp rules
%new_dir( 'test_new_dir', '\.m$', '-open' ); 
%to display only directories
%[the_dir, number] = new_dir( 'test_new_dir', '-dir' );
%to display only not directories
%[the_dir, number] = new_dir( 'test_new_dir', '~-dir' );
%
%2009 by S. Koehler
%%
% varargin = {'/users/skoehler/MATLAB/NORA/delta/C2F6', '-dir', '~^b' }
if any( strcmp( varargin, '-generate' ) )
    mkdir test_new_dir;
    fprintf( 'generating test data files\n' );
    for i=10:30;
        disp( fullfile( 'test_new_dir', sprintf( 'data%i.mat', i ) ) );
        save( fullfile( 'test_new_dir', sprintf( 'data%i.mat', i ) ), 'i' );
    end;
    str = 'hi mom';
    save( fullfile( 'test_new_dir', sprintf( 'data%i.m', i ) ), 'str', '-ASCII' );
    str = 'hello world';
    save( fullfile( 'test_new_dir', '~example.m' ), 'str', '-ASCII' );    
    return;
end;
if any( strcmp( varargin, '-find' ) )
    %%
    if ~isunix, error('-find works on Unix only'); end
    [status, tmp] = system(  sprintf( 'find %s -regex "%s"', varargin{1}, varargin{[false strcmp( varargin, '-find' )]} ) );
    if status
        error( tmp );
    end
    j = [0 find( tmp == 10 )];
    d = [];
    the_path = {};
    for i=1:numel(j)-1
        if ~strncmp( tmp(j(i)+1:j(i+1)-1), 'find: ', 5 )
            if isempty( d)
                d = dir( tmp(j(i)+1:j(i+1)-1) );
            else
                d(size(d,1)+1,1) = dir( tmp(j(i)+1:j(i+1)-1) );
            end;
            the_path{end+1} = fileparts( tmp(j(i)+1:j(i+1)-1) );
        end
    end
    %%
    varargin( [false strcmp( varargin, '-find' )] ) = [];
    varargin( strcmp( varargin, '-find' ) ) = [];
else
    d = dir( varargin{1} );
    the_path = repmat( {get_path(varargin{1})}, size( d,1 ), 1);
end
varargin(1) = [];
if any( strcmp( varargin, '-dir' ) )
    d(~[d.isdir]) = [];
    varargin( strcmp( varargin, '-dir' ) ) = [];
end
if any( strcmp( varargin, '~-dir' ) )
    d(~[d.isdir]) = [];
    varargin( strcmp( varargin, '~-dir' ) ) = [];
end
try
    names = {d.name};
catch
    error( 'could not locate any such files' );
end
str = regexprep( names,  {'\.+[^\.]*$', '^\D*', '.*\D(\d+)\D*$'}, {'nan', 'nan', '$1 '} );
numbers = sscanf( strvcat(str)', '%f' );
[dummy, ind] = sort( numbers);
search_numbers = [varargin{ cellfun( 'isclass', varargin, 'double' ) }];
if ~isempty( search_numbers )
    ind =ind(ismember( numbers(ind), search_numbers ));
    varargin( cellfun( 'isclass', varargin, 'double' ) ) = [];
end;
if any( strcmp( varargin, '-first' ) )
    first_flag = true;
    varargin(strcmp( varargin, '-first' ) ) = [];
else
    first_flag = false;
end

to_do = strncmp( varargin, '-', 1 );
pos_queery = ~to_do&~strncmp(varargin, '~', 1 );
neg_queery = strncmp(varargin, '~', 1 );
varargin = regexprep( varargin, {'^\\-', '-', '^~', '^\\~', '\\'}, {'', '', '', '~', '\'} );
for str = varargin(pos_queery)
	ind = ind( ~cellfun( 'isempty', regexp( names(ind), str{1} ) ) );
end;
for str = varargin(neg_queery)
	ind = setdiff( ind, ind( ~cellfun( 'isempty', regexp( names(ind), str{1} ) ) ) );
end;

if isempty( ind )
    varargout{1} = struct( 'name', {}, 'date', {}, 'bytes', {}, 'datenum', {}, 'number', {}, 'path', {} );
    return;
end;
d = struct( 'name', {d(ind).name}', 'date', {d(ind).date}', 'bytes', {d(ind).bytes}', 'datenum', {d(ind).datenum}', 'number', num2cell(numbers(ind)), ...
    'path', shiftdim(the_path(ind) ) );
for str = varargin(to_do)
    for i=1:numel(d)
        if ~strcmp( str{1}, 'display' )
            fprintf( '%s %s\n', str{1}, fullfile(  d(i).path, d(i).name ) );
        end;
        try
            evalin( 'caller', sprintf( '%s %s', str{1}, fullfile( d(i).path, d(i).name ) ) );
        catch
            warning( 'failed: %s %s', str{1}, fullfile( d(i).path, d(i).name ) );
        end
        if any( strcmp( varargin, 'first' ) )
            break
        end
    end
end;
%%
switch nargout
    case 0
        disp( shiftdim({d.name} ));
    case 1
        varargout = {d};
    case 2
        varargout = {d, numbers(ind)};
end
return;

function path = get_path(varargin1)
orig_dir = pwd;
if isdir( varargin1 )
    cd( varargin1 )
elseif isdir( fileparts(varargin1 ) );
    cd( fileparts(varargin1 ) );
end
path = pwd;
cd( orig_dir );
return;