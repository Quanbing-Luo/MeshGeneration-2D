function meshgeneration(geofile, meshfile, NameValueArgs)

arguments (Input)    
    geofile  {mustBeFile}
    meshfile  {mustBeTextScalar, mustBeNonzeroLengthText}
    NameValueArgs.Hmax {mustBePositive, mustBeScalarOrEmpty} = realmax
    NameValueArgs.Hgrad {mustBeFloat, mustBeScalarOrEmpty, ...      
        mustBeInRange(NameValueArgs.Hgrad,1,2, 'exclude-lower')} = 1.1
    NameValueArgs.Hvertex { mustBeRequiredCellData} 
    NameValueArgs.Hedge { mustBeRequiredCellData} 
end

currentfolder=cd;
geofilefull =  ['"' fullfile(currentfolder,geofile) '"'];
meshfilefull =  ['"' fullfile(currentfolder,meshfile) '"'];
meshFeaturefilefullchar =  fullfile(currentfolder,'meshFeature.txt') ;
meshFeaturefilefull =  ['"' meshFeaturefilefullchar '"'];

fileID = fopen(meshFeaturefilefullchar,'w');
fprintf(fileID, 'Hmax: %g\n', NameValueArgs.Hmax);
fprintf(fileID, 'Hgrad: %g\n', NameValueArgs.Hgrad);

if isfield(NameValueArgs,"Hvertex")
    lenv=length(NameValueArgs.Hvertex);
    nv=0;
    for i=1:2:lenv
        nv =nv+length(NameValueArgs.Hvertex{i});
    end
    fprintf(fileID, 'Hvertex: %u\n', nv);
    ii=0;
    for i=1:2:lenv
        first =NameValueArgs.Hvertex{i}-1;
        second =NameValueArgs.Hvertex{i+1};
        for j= 1:length(first)
            fprintf(fileID, '%u %u %g\n', ii, first(j), second);
            ii=ii+1;
        end
    end
else
    fprintf(fileID, 'Hvertex: 0\n');
end


if isfield(NameValueArgs,"Hedge")
    lene=length(NameValueArgs.Hedge);
    ne=0;
    for i=1:2:lene
        ne =ne+length(NameValueArgs.Hedge{i});
    end
    fprintf(fileID, 'Hedge: %u\n', ne);
    ii=0;
    for i=1:2:lene
        first =NameValueArgs.Hedge{i}-1;
        second =NameValueArgs.Hedge{i+1};
        for j= 1:length(first)
            fprintf(fileID, '%u %u %g\n', ii, first(j), second);
            ii=ii+1;
        end
    end
else
    fprintf(fileID, 'Hedge: 0\n');
end

fclose(fileID);


[funpath,name] = fileparts(mfilename('fullpath'));
command0= ['cd "' funpath '"'];
command1= ['.\' name ' ' geofilefull ' ' meshFeaturefilefull ' ' meshfilefull];

% status0 = system(command0);
% if status0 
%     error("Fail to set current run path.");
% end

status1 = system([command0 '&&' command1]);
if status1 
    error("Fail to run the command of meshgeneration.exe.");
end

end


function mustBeRequiredCellData(a)

if ~iscell(a) ||  ~isvector(a)  || rem(length(a) , 2) ~= 0
    eidType = 'meshgeneration:notEvenCellVector';
    msgType = 'Input must be a cell vector containing an even number of elements.';
    throwAsCaller(MException(eidType,msgType))
end

for i=1:2:length(a)
    if  ~isvector(a{i})  || any(a{i}<=0)  ||  any(a{i} ~= floor(a{i}))
        eidType = 'meshgeneration:notInteger';
        msgType = 'Odd value of input must be positive integer vector.';
        throwAsCaller(MException(eidType,msgType))
    end
end

for i=2:2:length(a)
    if  ~isscalar(a{i}) || any(a{i}<=0) 
        eidType = 'meshgeneration:notScalar';
        msgType = 'Even value of input must be a positive scalar.';
        throwAsCaller(MException(eidType,msgType))
    end
end

end



