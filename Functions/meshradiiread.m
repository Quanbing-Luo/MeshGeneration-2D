function  [nodes, elements, radii] = meshradiiread(filename)

arguments (Input)
    filename  {mustBeFile}
end

arguments (Output)
    nodes  {mustBeNumeric,mustBeRow(nodes,2)}
    elements{mustBeInteger,mustBePositive,mustBeRow(elements,3)}
    radii  {mustBeNumeric,mustBeRow(radii,1)}
end

fileID = fopen(filename,'r');

sznds = fread(fileID, 1, 'uint64');
array = fread(fileID, [3,sznds], 'double');
nodes=array(1:2,:); radii=array(3,:); 

szels = fread(fileID, 1, 'uint64');
elements = fread(fileID, [3,szels], 'uint64')+1;

fclose(fileID);

end


% Custom validation function
function mustBeRow(mx,n)
    % Test for equal size
    if size(mx,1)~=n
        eid = 'georead:wrongRow';
        msg = ['Array must have' char(n) 'row.'];
        throwAsCaller(MException(eid,msg))
    end
end

