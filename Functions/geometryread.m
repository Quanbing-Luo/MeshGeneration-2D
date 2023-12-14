function  [nodes, elements, edges] = geometryread(filename)

arguments (Input)
    filename  {mustBeFile}
end

arguments (Output)
    nodes  {mustBeNumeric,mustBeRow(nodes,2)}
    elements{mustBeInteger,mustBePositive,mustBeRow(elements,3)}
    edges  {mustBeInteger,mustBePositive,mustBeRow(edges,2)}
end

fileIDgeo = fopen(filename,'r');
sznds = fread(fileIDgeo, 1, 'uint64');
nodes = fread(fileIDgeo, [2,sznds], 'double');

szels = fread(fileIDgeo, 1, 'uint64');
elements = fread(fileIDgeo, [3,szels], 'uint64')+1;

szegs = fread(fileIDgeo, 1, 'uint64');
edges = fread(fileIDgeo, [2,szegs], 'uint64')+1;

fclose(fileIDgeo);

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

