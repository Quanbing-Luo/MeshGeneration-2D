function  geometrywrite(nodes, elements, edges, filename)
arguments (Input)
    nodes  {mustBeNumeric,mustBeRow(nodes,2)}
    elements{mustBeInteger,mustBePositive,mustBeRow(elements,3)}
    edges  {mustBeInteger,mustBePositive,mustBeRow(edges,2)}
    filename  {mustBeTextScalar, mustBeNonzeroLengthText}
end

fileIDgeo = fopen(filename,'w');
fwrite(fileIDgeo, size(nodes,2),'uint64');
fwrite(fileIDgeo, nodes,'double');

fwrite(fileIDgeo, size(elements,2),'uint64');
fwrite(fileIDgeo, elements - 1,'uint64');

fwrite(fileIDgeo, size(edges,2),'uint64');
fwrite(fileIDgeo, edges - 1,'uint64');
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

