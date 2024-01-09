function  meshradiiwrite(nodes, elements, radii, filename)
arguments (Input)
    nodes  {mustBeNumeric,mustBeRow(nodes,2)}
    elements{mustBeInteger,mustBePositive,mustBeRow(elements,3)}
    radii  {mustBeNumeric,mustBePositive,mustBeRow(radii,1)}
    filename  {mustBeTextScalar, mustBeNonzeroLengthText}
end

fileIDgeo = fopen(filename,'w');
fwrite(fileIDgeo, size(nodes,2),'uint64');
fwrite(fileIDgeo, [nodes; radii],'double');

fwrite(fileIDgeo, size(elements,2),'uint64');
fwrite(fileIDgeo, elements - 1,'uint64');
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

