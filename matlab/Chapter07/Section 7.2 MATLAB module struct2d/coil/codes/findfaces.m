function varargout = findfaces(simplicial_complex, i)
%   SYNTAX
%   IFACES = findfaces(SIMPLICIAL_COMPLEX, I);
%   [IFACES NSIMPLEX_INDICES] = findfaces(SIMPLICIAL_COMPLEX, I);
%   DESCRIPTION
%   FINDFACES(SIMPLICIAL_COMPLEX, I) returns a list of the I-faces
%   of the N-simplices in the (N+1)-by-M matrix SIMPLICIAL_COMPLEX.  An
%   N-simplex is the N-dimensional generalization of a line segment,
%   triangle, or tetrahedron, which are the 1-, 2-, and 3-simplex,
%   respectively.  An I-face of an N-simplex is an I-simplex that forms
%   part of the boundary of the N-simplex.  A vertex is considered a
%   0-face, so for purposes of this program, a point is considered a
%   0-simplex.
%
%   IFACES = findfaces(SIMPLICIAL_COMPLEX, I) returns IFACES with each
%   I-face appearing exactly once.
%
%   [IFACES NSIMPLEX_INDICES] = findfaces(SIMPLICIAL_COMPLEX,I) returns
%   IFACES with each face appearing as many times in IFACES as there are
%   N-simplices in SIMPLICIAL_COMPLEX that it bounds.
%
% Inputs:
%   SIMPLICIAL_COMPLEX, an array of simplices, represents each N-simplex as
%   one of its M columns. Each simplex consists of N+1 distinct positive
%   integers, each of which represents a vertex of the simplex.  Different
%   numbers are treated as representing distinct points, and vice versa.
%
%   I specifies the number of dimensions the faces have.  It may be a
%   nonnegative integer or a string; if it is the latter, the case of the
%   characters is ignored.  Allowed values of the string, with their
%   corresponding integers, are as follows:
%       Integer         Strings
%          0        'vertex' 'vertexes' 'vertices'
%          1         'edge'  'edges'
%         N-2       'ridge'  'ridges'
%         N-1       'facet'  'facets'
%   When N=3, 'face' and 'faces' are equivalent to 2.  Otherwise, those
%   strings may not be used because they are too similar to 'facet' and
%   'facets', respectively.
%
% Outputs:
%   IFACES is an (I+1)-by-(M*NCHOOSEK(N+1,I+1)) array of the faces, since
%   each of the M simplices has NCHOOSEK(N+1,I+1) I-faces.  IFACES
%   represents each face in a way similar to the way SIMPLICIAL_COMPLEX
%   represents simplices: each face consists of I+1 distinct positive
%   integers, each of which represents a distinct point which acts as one
%   vertex of the face.  If SIMPLICIAL_COMPLEX is in canonical form (i.e.,
%   each column is sorted internally, and the columns are sorted as a
%   group), then IFACES is in canonical form.
%
%   NSIMPLEX_INDICES is a row vector with as many columns as IFACES; each
%   element in NSIMPLEX_INDICES corresponding to an occurrence of an I-face
%   is a unique column index of SIMPLICIAL_COMPLEX where that face was
%   found.
%
%   EXAMPLE:  
%   Finding the edges in a rectangle consisting of two triangles.
%   If T = [1 1                                    (1) +-+ (2)
%           2 3                                        |\|
%           3 4], which represents the rectangle   (4) +-+ (3)
%
%   then findfaces(T,'edges') is [1 1 1 2 3
%                                 2 3 4 3 4]
%
%   and [E Te] = findfaces(T,1) results in 
%       E  = [1 1 1 1 2 3
%             2 3 3 4 3 4]
%
%       Te = [1 1 2 2 1 2];
%
% Terminology references:
%   Eric W. Weisstein. "Facet." From MathWorld--A Wolfram Web Resource.
%       http://mathworld.wolfram.com/Facet.html  
%   Eric W. Weisstein. "Ridge." From MathWorld--A Wolfram Web Resource.
%       http://mathworld.wolfram.com/Ridge.html 
%   Eric W. Weisstein. "Simplex." From MathWorld--A Wolfram Web Resource.
%       http://mathworld.wolfram.com/Simplex.html
%   Eric W. Weisstein. "Simplicial Complex." From MathWorld--A Wolfram Web
%       Resource. http://mathworld.wolfram.com/SimplicialComplex.html
%   Eric W. Weisstein. "Simplicial Subcomplex." From MathWorld--A Wolfram
%       Web Resource. http://mathworld.wolfram.com/SimplicialSubcomplex.html
%
%   Authors: A. Marut, S.N. Makarov
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

% Check validity of input, and make i an integer:
if nargout > 2
    error 'Too many output arguments.';
end

if ndims(simplicial_complex) > 2
    error 'FINDFACES: SIMPLICIAL_COMPLEX must be a 2-D array!';
end
[n m] = size(simplicial_complex);
n = n-1;

% if support for more strings is added, update the documentation for I
if isstr(i)
    if any(strcmpi(i, {'vertex', 'vertexes', 'vertices'}))
        i = 0;
    elseif any(strcmpi(i, {'edge', 'edges'}))
        i = 1;
    elseif n==3 && any(strcmpi(i, {'face', 'faces'}))
        i = 2;
    elseif any(strcmpi(i, {'ridge', 'ridges'}))
        i = n-2;
    elseif any(strcmpi(i, {'facet', 'facets'}))
        i = n-1;
    else
        error('FINDFACES: Unrecognized string parameter passed in I!');
    end
else
    if i ~= floor(i) || i < 0 || n < i
        error 'FINDFACES: I must be a nonnegative integer not greater than N!';
    end
end

faces_per_simplex = nchoosek(n+1,i+1);
range_to_store    = 1:faces_per_simplex;
iFaces            = zeros(i+1,m*faces_per_simplex);
selection_matrix  = nchoosek(1:n+1, i+1)';
if nargout==2
    nSimplex_indices  = zeros(1,m*faces_per_simplex);
end

% sort the columns of selection_matrix internally
selection_matrix = sort(selection_matrix);
% don't sort the columns as a group because the enumerated faces won't
% automatically be in canonical form.

for simplex_number = 1:m
    % get each simplex
    simplex = simplicial_complex(:,simplex_number);
    if i == 0
        % Workaround for the fact that column_vector(row_vector) = column_vector
        % simplex, a column vector, is one simplex.
        % simplex', a row vector, is n+1 points.
        current_faces = simplex';
    else
        % all possible combinations of i+1 vertices
        current_faces = simplex(selection_matrix);
    end
    iFaces(:,range_to_store) = current_faces;

    % these faces belong to the simplex_number-th simplex.
    if nargout==2
        nSimplex_indices(range_to_store) = simplex_number;
    end
    range_to_store = range_to_store + faces_per_simplex;
end

% The columns of iFaces are sorted internally if the input was so.
% Therefore, just sort them as a group.
if nargout == 2
    [iFaces I] = sortrows(iFaces');
    iFaces = iFaces';
    varargout{2} = nSimplex_indices(I);
else
    iFaces = unique(iFaces', 'rows')';
end

varargout{1} = iFaces;