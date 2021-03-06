function [BC_branchends, BC_nodecurrents, BC_nodeequality ] = setNodes(numsections);

%
% This file specifies the tree structure.
%
% INPUT:
%
% Number of sections in the tree
%
% OUTPUT:
%
% BC_branchends
% Specifies those sections (i.e. branches) which have sealed ends.
% Each row represents a section which has a sealed end.
% BC_branchends( sectionnumber, whichendissealed)
% sectionnumber = the assigned number of the section in the tree
% whichendissealed = 0 (left end), 1 (right end)
%
% BC_nodecurrents
% Specifies which branches meet at the nodes.
% Each row specifies the information for one node.
% BC_nodecurrents( NumberOfSectionsWhichMeet, a1, a2, b1, b2, c1, c3)
% In this example NumberOfSectionsWhichMeet = 3 and the sections are denoted a, b, c.
% NumberOfSectionsWhichMeet = the number of sections that join together at this node.
% a1 = the assigned section number of the section that meets at this node
% a2 = specifies which end joins the node (0: left, 1: right)
% Similarly for b1, b2, c1, c2
% Example: If a row contains the values (3,11,1,12,0,13,0),
% this indicates that 3 sections meet at this node. The first section has
% the number 11 and meets the node on its right, the other sections are 12
% and 13 and meet the node on their left.
%
% BC_nodeequality
% Specifies all the pairs of sections that meet (and hence have the same end point membrane potential).
% Every row is a pair of sections. E.g. a node where three sections meet has two such pairs.
% Note: in principle the information in this matrix can be extracted from BC_nodecurrents.
% But it is separately specified in this version of the code.
% BC_nodeequality( numberofsection1, endofsection1, numberofsection2, endofsection2)
% numberofsection1 = assigned number of the first section in the pair
% endofsection1 = which end of section is joined to the other section (0: left, 1: right)
% similarly for the second section of the pair.
% Example: If a row contains the values (13,1,14,0),
% this indicates that section 13 is joined on its right end to the left end
% of section 14.
%

if numsections==1,
    BC_branchends(1,1) = 1; BC_branchends(1,2) = 0; % section 1, sealed at left
    BC_branchends(2,1) = 1; BC_branchends(2,2) = 1; % section 2, sealed at right
    BC_nodecurrents = [];
    BC_nodeequality = [];
else

    fprintf('*** ERROR: [%s] This file supports one section only\n',mfilename);
    return
end



