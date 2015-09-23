# End-members-
calculate the number and permutations of unique end members for sublattice model used in CALPHAD method.

Case 1 All the sublattices are identical (FCC structure belongs to this case)

number of unique endmembers
recurrence formula for a(m,n)
 n>=m
a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)
 n<m
a(m,n) = C(n,1)*a(m-1,1)+...+C(n,k)*a(m-k,k)+...C(n,m)*a(0,m)

Permutations  of unique endmembers are also calculated in recurrence formula and 
coorrespond to the calculation of number.
e.g. for permutations p(m,n), the C(n,1)*a(m-1,1) term 
in the first sublattice is one of the n constituents, the rest of the m-1 sublattices
is the permutations of p(m-1,1)

Case 2 Two kinds of sublattices in symmetric position (BCC structure belongs to this case)

Suppose we have 2*m sublattices (2 kinds of sublattices in symmetric position, each 
kind have m sublattices) and n constituents in each sublattices, and the number of unique endmembers 
for this structure is b(2*m,n). The formula of general term of b(2*m,n) can be written into a 
ecurrence formula as follows:
 b(2*m,n) = a(2,a(m,n))

Suppose p(m:m,n) are the permutations of structure with 2m sublattices 
(2 kinds of sublattices in symmetric position, each kind have m sublattices) 
and n constituents in each, the strategy to obtain p(m:m,n) is as follows:
p(m:m,n) are the permutations of 2 sublattices, in each sublattices is the 
permutations of unique endmembers in m identical sublattices and n constituents in each.
