# Cofactor promiscuity (i.e., when a protein can perform the same function with different cofactors or combinations thereof)
- Examples: proteins that accept a "Cu cation" (Cu2+, Cu3+, etc.), a "divalent metal cation" (e.g., Zn2+, Fe2+, etc.), or that bind "2 divalent metal cations per subunit. Magnesium or manganese."
- Possible solutions:
  - Adding irreversible reactions converting each cofactor into a pseudometabolite representing the cofactor (e.g., "mg2_or_mn2_c" formed from either mg2_c or mn2_c).

NEW nomenclature:
PTM-rt2761+2mn2
PTM-rt2761+2mg2
PTM-rt2761+1mg2+1mn2

Old PTM options:
mg2ormn2_c:2

New PTM options:
(mg2 OR mn2):2
equivalent to "mg2:2 OR (mg2:1 AND mn2:1) OR mn2:2"
"metal2" for all metal divalent cations

- This would create 3 versions of rxn.
- Why not use separate cofactor allocation rxns representing possible demands for them, more easily handling permutations? B/c this wouldn't allow formation of different species with different functions. (SINK-mn2_c-for-rt2761_c and SINK-mg2_c-for-rt2761_c) 
- Alternatively, could just make cofactor synthesis occur at enzyme formation, cutting down on bloat.