#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
It implements the ssGSEA method (Subramanian et al. PNAS 2005).

Author: Daniel Nicorici, Daniel.Nicorici@gmail.com
Copyright (c) 2023 Daniel Nicorici

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""





def ssgsea(x, x_rows_genes, geneset, alpha = 0.25, scale = True, norm = False, min_genes=15, max_genes=500):
    """
    x = matrix as list of lists, ie. list of rows (rows = genes, columns = samples)
    x_rows_genes = labels for rows of x that are the gene names (eg. genes)
    geneset = one set of genes
    """
    genes = set(geneset).intersection(x_rows_genes)

    n = len(genes)
    if n < min_genes or n > max_genes or (not x):
        score = 0
    else:
        
        if type(x[0]) != list:
            # just one list as input instead list of lists => convert it to list of lists
            x = [[e] for e in x]

        # ranks (the most expressed genes get the highest rank)
        z = [[s.index(e)+1 for e in v] for v in list(zip(*x)) if (s:=sorted(v))]
        # sort decreasing the genes based on the ranks (most expressed gene ends in top of list)
        z = [ sorted(zip(x_rows,v),key=lambda k: k[1],reverse=True) for v in z]

    
        m = len(z[0]) # number of genes in the matrix
        # PG
        pg = [[e[1]**alpha if e[0] in genes else 0 for e in v] for v in z]
        png = [[1 if e==0 else 0 for e in v] for v in pg] # this is faster to be computed like this
        pg = [ list(map(lambda k,s: k/s,v,[sum(v)]*m)) for v in pg]
        png = [ list(map(lambda k,s: k/s,v,[sum(v)]*m)) for v in png]
        pg = [ (u:=v[0])+sum((u:=e+u) for e in v[1:]) for v in pg]
        png = [ (u:=v[0])+sum((u:=e+u) for e in v[1:]) for v in png]

        score = [e[0]-e[1] for e in zip(pg,png)]
        
        if scale:
            score = score / len(genes)
        
        if norm:
            score = score / (max(score)-min(score))
    
    return score

if __name__ == "__main__":
    x=[5.5,6.3,5.4,7.3,5.9,6.8,11.2,6.4,8.2,7.5]
    
    x=list(zip(x,x))
    x_rows=["Gene1","Gene2","Gene3","Gene4","Gene5","Gene6","Gene7","Gene8","Gene9","Gene10"]

    genesetA=["Gene1", "Gene3", "Gene5","Gene7"]
    genesetB=["Gene2", "Gene4", "Gene8"]

    ssgsea(x, x_rows, genesetA, scale=False, min_genes=1)
    # [-1.7476307995924376, -1.7476307995924376]

    ssgsea(x, x_rows, genesetB, scale=False, min_genes=1)
    # [-0.16603102972695538, -0.16603102972695538]
    
