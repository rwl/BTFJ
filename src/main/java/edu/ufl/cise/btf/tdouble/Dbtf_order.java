/**
 * BTF, by Timothy A. Davis, Copyright (C) 2004-2011, University of Florida,
 * with support from Sandia National Laboratories.  All Rights Reserved.
 * Copyright (C) 2011 Richard Lincoln
 *
 * BTF is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This Module is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this Module; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *
 */

package edu.ufl.cise.btf.tdouble;

/**
 * Find a permutation P and Q to permute a square sparse matrix into upper block
 * triangular form.  A(P,Q) will contain a zero-free diagonal if A has
 * structural full-rank.  Otherwise, the number of nonzeros on the diagonal of
 * A(P,Q) will be maximized, and will equal the structural rank of A.
 *
 * Q[k] will be "flipped" if a zero-free diagonal was not found.  Q[k] will be
 * negative, and j = BTF_UNFLIP (Q [k]) gives the corresponding permutation.
 *
 * R defines the block boundaries of A(P,Q).  The kth block consists of rows
 * and columns R[k] to R[k+1]-1.
 *
 * If maxwork > 0 on input, then the work performed in btf_maxtrans is limited
 * to maxwork*nnz(A) (excluding the "cheap match" phase, which can take another
 * nnz(A) work).  On output, the work parameter gives the actual work performed,
 * or -1 if the limit was reached.  In the latter case, the diagonal of A(P,Q)
 * might not be zero-free, and the number of nonzeros on the diagonal of A(P,Q)
 * might not be equal to the structural rank.
 *
 * See btf.h for more details.
 *
 * Copyright (c) 2004-2007.  Tim Davis, University of Florida,
 * with support from Sandia National Laboratories.  All Rights Reserved.
 */
public class Dbtf_order extends Dbtf_internal {

	/**
	 * BTF_ORDER permutes a square matrix into upper block triangular form.  It
	 * does this by first finding a maximum matching (or perhaps a limited matching
	 * if the work is limited), via the btf_maxtrans function.  If a complete
	 * matching is not found, BTF_ORDER completes the permutation, but flags the
	 * columns of P*A*Q to denote which columns are not matched.  If the matrix is
	 * structurally rank deficient, some of the entries on the diagonal of the
	 * permuted matrix will be zero.  BTF_ORDER then calls btf_strongcomp to find
	 * the strongly-connected components.
	 *
	 * On output, P and Q are the row and column permutations, where i = P[k] if
	 * row i of A is the kth row of P*A*Q, and j = BTF_UNFLIP(Q[k]) if column j of
	 * A is the kth column of P*A*Q.  If Q[k] < 0, then the (k,k)th entry in P*A*Q
	 * is structurally zero.
	 *
	 * The vector R gives the block boundaries, where block b is in rows/columns
	 * R[b] to R[b+1]-1 of the permuted matrix, and where b ranges from 1 to the
	 * number of strongly connected components found.
	 *
	 * This function only operates on square matrices (either structurally full-
	 * rank, or structurally rank deficient).
	 *
	 * @param n A is n-by-n in compressed column form
	 * @param Ap size n+1
	 * @param Ai size nz = Ap [n]
	 * @param maxwork do at most maxwork*nnz(A) work in the maximum
	 * transversal; no limit if <= 0
	 * @param work work performed in maxtrans, or -1 if limit reached
	 * @param P size n, row permutation
	 * @param Q size n, column permutation
	 * @param R size n+1.  block b is in rows/cols R[b] ... R[b+1]-1
	 * @param nmatch # nonzeros on diagonal of P*A*Q
	 * @param Work size 5n
	 * @return number of blocks found
	 */
	public static int btf_order(int n, int[] Ap, int[] Ai, double maxwork,
		    double work, int[] P, int[] Q, int[] R, int nmatch, int[] Work)
	{
		int[] Flag ;
	    int nblocks, i, j, nbadcol ;

	    /* ------------------------------------------------------------------ */
	    /* compute the maximum matching */
	    /* ------------------------------------------------------------------ */

	    /* if maxwork > 0, then a maximum matching might not be found */

	    nmatch = Dbtf_maxtrans.btf_maxtrans (n, n, Ap, Ai, maxwork, work, Q,
	    		Work) ;

	    /* ------------------------------------------------------------------ */
	    /* complete permutation if the matrix is structurally singular */
	    /* ------------------------------------------------------------------ */

	    /* Since the matrix is square, ensure BTF_UNFLIP(Q[0..n-1]) is a
	     * permutation of the columns of A so that A has as many nonzeros on the
	     * diagonal as possible.
	     */

	    if (nmatch < n)
	    {
	        /* get a size-n work array */
	        Flag = Work + n ;
	        for (j = 0 ; j < n ; j++)
	        {
	            Flag [j] = 0 ;
	        }

	        /* flag all matched columns */
	        for (i = 0 ; i < n ; i++)
	        {
	            j = Q [i] ;
	            if (j != EMPTY)
	            {
	                /* row i and column j are matched to each other */
	                Flag [j] = 1 ;
	            }
	        }

	        /* make a list of all unmatched columns, in Work [0..nbadcol-1]  */
	        nbadcol = 0 ;
	        for (j = n-1 ; j >= 0 ; j--)
	        {
	            if (!(Flag [j] != 0))
	            {
	                /* j is matched to nobody */
	                Work [nbadcol++] = j ;
	            }
	        }
	        ASSERT (nmatch + nbadcol == n) ;

	        /* make an assignment for each unmatched row */
	        for (i = 0 ; i < n ; i++)
	        {
	            if (Q [i] == EMPTY && nbadcol > 0)
	            {
	                /* get an unmatched column j */
	                j = Work [--nbadcol] ;
	                /* assign j to row i and flag the entry by "flipping" it */
	                Q [i] = BTF_FLIP (j) ;
	            }
	        }
	    }

	    /* The permutation of a square matrix can be recovered as follows: Row i is
	     * matched with column j, where j = BTF_UNFLIP (Q [i]) and where j
	     * will always be in the valid range 0 to n-1.  The entry A(i,j) is zero
	     * if BTF_ISFLIPPED (Q [i]) is true, and nonzero otherwise.  nmatch
	     * is the number of entries in the Q array that are non-negative. */

	    /* ------------------------------------------------------------------ */
	    /* find the strongly connected components */
	    /* ------------------------------------------------------------------ */

	    nblocks = Dbtf_strongcomp.btf_strongcomp (n, Ap, Ai, Q, P, R, Work) ;
	    return (nblocks) ;
	}

}
