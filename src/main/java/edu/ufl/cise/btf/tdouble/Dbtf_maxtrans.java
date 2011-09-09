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
 * Finds a column permutation that maximizes the number of entries on the
 * diagonal of a sparse matrix.  See btf.h for more information.
 *
 * This function is identical to cs_maxtrans in CSparse, with the following
 * exceptions:
 *
 *  (1) cs_maxtrans finds both jmatch and imatch, where jmatch [i] = j and
 *      imatch [j] = i if row i is matched to column j.  This function returns
 *      just jmatch (the Match array).  The MATLAB interface to cs_maxtrans
 *      (the single-output cs_dmperm) returns imatch, not jmatch to the MATLAB
 *      caller.
 *
 *  (2) cs_maxtrans includes a pre-pass that counts the number of non-empty
 *      rows and columns (m2 and n2, respectively), and computes the matching
 *      using the transpose of A if m2 < n2.  cs_maxtrans also returns quickly
 *      if the diagonal of the matrix is already zero-free.  This pre-pass
 *      allows cs_maxtrans to be much faster than maxtrans, if the use of the
 *      transpose is warranted.
 *
 *      However, for square structurally non-singular matrices with one or more
 *      zeros on the diagonal, the pre-pass is a waste of time, and for these
 *      matrices, maxtrans can be twice as fast as cs_maxtrans.  Since the
 *      maxtrans function is intended primarily for square matrices that are
 *      typically structurally nonsingular, the pre-pass is not included here.
 *      If this maxtrans function is used on a matrix with many more columns
 *      than rows, consider passing the transpose to this function, or use
 *      cs_maxtrans instead.
 *
 *  (3) cs_maxtrans can operate as a randomized algorithm, to help avoid
 *      rare cases of excessive run-time.
 *
 *  (4) this maxtrans function includes an option that limits the total work
 *      performed.  If this limit is reached, the maximum transveral might not
 *      be found.
 *
 * Thus, for general usage, cs_maxtrans is preferred.  For square matrices that
 * are typically structurally non-singular, maxtrans is preferred.  A partial
 * maxtrans can still be very useful when solving a sparse linear system.
 *
 * Copyright (c) 2004-2007.  Tim Davis, University of Florida,
 * with support from Sandia National Laboratories.  All Rights Reserved.
 */
public class Dbtf_maxtrans extends Dbtf_internal {

	/**
	 * Perform a depth-first-search starting at column k, to find an augmenting
	 * path.  An augmenting path is a sequence of row/column pairs (i1,k), (i2,j1),
	 * (i3,j2), ..., (i(s+1), js), such that all of the following properties hold:
	 *
	 *      * column k is not matched to any row
	 *      * entries in the path are nonzero
	 *      * the pairs (i1,j1), (i2,j2), (i3,j3) ..., (is,js) have been
	 *          previously matched to each other
	 *      * (i(s+1), js) is nonzero, and row i(s+1) is not matched to any column
	 *
	 * Once this path is found, the matching can be changed to the set of pairs
	 * path.  An augmenting path is a sequence of row/column pairs
	 *
	 *      (i1,k), (i2,j1), (i3,j2), ..., (i(s+1), js)
	 *
	 * Once a row is matched with a column it remains matched with some column, but
	 * not necessarily the column it was first matched with.
	 *
	 * In the worst case, this function can examine every nonzero in A.  Since it
	 * is called n times by maxtrans, the total time of maxtrans can be as high as
	 * O(n*nnz(A)).  To limit this work, pass a value of maxwork > 0.  Then at
	 * most O((maxwork+1)*nnz(A)) work will be performed; the maximum matching might
	 * not be found, however.
	 *
	 * This routine is very similar to the dfs routine in klu_kernel.c, in the
	 * KLU sparse LU factorization package.  It is essentially identical to the
	 * cs_augment routine in CSparse, and its recursive version (augment function
	 * in cs_maxtransr_mex.c), except that this routine allows for the search to be
	 * terminated early if too much work is being performed.
	 *
	 * The algorithm is based on the paper "On Algorithms for obtaining a maximum
	 * transversal" by Iain Duff, ACM Trans. Mathematical Software, vol 7, no. 1,
	 * pp. 315-330, and "Algorithm 575: Permutations for a zero-free diagonal",
	 * same issue, pp. 387-390.  The code here is a new implementation of that
	 * algorithm, with different data structures and control flow.  After writing
	 * this code, I carefully compared my algorithm with MC21A/B (ACM Algorithm 575)
	 * Some of the comparisons are partial because I didn't dig deeply into all of
	 * the details of MC21A/B, such as how the stack is maintained.  The following
	 * arguments are essentially identical between this code and MC21A:
	 *
	 * maxtrans     MC21A,B
	 * --------     -------
	 * n            N           identical
	 * k            JORD        identical
	 * Ap           IP          column / row pointers
	 * Ai           ICN         row / column indices
	 * Ap[n]        LICN        length of index array (# of nonzeros in A)
	 * Match        IPERM       output column / row permutation
	 * nmatch       NUMNZ       # of nonzeros on diagonal of permuted matrix
	 * Flag         CV          mark a node as visited by the depth-first-search
	 *
	 * The following are different, but analogous:
	 *
	 * Cheap        ARP         indicates what part of the a column / row has
	 *                          already been matched.
	 *
	 * The following arguments are very different:
	 *
	 * -            LENR        # of entries in each row/column (unused in maxtrans)
	 * Pstack       OUT         Pstack keeps track of where we are in the depth-
	 *                          first-search scan of column j.  I think that OUT
	 *                          plays a similar role in MC21B, but I'm unsure.
	 * Istack       PR          keeps track of the rows in the path.  PR is a link
	 *                          list, though, whereas Istack is a stack.  Maxtrans
	 *                          does not use any link lists.
	 * Jstack       OUT? PR?    the stack for nodes in the path (unsure)
	 *
	 * The following control structures are roughly comparable:
	 *
	 * maxtrans                     MC21B
	 * --------                     -----
	 * for (k = 0 ; k < n ; k++)    DO 100 JORD=1,N
	 * while (head >= 0)            DO 70 K=1,JORD
	 * for (p = Cheap [j] ; ...)    DO 20 II=IN1,IN2
	 * for (p = head ; ...)         DO 90 K=1,JORD
	 *
	 * @param k which stage of the main loop we're in
	 * @param Ap column pointers, size n+1
	 * @param Ai row indices, size nz = Ap [n]
	 * @param Match size n,  Match [i] = j if col j matched to i
	 * @param Cheap rows Ai [Ap [j] .. Cheap [j]-1] alread matched
	 * @param Flag Flag [j] = k if j already visited this stage
	 * @param Istack size n.  Row index stack.
	 * @param Jstack size n.  Column index stack.
	 * @param Pstack size n.  Keeps track of position in adjacency list
	 * @param work work performed by the depth-first-search
	 * @param maxwork maximum work allowed
	 * @return
	 */
	public static int augment(int k, int[] Ap, int[] Ai, int[] Match,
			int[] Cheap, int[] Flag, int[] Istack, int[] Jstack, int[] Pstack,
		    double work, double maxwork)
	{
		/* local variables, but "global" to all DFS levels: */
		int found ; /* true if match found.  */
		int head ;  /* top of stack */

		/* variables that are purely local to any one DFS level: */
		int j2 ;    /* the next DFS goes to node j2 */
		int pend ;  /* one past the end of the adjacency list for node j */
		int pstart ;
		int quick ;

		/* variables that need to be pushed then popped from the stack: */
		int i ;     /* the row tentatively matched to i if DFS successful */
		int j ;     /* the DFS is at the current node j */
		int p ;     /* current index into the adj. list for node j */
		/* the variables i, j, and p are stacked in Istack, Jstack, and Pstack */

		quick = (maxwork > 0) ? 1 : 0 ;

		/* start a DFS to find a match for column k */
		found = FALSE ;
		i = EMPTY ;
		head = 0 ;
		Jstack [0] = k ;
		ASSERT (Flag [k] != k) ;

		while (head >= 0)
		{
			j = Jstack [head] ;
			pend = Ap [j+1] ;

			if (Flag [j] != k)          /* a node is not yet visited */
			{

				/* ---------------------------------------------------------- */
				/* prework for node j */
				/* ---------------------------------------------------------- */

				/* first time that j has been visited */
				Flag [j] = k ;
				/* cheap assignment: find the next unmatched row in col j.  This
				 * loop takes at most O(nnz(A)) time for the sum total of all
				 * calls to augment. */
				for (p = Cheap [j] ; p < pend && !(found != 0); p++)
				{
					i = Ai [p] ;
					found = (Match [i] == EMPTY) ? 1 : 0 ;
				}
				Cheap [j] = p ;

				/* ---------------------------------------------------------- */

				/* prepare for DFS */
				if (found != 0)
				{
					/* end of augmenting path, column j matched with row i */
					Istack [head] = i ;
					break ;
				}
				/* set Pstack [head] to the first entry in column j to scan */
				Pstack [head] = Ap [j] ;
			}

			/* -------------------------------------------------------------- */
			/* quick return if too much work done */
			/* -------------------------------------------------------------- */

			if (quick != 0 && work > maxwork)
			{
				/* too much work has been performed; abort the search */
				return (EMPTY) ;
			}

			/* -------------------------------------------------------------- */
			/* DFS for nodes adjacent to j */
			/* -------------------------------------------------------------- */

			/* If cheap assignment not made, continue the depth-first search.  All
			 * rows in column j are already matched.  Add the adjacent nodes to the
			 * stack by iterating through until finding another non-visited node.
			 *
			 * It is the following loop that can force maxtrans to take
			 * O(n*nnz(A)) time. */

			pstart = Pstack [head] ;
			for (p = pstart ; p < pend ; p++)
			{
				i = Ai [p] ;
				j2 = Match [i] ;
				ASSERT (j2 != EMPTY) ;
				if (Flag [j2] != k)
				{
					/* Node j2 is not yet visited, start a depth-first search on
					 * node j2.  Keep track of where we left off in the scan of adj
					 * list of node j so we can restart j where we left off. */
					Pstack [head] = p + 1 ;
					/* Push j2 onto the stack and immediately break so we can
					 * recurse on node j2.  Also keep track of row i which (if this
					 * search for an augmenting path works) will be matched with the
					 * current node j. */
					Istack [head] = i ;
					Jstack [++head] = j2 ;
					break ;
				}
			}

			/* -------------------------------------------------------------- */
			/* determine how much work was just performed */
			/* -------------------------------------------------------------- */

			work += (p - pstart + 1) ;

			/* -------------------------------------------------------------- */
			/* node j is done, but the postwork is postponed - see below */
			/* -------------------------------------------------------------- */

			if (p == pend)
			{
				/* If all adjacent nodes of j are already visited, pop j from
				 * stack and continue.  We failed to find a match. */
				head-- ;
			}
		}

		/* postwork for all nodes j in the stack */
		/* unwind the path and make the corresponding matches */
		if (found != 0)
		{
			for (p = head ; p >= 0 ; p--)
			{
				j = Jstack [p] ;
				i = Istack [p] ;

				/* ---------------------------------------------------------- */
				/* postwork for node j */
				/* ---------------------------------------------------------- */
				/* if found, match row i with column j */
				Match [i] = j ;
			}
		}
		return (found) ;
	}

	/**
	 * BTF_MAXTRANS: finds a permutation of the columns of a matrix so that it has a
	 * zero-free diagonal.  The input is an m-by-n sparse matrix in compressed
	 * column form.  The array Ap of size n+1 gives the starting and ending
	 * positions of the columns in the array Ai.  Ap[0] must be zero. The array Ai
	 * contains the row indices of the nonzeros of the matrix A, and is of size
	 * Ap[n].  The row indices of column j are located in Ai[Ap[j] ... Ap[j+1]-1].
	 * Row indices must be in the range 0 to m-1.  Duplicate entries may be present
	 * in any given column.  The input matrix  is not checked for validity (row
	 * indices out of the range 0 to m-1 will lead to an undeterminate result -
	 * possibly a core dump, for example).  Row indices in any given column need
	 * not be in sorted order.  However, if they are sorted and the matrix already
	 * has a zero-free diagonal, then the identity permutation is returned.
	 *
	 * The output of btf_maxtrans is an array Match of size n.  If row i is matched
	 * with column j, then A(i,j) is nonzero, and then Match[i] = j.  If the matrix
	 * is structurally nonsingular, all entries in the Match array are unique, and
	 * Match can be viewed as a column permutation if A is square.  That is, column
	 * k of the original matrix becomes column Match[k] of the permuted matrix.  In
	 * MATLAB, this can be expressed as (for non-structurally singular matrices):
	 *
	 *      Match = maxtrans (A) ;
	 *      B = A (:, Match) ;
	 *
	 * except of course here the A matrix and Match vector are all 0-based (rows
	 * and columns in the range 0 to n-1), not 1-based (rows/cols in range 1 to n).
	 * The MATLAB dmperm routine returns a row permutation.  See the maxtrans
	 * mexFunction for more details.
	 *
	 * If row i is not matched to any column, then Match[i] is == -1.  The
	 * btf_maxtrans routine returns the number of nonzeros on diagonal of the
	 * permuted matrix.
	 *
	 * In the MATLAB mexFunction interface to btf_maxtrans, 1 is added to the Match
	 * array to obtain a 1-based permutation.  Thus, in MATLAB where A is m-by-n:
	 *
	 *      q = maxtrans (A) ;      % has entries in the range 0:n
	 *      q                       % a column permutation (only if sprank(A)==n)
	 *      B = A (:, q) ;          % permuted matrix (only if sprank(A)==n)
	 *      sum (q > 0) ;           % same as "sprank (A)"
	 *
	 * This behaviour differs from p = dmperm (A) in MATLAB, which returns the
	 * matching as p(j)=i if row i and column j are matched, and p(j)=0 if column j
	 * is unmatched.
	 *
	 *      p = dmperm (A) ;        % has entries in the range 0:m
	 *      p                       % a row permutation (only if sprank(A)==m)
	 *      B = A (p, :) ;          % permuted matrix (only if sprank(A)==m)
	 *      sum (p > 0) ;           % definition of sprank (A)
	 *
	 * This algorithm is based on the paper "On Algorithms for obtaining a maximum
	 * transversal" by Iain Duff, ACM Trans. Mathematical Software, vol 7, no. 1,
	 * pp. 315-330, and "Algorithm 575: Permutations for a zero-free diagonal",
	 * same issue, pp. 387-390.  Algorithm 575 is MC21A in the Harwell Subroutine
	 * Library.  This code is not merely a translation of the Fortran code into C.
	 * It is a completely new implementation of the basic underlying method (depth
	 * first search over a subgraph with nodes corresponding to columns matched so
	 * far, and cheap matching).  This code was written with minimal observation of
	 * the MC21A/B code itself.  See comments below for a comparison between the
	 * maxtrans and MC21A/B codes.
	 *
	 * This routine operates on a column-form matrix and produces a column
	 * permutation.  MC21A uses a row-form matrix and produces a row permutation.
	 * The difference is merely one of convention in the comments and interpretation
	 * of the inputs and outputs.  If you want a row permutation, simply pass a
	 * compressed-row sparse matrix to this routine and you will get a row
	 * permutation (just like MC21A).  Similarly, you can pass a column-oriented
	 * matrix to MC21A and it will happily return a column permutation.
	 *
	 * @param nrow A is nrow-by-ncol in compressed column form
	 * @param ncol
	 * @param Ap size ncol+1
	 * @param Ai size nz = Ap [ncol]
	 * @param maxwork do at most maxwork*nnz(A) work; no limit if <= 0. This
	 * work limit excludes the O(nnz(A)) cheap-match phase.
	 * @param work work = -1 if maxwork > 0 and the total work performed
	 * reached the maximum of maxwork*nnz(A)).
	 * Otherwise, work = the total work performed.
	 * @param Match size nrow.  Match [i] = j if column j matched to row i
	 * @param Work size 5*ncol
	 * @return # of columns in the matching
	 */
	public static int btf_maxtrans(int nrow, int ncol, int[] Ap, int[] Ai,
		    double maxwork, double work, int[] Match, int[] Work)
	{
		int[] Cheap, Flag, Istack, Jstack, Pstack ;
		int i, j, k, nmatch, work_limit_reached, result ;

		/* ------------------------------------------------------------------ */
		/* get workspace and initialize */
		/* ------------------------------------------------------------------ */

		//Cheap  = Work ; Work += ncol ;
		Cheap = new int [ncol] ;
		//Flag   = Work ; Work += ncol ;
		Flag = new int [ncol] ;

		/* stack for non-recursive depth-first search in augment function */
		//Istack = Work ; Work += ncol ;
		Istack = new int [ncol] ;
		//Jstack = Work ; Work += ncol ;
		Jstack = new int [ncol] ;
		//Pstack = Work ;
		Pstack = new int [ncol] ;

		/* in column j, rows Ai [Ap [j] .. Cheap [j]-1] are known to be matched */
		for (j = 0 ; j < ncol ; j++)
		{
			Cheap [j] = Ap [j] ;
			Flag [j] = EMPTY ;
		}

		/* all rows and columns are currently unmatched */
		for (i = 0 ; i < nrow ; i++)
		{
			Match [i] = EMPTY ;
		}

		if (maxwork > 0)
		{
			maxwork *= Ap [ncol] ;
		}
		work = 0 ;

		/* ------------------------------------------------------------------ */
		/* find a matching row for each column k */
		/* ------------------------------------------------------------------ */

		nmatch = 0 ;
		work_limit_reached = FALSE ;
		for (k = 0 ; k < ncol ; k++)
		{
			/* find an augmenting path to match some row i to column k */
			result = augment (k, Ap, Ai, Match, Cheap, Flag, Istack, Jstack, Pstack,
					work, maxwork) ;
			if (result == TRUE)
			{
				/* we found it.  Match [i] = k for some row i has been done. */
				nmatch++ ;
			}
			else if (result == EMPTY)
			{
				/* augment gave up because of too much work, and no match found */
				work_limit_reached = TRUE ;
			}
		}

		/* ------------------------------------------------------------------ */
		/* return the Match, and the # of matches made */
		/* ------------------------------------------------------------------ */

		/* At this point, row i is matched to j = Match [i] if j >= 0.  i is an
		 * unmatched row if Match [i] == EMPTY. */

		if (work_limit_reached != 0)
		{
			/* return -1 if the work limit of maxwork*nnz(A) was reached */
			work = EMPTY ;
		}

		return (nmatch) ;
	}

}
