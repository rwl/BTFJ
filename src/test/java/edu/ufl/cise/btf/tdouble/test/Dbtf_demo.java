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

package edu.ufl.cise.btf.tdouble.test;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;

import static edu.ufl.cise.btf.tdouble.Dbtf_order.btf_order;
import static edu.ufl.cise.btf.tdouble.Dbtf_strongcomp.btf_strongcomp;
import static edu.ufl.cise.btf.tdouble.Dbtf_maxtrans.btf_maxtrans;

import edu.ufl.cise.btf.tdouble.io.MatrixInfo;
import edu.ufl.cise.btf.tdouble.io.MatrixSize;
import edu.ufl.cise.btf.tdouble.io.MatrixVectorReader;

import junit.framework.TestCase;

public class Dbtf_demo extends TestCase {

	static final String WEST_0479 = "west0479.mtx";

	int n;
	int [] Ai, Ap;

	protected void setUp() throws Exception {
		int nnz;
		double [] x ;
		FileReader fileReader;
		MatrixVectorReader reader;
		MatrixInfo info;
		MatrixSize size;

		try
		{
			fileReader = new FileReader (WEST_0479) ;
			reader = new MatrixVectorReader (fileReader) ;

			info = reader.readMatrixInfo () ;
			size = reader.readMatrixSize (info) ;

			assertTrue("btf: A must be sparse, square, and non-empty",
					size.isSquare() ) ;

			n = size.numRows () ;
			nnz = size.numEntries () ;

			/* get sparse matrix A */
			Ai = new int [nnz] ;
			Ap = new int [nnz] ;
			x = new double [nnz] ;
			reader.readCoordinate (Ai, Ap, x) ;

		}
		catch (FileNotFoundException e)
		{
			fail () ;
			e.printStackTrace () ;
		}
		catch (IOException e)
		{
			fail () ;
			e.printStackTrace () ;
		}

		super.setUp();
	}

	public void btf_test() {
		int maxwork, nblocks ;
		int [] Q, P, R, Work ;
		MutableDouble work = new MutableDouble () ;
		MutableInt nmatch = new MutableInt () ;

		/* get output arrays */
		P = new int [n] ;
		Q = new int [n] ;
		R = new int [n+1] ;

		/* get workspace */
		Work = new int [5*n] ;

		maxwork = 0 ;
		work.setValue( 0 ) ;

		/* find the permutation to BTF */
		nblocks = btf_order (n, Ap, Ai, maxwork, work, P, Q, R, nmatch, Work) ;

		// TODO BTF assertions
	}

	public void strongcomp_test() {

		int nblocks ;
		int [] Q, P, R, Work ;

		/* get output arrays */
		P = new int [n] ;
		Q = null ;
		R = new int [n+1] ;

		/* get workspace of size 4n */
		Work = new int [4*n] ;

		/* find the strongly-connected components of A */
		nblocks = btf_strongcomp (n, Ap, Ai, Q, P, R, Work) ;

		// TODO strongcomp assertions
	}

	public void maxtrans_test() {

		int maxwork, nmatch ;
		int [] Work, Match ;
		MutableDouble work = new MutableDouble () ;

		/* get output array */
		Match = new int [n] ;

		/* get workspace of size 4n */
		Work = new int [4*n] ;

		maxwork = 0 ;
		work.setValue( 0 ) ;

		/* perform the maximum transversal */
		nmatch = btf_maxtrans (n, n, Ap, Ai, maxwork, work, Match, Work) ;

		// TODO maxtrans assertions

	}

}
