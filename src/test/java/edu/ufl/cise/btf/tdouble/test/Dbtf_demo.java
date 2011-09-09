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

import edu.ufl.cise.btf.tdouble.io.MatrixInfo;
import edu.ufl.cise.btf.tdouble.io.MatrixSize;
import edu.ufl.cise.btf.tdouble.io.MatrixVectorReader;

import junit.framework.TestCase;

public class Dbtf_demo extends TestCase {

	static final String WEST_0479 = "west0479.mtx";

	public void btf_demo_test() {

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

			int nnz = size.numEntries () ;
			int [] i = new int [nnz] ;
			int [] p = new int [nnz] ;
			double [] x ;

			x = new double [nnz] ;
			reader.readCoordinate (i, p, x) ;

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
	}

}
