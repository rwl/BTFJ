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

public class Dbtf_internal extends Dbtf {

	/**
	 * Enable debugging and assertions.
	 */
	public static boolean NDEBUG = true ;

	protected static void ASSERT (boolean a)
	{
		if (!NDEBUG)
		{
			assert a ;
		}
	}

	protected static void ASSERT (int a)
	{
		ASSERT (a != 0) ;
	}

	/**
	 * Enable diagnostic printing.
	 */
	public static boolean NPRINT = true ;

	protected static void PRINTF (String format, Object... args)
	{
		if (!NPRINT)
		{
			System.out.printf (format, args) ;
		}
	}

	protected static final int TRUE = 1 ;
	protected static final int FALSE = 0 ;
	protected static final int EMPTY = (-1) ;

	protected static int MIN (int a, int b)
	{
		return (((a) < (b)) ?  (a) : (b)) ;
	}

}
