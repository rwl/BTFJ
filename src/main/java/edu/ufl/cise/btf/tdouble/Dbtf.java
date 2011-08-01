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

public class Dbtf {

	/* ====================================================================== */
	/* === BTF marking of singular columns ================================== */
	/* ====================================================================== */

	/* BTF_FLIP is a "negation about -1", and is used to mark an integer j
	 * that is normally non-negative.  BTF_FLIP (-1) is -1.  BTF_FLIP of
	 * a number > -1 is negative, and BTF_FLIP of a number < -1 is positive.
	 * BTF_FLIP (BTF_FLIP (j)) = j for all integers j.  UNFLIP (j) acts
	 * like an "absolute value" operation, and is always >= -1.  You can test
	 * whether or not an integer j is "flipped" with the BTF_ISFLIPPED (j)
	 * macro.
	 */

	protected static int BTF_FLIP(int j)
	{
		return (-(j)-2) ;
	}

	protected static boolean BTF_ISFLIPPED(int j)
	{
		return ((j) < -1) ;
	}

	protected static int BTF_UNFLIP(int j)
	{
		return ((BTF_ISFLIPPED (j)) ? BTF_FLIP (j) : (j)) ;
	}

	/* ====================================================================== */
	/* === BTF version ====================================================== */
	/* ====================================================================== */

	/* All versions of BTF include these definitions.
	 * As an example, to test if the version you are using is 1.2 or later:
	 *
	 *      if (BTF_VERSION >= BTF_VERSION_CODE (1,2)) ...
	 *
	 * This also works during compile-time:
	 *
	 *      #if (BTF >= BTF_VERSION_CODE (1,2))
	 *          printf ("This is version 1.2 or later\n") ;
	 *      #else
	 *          printf ("This is an early version\n") ;
	 *      #endif
	 */

	protected static final String BTF_DATE = "Jan 25, 2011" ;
	protected static int BTF_VERSION_CODE(int main, int sub)
	{
		return ((main) * 1000 + (sub));
	}
	protected static final int BTF_MAIN_VERSION = 1 ;
	protected static final int BTF_SUB_VERSION = 1 ;
	protected static final int BTF_SUBSUB_VERSION = 2 ;
	protected static final int BTF_VERSION = BTF_VERSION_CODE(BTF_MAIN_VERSION,
			BTF_SUB_VERSION) ;

}
