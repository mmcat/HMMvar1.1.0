AC_INIT([hmmvar],[1.0],[mingml@vt.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
		 Makefile
		 src/Makefile
		 examples/Makefile
		 ])

# check psiblast
AC_PATH_PROG([PSIBLAST],[psiblast],[null])
if test "$PSIBLAST" = "null"
then
	AC_MSG_ERROR(NCBI psiblast missing from path. Please specify full path. (e.g. ./configure PSIBLAST=/path/to/psiblast/exe))
fi

# check muscle
AC_PATH_PROG([MUSCLE],[muscle],[null])
if test "$MUSCLE" = "null"
then
	AC_MSG_ERROR(Multialigner mulscel missing from path. Please specify full path.(e.g. ./configure MUSCLE=/path/to/muscle/exe) )
fi

#check hmmbuild/hmmsearch
AC_PATH_PROG([HMMER],[hmmbuild hmmsearch],[null])
if test "$HMMER" = "null"
then
	AC_MSG_ERROR(HMMER missing from path. Please specify full path. (e.g. ./configure HMM=/path/to/hmm/exe))
fi

#check blastdbcmd
AC_PATH_PROG([BLASTDBCMD],[blastdbcmd],[null])
if test "$BLASTDBCMD" = "null"
then
	AC_MSG_ERROR(blastdbcmd missing from path. Please specify full path. (e.g. ./configure BLASTDBCMD=/path/to/blastdbcmd/exe))
fi

AC_SUBST([BLAST_DB])

AC_OUTPUT


