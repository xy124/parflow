dnl Define a macro for supporting AMPS

AC_DEFUN([CASC_SUPPORT_AMPS],[

# Begin CASC_SUPPORT_AMPS
# Determine which AMPS layer to support.
# Defines AMPS and AMPS_SPLIT_FILE
AC_ARG_WITH(amps,
[ --with-amps=AMPS_TYPE  Set the version of AMPS to use: seq, mpi1, smpi, win32],
, with_amps=seq)

case "$with_amps" in
  no)
    AMPS=seq
  ;;
  seq)
    AMPS=seq
  ;;
  mpi1)
    AMPS=mpi1
  ;;
  smpi)
    AMPS=smpi
  ;;
  win32)
    AMPS=win32
  ;;
  *)
    AC_MSG_ERROR([Invalid AMPS version specified, use seq, mpi1, smpi, win32])
  ;;    
esac
AC_MSG_RESULT([configuring AMPS $AMPS support])
AC_SUBST(AMPS)
AC_DEFINE_UNQUOTED(AMPS,$AMPS,AMPS porting layer)

AC_ARG_WITH(amps,
[ --with-amps-sequential-io  Use AMPS sequentail I/O],
  AC_DEFINE(AMPS_SPLIT_FILE),
)

# END CASC_SUPPORT_AMPS

])dnl End definition of CASC_SUPPORT_AMPS
