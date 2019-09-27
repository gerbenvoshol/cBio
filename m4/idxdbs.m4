AC_DEFUN([CHECK_IDXDBS],
[
AC_MSG_CHECKING(for CBIO pre-indexed databases)


if test -f ./cbio/index/edam.xac; then
AC_MSG_RESULT(yes)
else
AC_MSG_RESULT(no)
echo ""
echo "Pre-indexed edam, taxon + drcat databases not found."
echo "Please download them from gerbenvoshol.nl, install them"
echo "and then repeat the configure step."
exit 1
fi
])
