echo Win9x/WinNT
echo off
rem JH Zhao 6/6/2001
rem TMPDIR is necessary for WinNT
rem otherwise pointr still needs /tmp (after env ... g77)
set TMPDIR=.
del poipro
del poiter
nucfama poidat poijob poipro poiter
move INTERMED.1 INTERMEDFILE
pointr  poijob poipro poiter
del INTERMED*
set TMPDIR=
echo done
