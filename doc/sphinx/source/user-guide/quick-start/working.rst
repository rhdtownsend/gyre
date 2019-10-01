<>

--------------

= Prerequisites =

To install GYRE you’ll first need to ensure that the following packages
are available on your system:

-  The [[http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk|MESA
   SDK]] software development kit plus its own prerequisites.

= Downloading & Unpacking the Source Code =

[[../downloads|Download]] the tarball containing the latest version of
GYRE, and unpack it using the ##tar## command. The example here will
unpack version 5.1 of the code into your home directory:

{{{ cd ~/ tar xvf gyre-5.1.tar.gz }}}

This places everything inside a new directory, ##~~/gyre##. Of course,
it’s possible to unpack the code anywhere you want; all that matters is
you remember where it is. A helpful way of doing this (and in fact, a
//mandatory step// in version 4.0 and later) is to store the full path
to the GYRE directory in the GYRE_DIR environment variable; for example,

{{{ export GYRE_DIR=~/gyre }}}

if you’re using the Bourne shell, or

{{{ setenv GYRE_DIR ~/gyre }}}

if you’re using the C shell. To avoid having to do this each time after
logging in, these commands can be added to the
##\ [STRIKEOUT:/.profile## (Bourne shell) and/or ##]/.cshrc## (C shell)
files.

With the environment variable set in this way, the top-level GYRE
directory can be referenced in shell commands as ##$GYRE_DIR##. The
examples in the following discussion (and elsewhere on this Wiki) make
extensive use of this shorthand.

= Updating the Source Code =

If at a later time you want to update GYRE to a more-recent release,
it’s simplest to delete everything from ##$GYRE_DIR## downwards and then
repeat the steps above.

= Compiling the Source Code =

To compile the source code, change into ##$GYRE_DIR## and build using
the ##make## tool:

{{{ cd $GYRE_DIR make }}}

To speed up the build process, you may wish to add the ##-j## flag to
##make##. Although some warning messages may be printed out during
compilation, these can safely be ignored. A successful build will
produce the following executables in the ##$GYRE_DIR/bin## subdirectory:

-  ##gyre## : the main adiabatic/non-adiabatic oscillation code
-  ##build_poly## : a support code for creating polytropic models
-  ##poly_to_fgong## : a support code for converting polytropic models
   to FGONG format

= Testing the Code =

To check that the main oscillation code has compiled correctly and give
reasonable results, change into ##$GYRE_DIR## and run tests using the
##make## tool:

{{{ cd $GYRE_DIR make test }}}

The output from the tests should look something like this:

{{{ TEST numerics (OpenMP)… …succeeded TEST numerics (band matrix)…
…succeeded }}}

If any of the tests fail, then you can submit a bug report in the
[[http://www.astro.wisc.edu/~townsend/gyre-forums/viewforum.php?f=3|GYRE
forums]]. Otherwise, if things look good, then it’s time to proceed to
the [[Running GYRE]] page.
