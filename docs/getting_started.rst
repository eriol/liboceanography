Getting Started
===============

Requirements
------------

In order to compile ``liboceanography`` you need:
    * a C compiler;
    * `cmake <http://www.cmake.org/>`_ 2.8 or greater.

For compiling unit tests you also need `Check unit test framework for C
<http://check.sourceforge.net/>`_.

Download
--------

You can download a source tarball of ``liboceanography`` from::

    http://downloads.mornie.org/liboceanography

or you can get latest code on the Mercurial repository::

    http://hg.mornie.org/oceanography/liboceanography

Installing liboceanography
--------------------------

``liboceanography`` has a cmake based configuration system that handles
build, test and installation.

The following instructions are for **unix**  type systems.

Once downloaded the tarball you have to do::

    $ tar xzf liboceanography-1.0.0.tar.gz
    $ cd liboceanography-1.0.0/
    $ mkdir build
    $ cmake ..
    $ make
    $ make test
    $ make install

