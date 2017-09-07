Installation
============

Manual Installation
--------------------

The following dependencies can be installed using `pip <https://pypi.python.org/pypi/pip>`_.

|  `ASE <https://wiki.fysik.dtu.dk/ase/>`_
|  `Flask <http://flask.pocoo.org/>`_
|  `peewee <http://docs.peewee-orm.com/en/latest/index.html>`_

Follow the installation instructions for the `QUIP <https://libatoms.github.io/QUIP/install.html>`_ package.

To use the imeall package to perform calculations with the mixed Fe-H EAM potential you must checkout the appropriate branch of QUIP.
  ``git checkout mod_eam``

The above should be sufficient to run the imeall package. 

To install the Imeall package clone the github repository:
		``git clone https://github.com/Montmorency/imeall.git``

Add the location of the ``imeall`` directory to your ``PYTHONPATH`` variable:
  ``export PYTHONPATH=$HOME/imeall::$PYTHONPATH``

Installation from Docker Image
------------------------------
Alternatively if you don't wish to do any development a 
`Docker Image <https://github.com/kcl-tscm/imeall-docker>`_ is available.


Environment
-----------
  In order to use imeall effectively the following modifications should be made.
The configuration file for the imeall package is located at ~/$HOME/imeall/instance/config.py.
The standard location for the imeall grain boundary directory tree database is in the
``imeall`` module (i.e. the directory with ``__init__.py``). If your local directory 
tree is in a non-standard location modify the appropriate configuration variable:
  ``GRAIN_DATABASE= /path/to/grain_database/``

A `representative subset <https://imeall.co.uk>`_ of the grain boundary database to download for 
testing is available from the ``imeall`` website hosted on the `NOMAD <https://nomad-coe.eu>`_ server.

As previously mentioned the imeall module location should be appended to your ``PYTHONPATH`` variable.

