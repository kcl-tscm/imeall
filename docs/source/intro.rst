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
Alternatively if you don't wish to do any development a docker image 
is available from `Docker Image <https://github.com/kcl-tscm/imeall-docker>`_.



