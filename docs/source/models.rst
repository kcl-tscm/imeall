Models and Queries
==================

These modules define the data models for :class:`GrainBoundary` and :class:`SubGrainBoundary` 
structures and the SQL and directory tree database.
They also contains classes for performing maintenance of the 
grain boundary database and querying and copying structures from the directory tree.

Maintenance
------------
:class:`imeall.models.GBMaintenance` contains methods for performing 
maintenance and query operations on the directory tree structure and 
modifying `gb.json` and `subgb.json` files.

.. autoclass:: imeall.models.GBMaintenance
  :members:

Queries
--------
:class:`imeall.models.GBQuery` contains methods for checking out and
copying structure files and json dictionaries from the 
directory tree grain boundary database.

.. autoclass:: imeall.models.GBQuery
  :members:

SQL Models
----------
These classes define the data models for the :class:`imeall.gb_models.GrainBoundary`
and :class:`imeall.gb_models.SubGrainBoundary models`.

.. autoclass:: imeall.gb_models.GrainBoundary
  :members:

.. autoclass:: imeall.gb_models.SubGrainBoundary
  :members:

.. autoclass:: imeall.gb_models.Fracture
  :members:

Serialization Methods
----------------------
These methods flip vector or list data back and forth
between comma separated value strings which can be stored
in the SQL database.

.. automodule:: imeall.gb_models
  :members: serialize_vector, deserialize_vector_float, deserialize_vector_int, create_tables

SQL Database Methods
---------------------
Methods for updating and checking the integrity of the SQL database. These methods ensure
the directory tree and the SQL database remain normalized and synchronized.

.. automodule:: imeall.gb_models
  :members: add_conv_key, gb_check_dir_integrity, gb_check_path, gb_check_conv, gb_check_force, change_json_key, insert_subgrain, populate_db

