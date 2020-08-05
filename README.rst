==================
potentials-library
==================

Introduction
------------

This repository provides a copy of the records on potentials.nist.gov for use
with the potentials Python package.  This is meant as a convenience for users
who wish to work with all potentials/crystals/etc hosted in the database.

Notes:

- Downloading from this repository is faster than downloading records through
  the potentials package.  As a git repository, subsequent updates are also
  quicker.
  
- The contained records represent a snapshop of potentials.nist.gov at a
  moment in time.  There is no guarantee that it will be 100% up to date
  (yet).

- Users are free to add their own records, clone the repository, etc.  Pull
  requests can be used to submit new records for adding to the NIST database.

Usage
-----

The directory for the potentials-library repository should correspond to the
library directory set for the potentials package:

    import potentials
    settings = potentials.Settings()
    
    # See the current library_directory
    print(settings.library_directory)
    
    # Change the library directory (optional)
    settings.set_library_directory(<PATH>)

Once the records are in the library directory, the potentials.Database class
should be able to find and load them.